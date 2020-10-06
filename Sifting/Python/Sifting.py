from typing import List, Any, Union

import cplex
import numpy as np
import time
from scipy import sparse


class CplexSifting(object):
    num_iter_not_in_basis: Union[List[int], Any]

    def __init__(self, A, b, c, ub, senses, filename):

        # Stores the original Cplex model instance
        self.__original_model = cplex.Cplex()

        # Stores the current restricted master problem (working problem)
        self.__rmp_model = cplex.Cplex()

        # Initialize callback functions
        self.__init_cols = None
        self.__add_cols = None
        self.__drop_cols = None

        # Initialize problem data
        self.__A = None
        self.__b = None
        self.__c = None
        self.__ub = None
        self.__senses = None

        # Incumbent objective value
        self.__incumbent_obj = 0.0

        # Incumbent dual solution
        self.__incumbent_dual = None

        # Iteration time
        self.__solution_time = []

        # Method for solving working problem, default is 0: cplex default selection
        self.__subproblem_method = 0

        # A sufficiently large number << cplex.infinity
        self.__big_M = 1e+10

        # Total number of working problems solved so far (also implies number of sifting iterations)
        self.rmp_iter_count = 0

        # Other arrays related to online algorithm
        self.__is_row_equal = None
        self.__online_sol = None
        self.__online_dual = None

        # Sifting Parameters

        self.__max_iter = 30

        # Online algorithm parameters
        self.__boosting_param = 3

        # RMP(Working problem size) parameters
        self.__max_rmp_size = 500000
        self.__num_column_in = 250000
        self.__num_column_out = 10000
        self.__allow_barrier = True

        # Column pricing paramerers
        self.__STANDARD_PRICE = 0
        self.__LAMBDA_PRICE = 1

        # Iteration status parameters
        self.__COLUMN_STATUS_NOT_IN_RMP = -2
        self.__COLUMN_STATUS_AT_UPPER_BOUND = -3

        self.__pricing_method = self.__STANDARD_PRICE
        
        # Set online algorithm implementation
        self.__online_alg = self.__sparse_online_LP

        # Read data from ".mps" format files
        # This method is very slow in practice and it is suggested that one should
        # directly use matrix information for initialization

        if filename:
            print("Reading original model from {0} ...".format(filename))
            self.__original_model.read(filename)

            print("Extracting model information: constraint matrix A")

            # This part has to be rewritten since the SparsePair returned by Cplex cannot get directly used

            Ac = []
            Ar = []
            Aval = []

            # Store problem size temporarily
            m = self.__original_model.linear_constraints.get_num()
            n = self.__original_model.variables.get_num()

            for i in range(m):
                rowLinExp = self.__original_model.linear_constraints.get_rows(i).unpack()
                idx = rowLinExp[0]
                vals = rowLinExp[1]
                nnz = len(idx)
                Ac = Ac + idx
                Ar = Ar + [i] * nnz
                Aval = Aval + vals
                if i % 100 == 0:
                    print("Extracted {0} out of {1} rows ...".format(i, m))

            A = sparse.coo_matrix((Aval, (Ar, Ac)), shape=(m, n))
            A = A.tocsc()

            print("Extracting model information: right hand side vector b")
            self.__b = np.asarray(self.__original_model.linear_constraints.get_rhs())

            print("Extracting model information: sense information")
            self.__senses = self.__original_model.linear_constraints.get_senses()

            print("Extracting model information: linear objective coefficient")
            self.__c = np.asarray(self.__original_model.objective.get_linear())

            print("Extracting model information: variable upperbounds")
            self.__ub = np.asarray(self.__original_model.variables.get_upper_bounds())

        # Directly initialize problem with data passed
        if not filename:
            self.__A = A
            self.__b = b.squeeze()
            self.__c = c.squeeze()
            self.__ub = ub.squeeze()
            self.__senses = senses

            self.__A = self.__A.tocsc()
            Abeg = self.__A.indptr
            Aind = self.__A.indices
            Aval = self.__A.data
            Acnt = Abeg[1:] - Abeg[:-1]
            Abeg = Abeg[:-1]

            m = A.shape[0]
            n = A.shape[1]

            print("Setting model from input data")

            self.__original_model.copylp(numcols=n,
                                         numrows=m,
                                         objsense=self.__original_model.objective.sense.minimize,
                                         obj=self.__c,
                                         rhs=self.__b,
                                         senses=self.__senses,
                                         matbeg=Abeg.tolist(),
                                         matcnt=Acnt.tolist(),
                                         matind=Aind.tolist(),
                                         matval=Aval,
                                         lb=[0.0] * n,
                                         ub=self.__ub,
                                         colnames=None,
                                         rownames=None)

            print("Model setting complete")

        # Initialize constraint matrix size
        self.nRow = A.shape[0]
        self.nCol = A.shape[1]

        # Initialize auxiliary array for determination of constraint type
        self.__is_row_equal = np.asarray([self.__senses[i] == "E" for i in range(self.nRow)])

        # Initialize basis information
        self.__current_aux_col_basis = [self.__original_model.start.status.at_lower_bound] * 2 * self.nRow
        self.__current_col_basis_status = [self.__original_model.start.status.at_lower_bound] * self.nCol
        self.__current_row_basis = [self.__original_model.start.status.at_lower_bound] * self.nRow
        self.__current_rmp_cols = []
        self.num_iter_not_in_basis = [self.__COLUMN_STATUS_NOT_IN_RMP] * self.nCol

        print("Sifting framework initialization complete")
        print("***************** Model Information **************")
        print("* Number of constraints: {0} ".format(self.nRow))
        print("* Number of variables: {0}".format(self.nCol))
        print("* Model lower bound: 0   Maximum upperbound: {0}".format(np.max(self.__ub)))
        print("**************************************************")

    @staticmethod
    def __sparse_online_LP(A, b, c, is_row_equal, K, b_sense):

        # Extract shape information
        m = A.shape[0]
        n = A.shape[1]
        is_constr_exist = not (m == 0)

        # Extract component arrays
        if is_constr_exist:
            A_p = A.indptr
            A_i = A.indices
            A_x = A.data
            d = b / n
        else:
            raise Warning("No constraint exists")

        # Initialize stepsize
        step = 1 / np.sqrt(n * K)

        # Initialize update information
        last_update = np.zeros(m)
        nUpdate = 0

        # Initialize solution array
        x = np.zeros(n)

        if b_sense == 0:
            y = np.ones(m)
        else:
            y = np.zeros(m)

        for k in range(K):
            p = np.random.permutation(n)
            for i in range(n):
                j = p[i]
                xk = 0
                if is_constr_exist:
                    data_idx_beg = A_p[j]
                    data_idx_end = A_p[j + 1]
                    row_idx = A_i[data_idx_beg:data_idx_end]
                    y[row_idx] = np.subtract(y[row_idx],
                                             step * np.multiply(d[row_idx], (nUpdate - last_update[row_idx])))
                    ineq_row_idx = row_idx[np.where(is_row_equal[row_idx] == 0)[0]]
                    y[ineq_row_idx] = np.maximum(0, y[ineq_row_idx])
                    last_update[row_idx] = nUpdate
                    xk += np.dot(A_x[data_idx_beg:data_idx_end], y[row_idx])

                xk = (c[j] >= xk)

                nUpdate += 1

                if xk and is_constr_exist:
                    data_idx_beg = A_p[j]
                    data_idx_end = A_p[j + 1]
                    row_idx = A_i[data_idx_beg:data_idx_end]
                    y[row_idx] = np.add(y[row_idx], step * np.subtract(A_x[data_idx_beg:data_idx_end], d[row_idx]))
                    ineq_row_idx = row_idx[np.where(is_row_equal[row_idx] == 0)[0]]
                    y[ineq_row_idx] = np.maximum(0, y[ineq_row_idx])
                    last_update[row_idx] = nUpdate

                x[j] += xk

        return x / K, y

    def __sub_lp_solve(self, cols, mode):

        # Get problem size
        m = self.nRow
        rmp_size = len(cols)

        # Select corresponding columns
        Asub = self.__A[:, cols]
        csub = self.__c[cols]
        ubsub = self.__ub[cols]

        # Create Cplex problem instance
        self.__rmp_model = cplex.Cplex()

        # Add auxiliary part
        extendc = np.concatenate([np.ones(2 * m) * self.__big_M, csub], axis=0)
        # extendc = np.concatenate([np.ones(2 * m) * cplex.infinity, csub], axis=0)
        extendub = np.concatenate([np.ones(2 * m) * cplex.infinity, ubsub], axis=0)

        idty = sparse.identity(m).tocsc()
        extendA = sparse.hstack([idty, -idty, Asub])

        Abeg = extendA.indptr
        Aind = extendA.indices
        Aval = extendA.data
        Acnt = Abeg[1:] - Abeg[:-1]
        Abeg = Abeg[:-1]

        ext_n_row = extendA.shape[0]
        ext_n_col = extendA.shape[1]

        assert (ext_n_row == m)
        assert (ext_n_col == rmp_size + 2 * m)

        self.__rmp_model.copylp(numcols=ext_n_col,
                                numrows=ext_n_row,
                                objsense=self.__rmp_model.objective.sense.minimize,
                                obj=extendc, rhs=self.__b,
                                senses=self.__senses,
                                matbeg=Abeg.tolist(),
                                matcnt=Acnt.tolist(),
                                matind=Aind.tolist(),
                                matval=Aval,
                                lb=[0.0] * ext_n_col,
                                ub=extendub,
                                colnames=None,
                                rownames=None)

        if mode == "L":
            extColBasis = self.__current_aux_col_basis + [self.__current_col_basis_status[i] for i in cols]
            initRowBasis = [self.__current_row_basis[i] for i in range(m)]
            self.__rmp_model.start.set_start(extColBasis, initRowBasis, [], [], [], [])

        self.__rmp_model.parameters.lpmethod.set(self.__subproblem_method)

        if not self.__allow_barrier:
            self.__rmp_model.parameters.barrier.limits.iteration.set(0)

        start_time = self.__rmp_model.get_dettime()
        self.__rmp_model.solve()
        elapsed_time = self.__rmp_model.get_dettime() - start_time

        basis_info = self.__rmp_model.solution.basis.get_basis()

        # Extract column basis information
        col_basis_new = basis_info[0]
        
        assert(len(col_basis_new) == rmp_size + 2 * m)

        # Extract row basis information
        row_basis_new = basis_info[1]
        
        assert(len(row_basis_new) == m)

        # Collect basis information

        # Collect basis information of auxiliary columns
        self.__current_aux_col_basis = col_basis_new[0: 2 * m]

        for i in range(rmp_size):
            self.__current_col_basis_status[cols[i]] = col_basis_new[i + 2 * m]
            # Update the array recording number of iterations not in basis for columns in working problem
            if col_basis_new[i + 2 * m] == self.__rmp_model.solution.basis.status.at_lower_bound:
                self.num_iter_not_in_basis[cols[i]] = max(self.num_iter_not_in_basis[cols[i]], 0) + 1
            elif col_basis_new[i + 2 * m] == self.__rmp_model.solution.basis.status.at_upper_bound:
                self.num_iter_not_in_basis[cols[i]] = self.__COLUMN_STATUS_AT_UPPER_BOUND

        for i in range(m):
            self.__current_row_basis[i] = row_basis_new[i]

        self.__incumbent_dual = np.asarray(self.__rmp_model.solution.get_dual_values()).reshape(1, -1)
        assert (self.__incumbent_dual.shape == (1, m))

        self.rmp_iter_count += 1
        print("************** RMP Solution Summary *************")
        print("* RMP Solution Count: {0}".format(self.rmp_iter_count))
        print("* RMP Solution Summary: ")
        print("* RMP Size: {0} + {1} "
              "(Auxiliary) ({2}) out of "
              "{3} columns".format(rmp_size, 2 * m, self.__rmp_model.variables.get_num(), self.nCol))
        print("* RMP Solution time: {0}".format(elapsed_time))

        self.__incumbent_obj = self.__rmp_model.solution.get_objective_value()

        if self.__incumbent_obj >= self.__big_M / 2:
            print("* RMP without auxiliary variables is infeasible")
        else:
            print("* RMP Objective: {0}".format(self.__rmp_model.solution.get_objective_value()))

        print("*************************************************")
        return

    def __default_init_cols(self, pricing, threshold):
        # Apply default online algorithm to get initial columns for sifting

        # The problem stored in the framework is to minimize c' * x and for online linear
        # program we need to invert the sign of c
        # b_sense = np.max(self.__b * (1 - self.__is_row_equal)) > 0
        
        b_sense = np.min(self.__b) > 0
        print("********** Column Initialization start ***********")
        if self.__online_alg == self.__sparse_online_LP:
            x_approx, y_approx = self.__online_alg(self.__A, self.__b, -self.__c,
                                                         self.__is_row_equal, self.__boosting_param, b_sense)
        else:
            x_approx = self.__online_alg(self.__A.indptr, self.__A.indices, self.__A.data, self.__b, -self.__c,
                                                        self.__is_row_equal, self.__boosting_param, b_sense)
            y_approx = None
        self.__online_sol = x_approx
        self.__online_dual = y_approx
        
        
        self.__current_rmp_cols = np.where(x_approx == 1)[0].squeeze().tolist()
        iter_not_in_basis = np.asarray(self.num_iter_not_in_basis)
        iter_not_in_basis[self.__current_rmp_cols] = 0
        self.num_iter_not_in_basis = iter_not_in_basis.tolist()
        '''
        self.__current_rmp_cols = np.where(pricing >= threshold)[0].squeeze().tolist()
        iter_not_in_basis = np.asarray(self.num_iter_not_in_basis)
        iter_not_in_basis[self.__current_rmp_cols] = 0
        self.num_iter_not_in_basis = iter_not_in_basis.tolist()
        '''
        # Do not forget to update num_iter_not_in_basis while adding columns
        
        print("* Summary: {0} out of {1} initial columns are selected".format(len(self.__current_rmp_cols), self.nCol))

    def __default_add_cols(self, candidate_cols, reduced_costs):

        # Compute the number of columns that can still be added without violating the tolerance
        num_cols_available = min(self.__max_rmp_size - len(self.__current_rmp_cols), self.__num_column_in)

        candidate_cols = np.asarray(candidate_cols)

        # Take the columns with negative reduced costs
        # This is an array of indices of columns in candidate_cols
        neg_rdc_cols_idx = np.where(reduced_costs < 0)[0]

        # Get the column indices in all columns with negative reduced costs
        neg_rdc_cols = candidate_cols[neg_rdc_cols_idx]

        # If we can add all columns
        if len(neg_rdc_cols) <= num_cols_available + self.__num_column_out:
            self.__current_rmp_cols += neg_rdc_cols.tolist()
        # If we cannot, we consider following pricing rules
        else:
            neg_rdcs = reduced_costs[neg_rdc_cols_idx]
            if self.__pricing_method == self.__STANDARD_PRICE:
                # Standard pricing rule
                neg_rdc_ord = np.argsort(neg_rdcs)
                neg_rdc_cols_idx_ord = neg_rdc_cols[neg_rdc_ord]
                self.__current_rmp_cols += neg_rdc_cols_idx_ord[0: num_cols_available + self.__num_column_out].tolist()
            else:
                # Lambda pricing rule
                lbd_prices = self.__c[neg_rdc_cols] / neg_rdcs
                lbd_price_ord = np.argsort(lbd_prices)
                neg_rdc_cols_idx_ord = neg_rdc_cols[lbd_price_ord]
                self.__current_rmp_cols += neg_rdc_cols_idx_ord[0: num_cols_available + self.__num_column_out].tolist()

        return

    def __default_drop_cols(self):

        # Get number of columns to drop
        num_cols_out = len(self.__current_rmp_cols) - self.__max_rmp_size

        # Do not drop columns if size of the working problem has not yet exceeded the tolerance
        if num_cols_out <= 0:
            return

        # Based on the principle of purging, we drop columns in working problem that have maximum
        # number of iterations not in basis
        rmp_cols = np.asarray(self.__current_rmp_cols)
        all_cols_iterations = np.asarray(self.num_iter_not_in_basis)
        rmp_iterations = all_cols_iterations[self.__current_rmp_cols]
        all_col_basis = np.asarray(self.__current_col_basis_status)

        # Since we never remove columns in basis and at upper bound, we exclude them
        # More in details, we take the columns with
        rmp_cols_removable_idx = np.where(rmp_iterations > 0)[0]
        rmp_cols_removable = rmp_cols[rmp_cols_removable_idx]

        # If we can remove all of these columns
        if len(rmp_cols_removable) <= num_cols_out:
            all_cols_iterations[rmp_cols_removable] = self.__COLUMN_STATUS_NOT_IN_RMP
        
            try:
                assert(sum(all_col_basis[rmp_cols_removable]) == 0)
            except:
                print(all_col_basis[rmp_cols_removable])
            all_col_basis[rmp_cols_removable] = self.__rmp_model.solution.basis.status.at_lower_bound
            self.num_iter_not_in_basis = all_cols_iterations.tolist()
            self.__current_rmp_cols = list(set(self.__current_rmp_cols).difference(set(rmp_cols_removable)))
            
        # If we cannot remove all of them, we remove the columns that are not in basis for
        # most iterations
        else:
            rmp_iterations_removable = rmp_iterations[rmp_cols_removable_idx]
            rmp_iterations_ord = np.argsort(rmp_iterations_removable)
            rmp_cols_removable_ord = rmp_cols_removable[rmp_iterations_ord]
            all_cols_iterations[rmp_cols_removable_ord[- num_cols_out - 1: -1]] = self.__COLUMN_STATUS_NOT_IN_RMP
            self.num_iter_not_in_basis = all_cols_iterations.tolist()
            all_col_basis[rmp_cols_removable_ord[- num_cols_out - 1: -1]] = \
                                                                self.__rmp_model.solution.basis.status.at_lower_bound
            self.__current_col_basis_status = all_col_basis.tolist()
            self.__current_rmp_cols = \
                list(set(self.__current_rmp_cols).difference(set(rmp_cols_removable_ord[- num_cols_out - 1: -1])))

    def __get_reduced_costs(self, cols):
        # This function returns the reduced costs of all columns with indices specified
        coefs = self.__c[cols]
        reduced_costs = coefs - (self.__incumbent_dual * self.__A[:, cols]).squeeze()
        return reduced_costs

    def ColGen(self):
        # This function integrates all the functionalities define above to do customized column generation
        # Check callbacks
        """
        if self.__init_cols:
            get_init_cols = self.__init_cols
        else:
            get_init_cols = self.__default_init_cols
            print("Using default initial column selection strategy")

        if self.__add_cols:
            add_cols = self.__add_cols
        else:
            add_cols = self.__default_add_cols
            print("Using default column adding strategy")

        if self.__drop_cols:
            drop_cols = self.__drop_cols
        else:
            drop_cols = self.__default_drop_cols
            print("Using default column dropping strategy")
        """

        print(" ***************************************************************\n")
        print(" ******************** Cplex Sifting Framework ******************\n")
        print(" ***************************************************************\n")

        # Get initial columns to start restricted master problem
        ollp_start = time.time()
        print("\n******   Initial Column Selection ******\n")
        
        pricing = None
        threshold = 0
        
        self.__default_init_cols(pricing, threshold)
        ollp_time = time.time() - ollp_start
        print("* Online algorithm takes {0} seconds".format(ollp_time))
        self.__solution_time.append(ollp_time)
        print("\n******   Initial Column Selection Done ******\n")

        # Solve initial RMS with extended columns
        iter_time = time.time()
        self.__sub_lp_solve(self.__current_rmp_cols, "I")
        # Get candidate columns and reduced costs
        candidate_cols = list(set(list(range(self.nCol))).difference(set(self.__current_rmp_cols)))
        reduced_costs = self.__get_reduced_costs(candidate_cols)

        # Begin sifting
        optimality = (np.min(reduced_costs) >= 0)
        num_iter = 1
        iter_time = time.time() - iter_time

        self.__solution_time.append(iter_time)
        print("\n* Summary: Iteration {0} takes {1} seconds\n".format(num_iter, iter_time))
        num_iter += 1

        while not optimality and num_iter <= self.__max_iter:
            iter_time = time.time()
            print("\n* Begin sifting iteration {0}".format(num_iter))
            print("* Current objective: {0}".format(self.__incumbent_obj))

            size_of_rmp = len(self.__current_rmp_cols)
            self.__default_add_cols(candidate_cols, reduced_costs)
            print("* Strategy added {0} "
                  "columns and there are {1} columns "
                  "in the working problem".format(len(self.__current_rmp_cols) - size_of_rmp,
                                                  len(self.__current_rmp_cols)))
            
            size_of_rmp = len(self.__current_rmp_cols)
            self.__default_drop_cols()
            print("* Strategy removed {0} "
                  "columns and there are {1} columns "
                  "in the working problem\n".format(size_of_rmp - len(self.__current_rmp_cols),
                                                  len(self.__current_rmp_cols)))

            self.__sub_lp_solve(self.__current_rmp_cols, "L")

            candidate_cols = list(set(list(range(self.nCol))).difference(set(self.__current_rmp_cols)))
            reduced_costs = self.__get_reduced_costs(candidate_cols)
            optimality = (np.min(reduced_costs) >= 0)

            num_iter += 1
            iter_time = time.time() - iter_time
            self.__solution_time.append(iter_time)
            print("\n* Summary: Iteration {0} takes {1} seconds".format(num_iter, iter_time))

        if num_iter < self.__max_iter:
            print("\n************* Sifting framework successfully ends *************")
            print("* Solution Summary: ")
            print("* Objective value                   : {0}".format(self.__rmp_model.solution.get_objective_value()))
            print("* Total number of sifting iterations: {0}".format(num_iter))
            print("***************************************************************")
            print("* Iteration Summary:")
            print("- Online Algorithm   : {0} seconds".format(ollp_time))
            for i in range(num_iter - 1):
                print("- Sifting iteration {0}: {1} seconds".format(i + 1, self.__solution_time[i + 1]))
            total_time = sum(self.__solution_time)
            sifting_time = total_time - ollp_time

            print("- Total solution time: {0} ({1} + {2}) seconds".format(total_time, ollp_time, sifting_time))
            print("***************************************************************")
        else:
            print("************* Iteration limit exceeded ****************")
            print("* Solution Summary ")
            print(
                "* Best Objective value                   : {0}".format(
                    self.__rmp_model.solution.get_objective_value()))
            
    def solve_original(self):
        self.__original_model.solve()

    def get_original_model(self):
        return self.__original_model

    def get_rmp_model(self):
        return self.__rmp_model

    def get_solution_time(self):
        return self.__solution_time

    def clear_basis(self):
        self.__current_aux_col_basis = []
        self.__current_col_basis_status = []
        self.__current_row_basis = []
        self.__current_rmp_cols = []
        
    def clear_info(self):
        self.__solution_time = []

    def set_subproblem_method(self, method):
        self.__subproblem_method = method

    def set_boosting_param(self, K):
        self.__boosting_param = K

    def set_barrier_switch(self, S):
        self.__allow_barrier = (S == 1)

    def set_max_sifting_iter(self, max_iter):
        self.__max_iter = max_iter

    def set_max_rmp_size(self, max_rmp_size):
        self.__max_rmp_size = max_rmp_size

    def set_num_cols_in(self, num_cols_in):
        self.__num_column_in = num_cols_in

    def set_num_cols_out(self, num_cols_out):
        self.__num_column_out = num_cols_out

    def set_pricing_method(self, method):
        if method == "standard":
            self.__pricing_method = self.__STANDARD_PRICE
        else:
            self.__pricing_method = self.__LAMBDA_PRICE
    
    def set_cython_method(self, func):
        self.__online_alg = func
        

    """
    def add_init_callback(self, func):
        self.__init_cols = func

    def add_add_col_callback(self, func):
        self.__add_cols = func

    def add_drop_col_callback(self, func):
        self.__drop_cols = func

    """