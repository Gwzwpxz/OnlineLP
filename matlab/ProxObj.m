function [fval] = ProxObj(a, c, pold, gamma, d, p)

fval = d'*p + max(c - a'*p, 0) + 0.5 * gamma * norm(p - pold)^2;

end % End function