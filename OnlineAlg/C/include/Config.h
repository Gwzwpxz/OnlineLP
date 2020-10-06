#ifndef CONFIG_H
#define CONFIG_H



/* Define Solution Status */
#define OLLP_SOLUTION_STATUS_OPTIMAL 0
#define OLLP_SOLUTION_STATUS_INVALID -1
#define OLLP_RETCODE_OK 0
#define OLLP_RETCODE_INVALID -1

/* Check data type */
#ifdef  OLLP_LONG
#define OLLP_int ptrdiff_t
#define OLLP_ID "%td"

#define OLLP_FASTLP FastLP_l

#else
#define OLLP_int int
#define OLLP_ID "%d"

#define OLLP_FASTLP FastLP

#define OLLP_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define OLLP_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define OLLP_SIGN(x) (((x) > 0.0) ? 1.0 : (((x) < 0.0) ? -1.0 : 0.0))
#define OLLP_FREE(var) do {free((var)); (var) = NULL;} while (0)

#define OLLP_DATE "July 18th, 2020"
#define OLLP_MAIN_VERSION 0
#define OLLP_SUB_VERSION 1
#endif
#endif
