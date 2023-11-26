#pragma once

#include <stdio.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifdef _OPENMP

#include <omp.h>

#define GET_TIME() (omp_get_wtime())

#else

#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC)

#endif

#if NDEBUG

#define DEBUG_PRINTF(fmt, ...)
#define DEBUG_PRINT(msg)

#else

#define DEBUG_PRINTF(fmt, ...)                                           \
    printf("[DEBUG][%s:%d] " fmt "\n", __FILE__, __LINE__, __VA_ARGS__); \
    fflush(stdout);
#define DEBUG_PRINT(msg)                                    \
    printf("[DEBUG][%s:%d] %s\n", __FILE__, __LINE__, msg); \
    fflush(stdout);

#endif

/**
 * @brief Copy the string passed in argument
 *
 * @param str [IN] the string to copy
 * @return char* a copy of the string passed in argument or NULL if memory
 * allocation for the copy failed
 */
char* copy_string(char* str);
