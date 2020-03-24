/**
 * @file print.c
 * @brief Contains the initialization of VERBOSE_LEVEL
 */
#include "print.h"

/**
 * @brief Controls the amount of information that is outputted to stdout
*/
const char VERBOSE_LEVEL = VERBOSE;

/*
#ifdef GPU
#pragma omp declare target
int print(const char *format,...)
{

return 0;

}
#pragma omp end declare target
#endif
*/
