#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "genetic_structure.h"

#ifdef GENETIC_ALGORITHM_C

void parallel_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars );

void sequential_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars );

#else

extern void parallel_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars );

extern void sequential_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars );

#endif

#endif