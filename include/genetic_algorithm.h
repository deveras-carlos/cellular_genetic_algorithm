#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#define FALSE 0
#define TRUE 1

#define MAX_POPULATION 40
#define AMT_NEIGHBORS 4
#define MAX_GEN 1000

#define PREBATI (rand()%101/100.F)
#define XALPHA  0.2
#define PMUTAC 10

#define GREEDY_START_RATE 0.3
#define MIN_ITER_LOCAL_SEARCH 50
#define MIN_ITERATIONS 100

#define MAX_WAIT_TRIALS 150

#define LOCAL_SEARCH_TYPE 2
#define RANDOM_TYPE 1
#define CROSSOVER_TYPE 0

#ifdef GENETIC_ALGORITHM_C
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <float.h>

typedef struct _chromossome_ {
    double* genes;
    unsigned int neighbors[ AMT_NEIGHBORS ];
    double fitness;
    char type; // 0 if crossover, 1 if random, 2 if local search
    char sels;
} Chromossome;

typedef struct _population_ {
        Chromossome     reading_individuals[ MAX_POPULATION ];
        Chromossome     writing_individuals[ MAX_POPULATION ];
        double          sum_fitness;
        int             population_size;
        int             individual_size;
        int             best;
        int             worst;
        int             mutation_amount;
        int             equals;
        int             best_generation;

        double          lower_limit;
        double          upper_limit;
} Population;

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit
);

int fix_unfeasible( double* xr, double lower_limit, double upper_limit );

void blend_crossover( Population* population, int father, int mother, int son, float alpha );

// void mutation( Population* population, unsigned int individual, unsigned int current_generation, int exponent );

// void local_search( Population* population, double ( *fitness_function )( double*, int n ) );

void genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit );

#else

extern void genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit );

#endif

#endif