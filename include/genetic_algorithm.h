#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#define READING 0
#define WRITING 1
#define SWITCH( i ) ( P.individuals[ i ].RW_type = !P.individuals[ i ].RW_type )

#define MAX_GENES 10
#define ROWS 6
#define COLS 4
#define MAX_POPULATION ROWS * COLS
#define AMT_NEIGHBORS 4
#define MAX_GEN 1000

#ifdef GENETIC_ALGORITHM_C
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct _chromossome_ {
    double genes[MAX_GENES];
    unsigned int neighbors[ AMT_NEIGHBORS ];
    double fitness;
    char RW_type; // Reading and Writing type - 0 for Reading, 1 for Writing
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
} Population;

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit
);

int fix_unfeasible( double* xr, double lower_limit, double upper_limit );

void blend_crossover( Population* population, int father, int mother, int son, float alpha );

void mutation( Population* population, unsigned int individual, unsigned int current_generation, int exponent );

void local_search( Population* population, double ( *fitness_function )( double*, int n ) );

void genetic_algorithm(  );

#else

extern void genetic_algorithm(  );

#endif

#endif