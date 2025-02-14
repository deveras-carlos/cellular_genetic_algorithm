#ifndef _GENETIC_STRUCTURE_H_
#define _GENETIC_STRUCTURE_H_

#define FALSE 0
#define TRUE 1

#define MAX_POPULATION 40
#define AMT_NEIGHBORS 4
#define MAX_GEN 1000

#define PREBATI (rand()%101/100.F)
#define XALPHA  0.2
#define PMUTAC 30

#define GREEDY_START_RATE 0.3
#define MIN_ITER_LOCAL_SEARCH 10
#define MIN_ITERATIONS 1000

#define MAX_WAIT_TRIALS 1000

#define CHANGE_TYPE_CHANCE 0.05
#define AMT_TYPE_IND 3
#define LOCAL_SEARCH_TYPE 2
#define CROSSOVER_TYPE 1
#define RANDOM_TYPE 0

typedef struct _chromossome_ {
    double* genes;
    unsigned int neighbors[ AMT_NEIGHBORS ];
    double fitness;
    char type;
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

#endif