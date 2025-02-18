#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <float.h>
#include <sys/stat.h> // For creating the "results" directory
#include <sys/types.h>

void save_results(const Population *population, double execution_time, 
                  unsigned int i, unsigned int max_wait_trials, 
                  unsigned int population_size, unsigned int j, unsigned int curr_amt_vars) {
    // Create the "results" folder if it doesn't exist
    struct stat st = {0};
    if (stat("results", &st) == -1) {
        mkdir("results", 0700); // Creates folder with permissions
    }

    // Generate the filename
    char filename[256];
    snprintf(filename, sizeof(filename), "results/results_i%d_vars_%d_trials%d_pop%d_j%d.txt", 
             i, curr_amt_vars, max_wait_trials, population_size, j);

    // Open the file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        return;
    }

    // Write metadata
    fprintf(file, "Execution Time: %.2f seconds\n", execution_time);
    fprintf(file, "Population Size: %d\n", population->population_size);
    fprintf(file, "Individual Size: %d\n", population->individual_size);
    fprintf(file, "Best Fitness: %.10f\n", 
            population->reading_individuals[population->best].fitness);
    fprintf(file, "Worst Fitness: %.10f\n", 
            population->reading_individuals[population->worst].fitness);
    fprintf(file, "Sum of Fitness: %.10f\n", population->sum_fitness);
    fprintf(file, "Best Generation: %d\n", population->best_generation);

    // Write population data
    fprintf(file, "\nPopulation Data:\n");
    for (int i = 0; i < population->population_size; i++) {
        fprintf(file, "Individual %d: Fitness = %.10f, Type = %d, Sels = %d\n", 
                i, population->reading_individuals[i].fitness, 
                population->reading_individuals[i].type, 
                population->reading_individuals[i].sels);
        fprintf(file, "Genes: ");
        for (int j = 0; j < population->individual_size; j++) {
            fprintf(file, "%.10f ", population->reading_individuals[i].genes[j]);
        }
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
    printf("Results saved to %s\n", filename);
}

float randgen(float lower_limit, float upper_limit){
	float fRandomVal;

	fRandomVal = ( rand(  ) % 101 ) / 100.;
	fRandomVal = ( rand(  ) % 101 ) / 100.;

	return (	lower_limit + (float) (fRandomVal * (upper_limit - lower_limit)) );
	return (	lower_limit + (float) (fRandomVal * (upper_limit - lower_limit)) );

}

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double (*fitness_function)(double*, int n), double lower_limit, double upper_limit
) {
    int i, j;
    int rows = (int)(population_size / 4); // Assuming 4 columns
    int cols = 4; // Fixed number of columns
    double step_size = (upper_limit - lower_limit) / (double)(population_size);
    double sum = 0;

    for (i = 0; i < population_size; i++) {
        population->reading_individuals[i].genes = (double*)malloc(individual_size * sizeof(double));
        population->writing_individuals[i].genes = (double*)malloc(individual_size * sizeof(double));

        // Compute grid position
        int row = i / cols;
        int col = i % cols;

        // Von Neumann neighbors with wrapping edges
        population->writing_individuals[i].neighbors[0] = ((row - 1 + rows) % rows) * cols + col; // Top neighbor
        population->writing_individuals[i].neighbors[1] = ((row + 1) % rows) * cols + col;       // Bottom neighbor
        population->writing_individuals[i].neighbors[2] = row * cols + (col - 1 + cols) % cols;  // Left neighbor
        population->writing_individuals[i].neighbors[3] = row * cols + (col + 1) % cols;        // Right neighbor

        // Alternate rows: swap left and right neighbors
        if (row % 2 == 1) {
            int temp = population->writing_individuals[i].neighbors[2];
            population->writing_individuals[i].neighbors[2] = population->writing_individuals[i].neighbors[3];
            population->writing_individuals[i].neighbors[3] = temp;
        }

        for (j = 0; j < AMT_NEIGHBORS; j++) {
            population->reading_individuals[i].neighbors[j] = population->writing_individuals[i].neighbors[j];
        }

        // Initialize genes with a value around a base value plus some random offset
        for (j = 0; j < individual_size; j++) {
            double base_value = lower_limit + step_size * i;
            double random_offset = ((rand() % 101) / 100.0) * step_size - (step_size / 2.0);
            double gene = base_value + random_offset;

            // Clamp gene value to limits
            if (gene < lower_limit) gene = lower_limit;
            if (gene > upper_limit) gene = upper_limit;

            population->reading_individuals[i].genes[j] = gene;
            population->writing_individuals[i].genes[j] = gene;
        }

        // Assign fitness
        double fit = fitness_function(population->reading_individuals[i].genes, individual_size);
        population->writing_individuals[i].fitness = population->reading_individuals[i].fitness = fit;

        // Assign individual type
        if (i > (int)(population_size * 0.6)) {
            population->reading_individuals[i].type = LOCAL_SEARCH_TYPE;
            population->writing_individuals[i].type = LOCAL_SEARCH_TYPE;
        } else if (i % (population_size / 5) == 0) {
            population->reading_individuals[i].type = RANDOM_TYPE;
            population->writing_individuals[i].type = RANDOM_TYPE;
        } else {
            population->reading_individuals[i].type = CROSSOVER_TYPE;
            population->writing_individuals[i].type = CROSSOVER_TYPE;
        }

        population->reading_individuals[i].sels = 0;
        population->writing_individuals[i].sels = 0;

        sum += fit;
    }

    population->population_size = population_size;
    population->individual_size = individual_size;
    population->sum_fitness = sum;
    population->best = 0;
    population->worst = 0;
    population->mutation_amount = 0;
    population->equals = 0;
}

int fix_unfeasible( double* xr, double lower_limit, double upper_limit ){
    if ( ( *xr ) > upper_limit )
        *xr = upper_limit - PREBATI * ( ( *xr ) - upper_limit ) / ( ( *xr ) - lower_limit );
    else if ( ( *xr ) < lower_limit )
        *xr = lower_limit + PREBATI * ( lower_limit - ( *xr ) ) / ( upper_limit - ( *xr ) );
    
    return ( 1 );
}

void blend_crossover( Population* population, int father, int mother, int son, float alpha ){
    double a, b, r;
    int i;

    a = -alpha;
    b = 1 + alpha;

    for ( i = 0; i < population->individual_size; i++ ){
        r = a + ( ( rand(  ) % 101 ) / 100. ) * ( b - a );
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ father ].genes[ i ] +
        r * ( population->reading_individuals[ mother ].genes[ i ] - population->reading_individuals[ father ].genes[ i ] );
        fix_unfeasible( &( population->writing_individuals[ son ].genes[ i ] ), population->lower_limit, population->upper_limit );
    }
}

void mutation(Population* population, unsigned int individual, unsigned int current_generation, 
              int exponent, double mutation_rate, int no_improvement_trials) {
    int i;

    double base_mutation_range = pow(2, -exponent);
    double adjustment_factor = 1.0 + (no_improvement_trials * 0.05);
    double mutation_range = base_mutation_range * adjustment_factor;

    for (i = 0; i < population->individual_size; i++) {
        if ((rand() % 101) / 100.0 < mutation_rate) {
            double original_gene = population->writing_individuals[individual].genes[i];
            
            double mutation_value = randgen(-mutation_range, mutation_range);
            population->writing_individuals[individual].genes[i] += mutation_value;

            fix_unfeasible(&(population->writing_individuals[individual].genes[i]), 
                            population->lower_limit, population->upper_limit);
        }
    }
}

void mutation_local_search(
    Chromossome* individual, int individual_size,
    double lower_limit, double upper_limit, double (*fitness_function)(double*, int),
    double mutation_rate, double scale_factor, double random_factor
) {
    int i;
    double current_fitness, best_fitness;
    double original_gene, mutated_gene, mutation_range;

    best_fitness = fitness_function(individual->genes, individual_size);

    for (i = 0; i < individual_size; i++) {
        if ((rand() % 101) / 100.0 < mutation_rate) {
            original_gene = individual->genes[i];

            mutation_range = fabs(original_gene) * scale_factor + randgen(-random_factor, random_factor);

            mutated_gene = original_gene + randgen(-mutation_range, mutation_range);

            fix_unfeasible(&mutated_gene, lower_limit, upper_limit);

            individual->genes[i] = mutated_gene;
            current_fitness = fitness_function(individual->genes, individual_size);

            if (current_fitness < best_fitness) {
                best_fitness = current_fitness;
            } else {
                individual->genes[i] = original_gene;
            }
        }
    }

    individual->fitness = best_fitness;
}

void parallel_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars ){
    
    clock_t start, end;

    Population population;
    int pa1, pa2;
    double fit, gene;
    int cfo = 0, generation = 0;
    int i, j;
    int no_improvement_trials = 0;
    int num_threads_crossover = (int)(population_size * 0.8);
    int num_threads_local_search = population_size - num_threads_crossover;

    char aux_type;

    int max_iter_local_search = MIN_ITER_LOCAL_SEARCH;

    Chromossome auxiliar_individual;
    double *auxiliar_genes;
    int num_threads = population_size;

    double mutation_rate = PMUTAC / 100.0;

    srand( ( unsigned ) time( 0 ) );

    start_population( &population, population_size, individual_size, fitness_function, lower_limit, upper_limit );
    population.lower_limit = lower_limit;
    population.upper_limit = upper_limit;

    start = clock();
    do {
        int improvement = 0;

        #pragma omp parallel sections
        {
            // Section 1: Local search operations
            #pragma omp section
            {
                #pragma omp parallel for num_threads(num_threads_local_search) \
                private(j) shared(population, improvement)
                for (int i = num_threads_crossover; i < population.population_size; i++) {
                    if (population.reading_individuals[i].type == LOCAL_SEARCH_TYPE) {
                        for (j = 0; j < max_iter_local_search; j++) {
                            mutation_local_search(
                                &(population.writing_individuals[i]),
                                population.individual_size,
                                population.lower_limit,
                                population.upper_limit,
                                fitness_function,
                                mutation_rate,
                                1,
                                (upper_limit - lower_limit) / 100.0
                            );
                        }
                        if (round(population.writing_individuals[i].fitness * 10000000) < 
                            round(population.reading_individuals[i].fitness * 10000000)) {
                            population.writing_individuals[i].sels++;
                            improvement = 1;
                        }
                    }
                }
            }

            // Section 2: Crossover and mutation operations
            #pragma omp section
            {
                #pragma omp parallel for num_threads(num_threads_crossover) \
                private(pa1, pa2, fit) shared(population, generation, cfo, no_improvement_trials)
                for (int i = 0; i < num_threads_crossover; i++) {
                    if (population.reading_individuals[i].type == RANDOM_TYPE) {
                        for (j = 0; j < population.individual_size; j++) {
                            gene = (double)randgen(lower_limit, upper_limit);
                            population.reading_individuals[i].genes[j] = gene;
                            population.writing_individuals[i].genes[j] = gene;
                        }
                        fit = fitness_function(population.reading_individuals[i].genes, 
                                            population.individual_size);
                        population.writing_individuals[i].fitness = fit;
                        population.reading_individuals[i].fitness = fit;
                    } else if (population.reading_individuals[i].type == CROSSOVER_TYPE) {
                        pa1 = population.reading_individuals[i].neighbors[rand() % AMT_NEIGHBORS];
                        pa2 = population.reading_individuals[i].neighbors[rand() % AMT_NEIGHBORS];

                        if (pa1 != pa2) {
                            blend_crossover(
                                &population,
                                pa1,
                                pa2,
                                i,
                                XALPHA * (1.0 - (double)generation / MAX_GEN)
                            );
                            population.writing_individuals[i].fitness = fitness_function(
                                population.writing_individuals[i].genes,
                                population.individual_size
                            );

                            mutation(&population, i, generation, 2, mutation_rate, no_improvement_trials);

                            if (round(population.writing_individuals[i].fitness * 10000000) <
                                round(population.reading_individuals[i].fitness * 10000000)) {
                                population.writing_individuals[i].sels++;
                                improvement = 1;
                            }
                            cfo++;
                        }
                    }
                }
            }
        }

        // Atualização de indivíduos e atualização de tipos
        #pragma omp parallel for num_threads(num_threads) \
        private(auxiliar_genes) shared(population)
        for (int i = 0; i < population.population_size; i++) {
            if (population.writing_individuals[i].sels > population.reading_individuals[i].sels) {
                population.reading_individuals[i].fitness = population.writing_individuals[i].fitness;
                auxiliar_genes = population.reading_individuals[i].genes;
                population.reading_individuals[i].genes = population.writing_individuals[i].genes;
                population.writing_individuals[i].genes = auxiliar_genes;
                population.reading_individuals[i].sels++;
            }

            double rand_change_type = ((rand() % 101) / 100.0);
            int new_pos;
            char old_type, new_type;
            if ( population.reading_individuals[i].type != RANDOM_TYPE && rand_change_type < CHANGE_TYPE_CHANCE) {
                if (population.reading_individuals[i].type == CROSSOVER_TYPE) {
                    new_pos = (rand() % num_threads_crossover);
                } else if (population.reading_individuals[i].type == LOCAL_SEARCH_TYPE) {
                    new_pos = ( ( rand() % num_threads_local_search ) + num_threads_crossover);
                }
                
                if ( population.reading_individuals[ new_pos ].type != RANDOM_TYPE ){
                    #pragma omp critical
                    {
                    Chromossome temp = population.reading_individuals[i];
                    population.reading_individuals[i].genes = population.reading_individuals[new_pos].genes;
                    population.reading_individuals[i].fitness = population.reading_individuals[new_pos].fitness;
                    population.reading_individuals[new_pos].genes = temp.genes;
                    population.reading_individuals[new_pos].fitness = temp.fitness;

                    temp = population.writing_individuals[i];
                    population.writing_individuals[i].genes = population.writing_individuals[new_pos].genes;
                    population.writing_individuals[i].fitness = population.writing_individuals[new_pos].fitness;
                    population.writing_individuals[new_pos].genes = temp.genes;
                    population.writing_individuals[new_pos].fitness = temp.fitness;
                    }
                }
            }
        }

        if (improvement || generation < MIN_ITERATIONS) {
            no_improvement_trials = 0;
            mutation_rate = PMUTAC / 100.0;
        } else {
            no_improvement_trials++;
            mutation_rate = fmin(mutation_rate + 0.05, 1.0);
        }


        generation++;
        printf( "Generation: %d\n", generation );
    } while ( generation < max_wait_trials );
    end = clock();

    // Calculate execution time
    double execution_time = ( (double)(end - start) / CLOCKS_PER_SEC ) / ( ( double )num_threads );

    // Update population statistics (e.g., best, worst, sum of fitness)
    population.best = 0;
    population.worst = 0;
    population.sum_fitness = 0.0;
    for (int i = 0; i < population.population_size; i++) {
        population.sum_fitness += population.reading_individuals[i].fitness;
        if (population.reading_individuals[i].fitness < 
            population.reading_individuals[population.best].fitness) {
            population.best = i;
        }
        if (population.reading_individuals[i].fitness > 
            population.reading_individuals[population.worst].fitness) {
            population.worst = i;
        }
    }
    population.best_generation = generation;

    // Save results
    save_results(&population, execution_time, curr_fitness, max_wait_trials, population_size, curr_exec, curr_amt_vars);

    // Free allocated memory
    free(population.reading_individuals->genes);
    free(population.writing_individuals->genes);
}

void sequential_genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit, int max_wait_trials,
    int curr_exec, int curr_fitness, int curr_amt_vars ){
    
    clock_t start, end;

    Population population;
    int pa1, pa2;
    double fit, gene;
    int cfo = 0, generation = 0;
    int i, j;
    int no_improvement_trials = 0;
    int num_threads_crossover = (int)(population_size * 0.8);
    int num_threads_local_search = population_size - num_threads_crossover;

    char aux_type;

    int max_iter_local_search = MIN_ITER_LOCAL_SEARCH;

    Chromossome auxiliar_individual;
    double *auxiliar_genes;
    int num_threads = population_size;

    double mutation_rate = PMUTAC / 100.0;

    srand( ( unsigned ) time( 0 ) );

    start_population( &population, population_size, individual_size, fitness_function, lower_limit, upper_limit );
    population.lower_limit = lower_limit;
    population.upper_limit = upper_limit;

    start = clock();
    do {
        int improvement = 0;

        // Section 1: Local search operations
        for (int i = num_threads_crossover; i < population.population_size; i++) {
            if (population.reading_individuals[i].type == LOCAL_SEARCH_TYPE) {
                for (j = 0; j < max_iter_local_search; j++) {
                    mutation_local_search(
                        &(population.writing_individuals[i]),
                        population.individual_size,
                        population.lower_limit,
                        population.upper_limit,
                        fitness_function,
                        mutation_rate,
                        1,
                        (upper_limit - lower_limit) / 100.0
                    );
                }
                if (round(population.writing_individuals[i].fitness * 10000000) < 
                    round(population.reading_individuals[i].fitness * 10000000)) {
                    population.writing_individuals[i].sels++;
                    improvement = 1;
                }
            }
        }

        // Section 2: Crossover and mutation operations
        for (int i = 0; i < num_threads_crossover; i++) {
            if (population.reading_individuals[i].type == RANDOM_TYPE) {
                for (j = 0; j < population.individual_size; j++) {
                    gene = (double)randgen(lower_limit, upper_limit);
                    population.reading_individuals[i].genes[j] = gene;
                    population.writing_individuals[i].genes[j] = gene;
                }
                fit = fitness_function(population.reading_individuals[i].genes, 
                                    population.individual_size);
                population.writing_individuals[i].fitness = fit;
                population.reading_individuals[i].fitness = fit;
            } else if (population.reading_individuals[i].type == CROSSOVER_TYPE) {
                pa1 = population.reading_individuals[i].neighbors[rand() % AMT_NEIGHBORS];
                pa2 = population.reading_individuals[i].neighbors[rand() % AMT_NEIGHBORS];

                if (pa1 != pa2) {
                    blend_crossover(
                        &population,
                        pa1,
                        pa2,
                        i,
                        XALPHA * (1.0 - (double)generation / MAX_GEN)
                    );
                    population.writing_individuals[i].fitness = fitness_function(
                        population.writing_individuals[i].genes,
                        population.individual_size
                    );

                    mutation(&population, i, generation, 2, mutation_rate, no_improvement_trials);

                    if (round(population.writing_individuals[i].fitness * 10000000) <
                        round(population.reading_individuals[i].fitness * 10000000)) {
                        population.writing_individuals[i].sels++;
                        improvement = 1;
                    }
                    cfo++;
                }
            }
        }

        // Atualização de indivíduos e atualização de tipos
        for (int i = 0; i < population.population_size; i++) {
            if (population.writing_individuals[i].sels > population.reading_individuals[i].sels) {
                population.reading_individuals[i].fitness = population.writing_individuals[i].fitness;
                auxiliar_genes = population.reading_individuals[i].genes;
                population.reading_individuals[i].genes = population.writing_individuals[i].genes;
                population.writing_individuals[i].genes = auxiliar_genes;
                population.reading_individuals[i].sels++;
            }

            double rand_change_type = ((rand() % 101) / 100.0);
            int new_pos;
            char old_type, new_type;
            if ( population.reading_individuals[i].type != RANDOM_TYPE && rand_change_type < CHANGE_TYPE_CHANCE) {
                if (population.reading_individuals[i].type == CROSSOVER_TYPE) {
                    new_pos = (rand() % num_threads_crossover);
                } else if (population.reading_individuals[i].type == LOCAL_SEARCH_TYPE) {
                    new_pos = ( ( rand() % num_threads_local_search ) + num_threads_crossover);
                }
                
                // Swap reading_individuals[i] with reading_individuals[new_pos]
                Chromossome temp = population.reading_individuals[i];
                population.reading_individuals[i].genes = population.reading_individuals[new_pos].genes;
                population.reading_individuals[i].fitness = population.reading_individuals[new_pos].fitness;
                population.reading_individuals[new_pos].genes = temp.genes;
                population.reading_individuals[new_pos].fitness = temp.fitness;

                temp = population.writing_individuals[i];
                population.writing_individuals[i].genes = population.writing_individuals[new_pos].genes;
                population.writing_individuals[i].fitness = population.writing_individuals[new_pos].fitness;
                population.writing_individuals[new_pos].genes = temp.genes;
                population.writing_individuals[new_pos].fitness = temp.fitness;
            }
        }

        if (improvement || generation < MIN_ITERATIONS) {
            no_improvement_trials = 0;
            mutation_rate = PMUTAC / 100.0;
        } else {
            no_improvement_trials++;
            mutation_rate = fmin(mutation_rate + 0.05, 1.0);
        }

        generation++;
    } while ( generation < max_wait_trials );


    end = clock();

    // Calculate execution time
    double execution_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Update population statistics (e.g., best, worst, sum of fitness)
    population.best = 0;
    population.worst = 0;
    population.sum_fitness = 0.0;
    for (int i = 0; i < population.population_size; i++) {
        population.sum_fitness += population.reading_individuals[i].fitness;
        if (population.reading_individuals[i].fitness > 
            population.reading_individuals[population.best].fitness) {
            population.best = i;
        }
        if (population.reading_individuals[i].fitness < 
            population.reading_individuals[population.worst].fitness) {
            population.worst = i;
        }
    }
    population.best_generation = generation;

    // Save results
    save_results(&population, execution_time, curr_fitness, max_wait_trials, population_size, curr_exec, curr_amt_vars);

    // Free allocated memory
    free(population.reading_individuals->genes);
    free(population.writing_individuals->genes);
}

#endif
