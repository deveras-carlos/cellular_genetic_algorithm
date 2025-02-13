#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

float randgen(float lower_limit, float upper_limit){
	float fRandomVal;

	fRandomVal = ( rand(  ) % 101 ) / 100.;
	fRandomVal = ( rand(  ) % 101 ) / 100.;

	return (	lower_limit + (float) (fRandomVal * (upper_limit - lower_limit)) );
	return (	lower_limit + (float) (fRandomVal * (upper_limit - lower_limit)) );

}

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit
){
    int i, j, k;
    double gene, base_value, random_offset;
    double fit, sum = 0;
    int quintil = population_size / 5;
    float greedy_n;

    double step_size = (upper_limit - lower_limit) / (double)(population_size);

    for ( i = 0; i < population_size; i++ ){
        population->reading_individuals[ i ].genes = ( double* ) malloc( individual_size * sizeof( double ) );
        population->writing_individuals[ i ].genes = ( double* ) malloc( individual_size * sizeof( double ) );

        if ( i > 1 && i < population_size - 2 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = i + 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = i + 2;
        } else if ( i == 0 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = population_size - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = population_size - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = i + 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = i + 2;
        } else if ( i == 1 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = population_size - 1;
            population->writing_individuals[ i ].neighbors[ 1 ] = 0;
            population->writing_individuals[ i ].neighbors[ 2 ] = i + 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = i + 2;
        } else if ( i == population_size - 1 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = 0;
            population->writing_individuals[ i ].neighbors[ 3 ] = 1;
        } else if ( i == population_size - 2 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = population_size - 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = 0;
        }


        for ( j = 0; j < 4; j++ )
            population->reading_individuals[ i ].neighbors[ j ] = population->writing_individuals[ i ].neighbors[ j ];
        

        for (j = 0; j < individual_size; j++) {
            base_value = lower_limit + step_size * i;
            random_offset = ((rand() % 101) / 100.0) * step_size - (step_size / 2.0);
            gene = base_value + random_offset;

            // Clamp gene value to limits
            if (gene < lower_limit) gene = lower_limit;
            if (gene > upper_limit) gene = upper_limit;

            population->reading_individuals[i].genes[j] = gene;
            population->writing_individuals[i].genes[j] = gene;
        }


        fit = fitness_function( population->reading_individuals[ i ].genes, individual_size );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fit;

        if ( i > ( int )( population->population_size * 0.6 ) ){
            population->reading_individuals[ i ].type = LOCAL_SEARCH_TYPE;
            population->writing_individuals[ i ].type = LOCAL_SEARCH_TYPE;
        } else if ( i == quintil || i == ( 2 * quintil ) || i == ( 3 * quintil ) || i == ( 4 * quintil )  ){
            population->reading_individuals[ i ].type = RANDOM_TYPE;
            population->writing_individuals[ i ].type = RANDOM_TYPE;
        } else {
            population->reading_individuals[ i ].type = CROSSOVER_TYPE;
            population->writing_individuals[ i ].type = CROSSOVER_TYPE;
        }

        population->reading_individuals[ i ].sels = 0;
        population->writing_individuals[ i ].sels = 0;

        sum += ( fit ); // fabs
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


void local_search(
    Chromossome* individual, int individual_size,
    double lower_limit, double upper_limit, double (*fitness_function)(double*, int),
    int cfo, int max_cfo
) {
    int j;
    double best_gene;
    double best_fitness;
    double current_fitness;

    // Adjust step size and neighborhood size based on cfo
    double step_size = (upper_limit - lower_limit) / (100.0 + cfo * 0.1);
    double neighborhood_size = (upper_limit - lower_limit) * (0.5 - (double)cfo / ( 2 *  max_cfo ) );

    for (j = 0; j < individual_size; j++) {
        double current_gene = individual->genes[j];
        double lower_bound = fmax(lower_limit, current_gene - neighborhood_size);
        double upper_bound = fmin(upper_limit, current_gene + neighborhood_size);

        best_gene = current_gene;
        best_fitness = fitness_function(individual->genes, individual_size);

        for (double candidate = lower_bound; candidate <= upper_bound; candidate += step_size) {
            individual->genes[j] = candidate;
            current_fitness = fitness_function(individual->genes, individual_size);

            if (current_fitness < best_fitness) {
                best_fitness = current_fitness;
                best_gene = candidate;
            }
        }

        individual->genes[j] = best_gene;
    }

    // Update fitness after the local search
    individual->fitness = fitness_function(individual->genes, individual_size);
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
    
void genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit ){
    
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

    printf("\n\tExecutado em tempo = %.4f,CFO=%d ", (double) ( (end - start) / (CLOCKS_PER_SEC*num_threads) ), cfo);

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\nFinal max_iter_local_search value: %d\n", max_iter_local_search );

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
                        if (round(population.writing_individuals[i].fitness * 10000) < 
                            round(population.reading_individuals[i].fitness * 10000)) {
                            population.writing_individuals[i].sels++;
                            improvement = 1;
                        }
                    }
                }
            }

            // Section 3: Crossover and mutation operations
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

                            if (round(population.writing_individuals[i].fitness * 10000) <
                                round(population.reading_individuals[i].fitness * 10000)) {
                                population.writing_individuals[i].sels++;
                                improvement = 1;
                            }
                            cfo++;
                        }
                    }
                }
            }
        }

        #pragma omp parallel for num_threads(num_threads) \
        private(auxiliar_genes) shared(population)
        for (int i = 0; i < population.population_size; i++) {
            if (population.writing_individuals[i].sels > population.reading_individuals[i].sels) {
                #pragma omp critical
                {
                    population.reading_individuals[i].fitness = population.writing_individuals[i].fitness;
                    auxiliar_genes = population.reading_individuals[i].genes;
                    population.reading_individuals[i].genes = population.writing_individuals[i].genes;
                    population.writing_individuals[i].genes = auxiliar_genes;
                    population.reading_individuals[i].sels++;
                }
            }

            double rand_change_type = ((rand() % 101) / 100.0);
            if (population.reading_individuals[i].type != RANDOM_TYPE && 
                rand_change_type < CHANGE_TYPE_CHANCE) {
                char new_type = (rand() % LOCAL_SEARCH_TYPE) + CROSSOVER_TYPE;

                #pragma omp critical
                {
                    population.reading_individuals[i].type = new_type;
                    population.writing_individuals[i].type = new_type;
                }
            }
        }

        // Update mutation rate and other parameters
        if (improvement || generation < MIN_ITERATIONS) {
            no_improvement_trials = 0;
            mutation_rate = PMUTAC / 100.0;
        } else {
            no_improvement_trials++;
            mutation_rate = fmin(mutation_rate + 0.05, 1.0); // Increase by 5%, max 100%
        }

        printf("Mutation rate: %lf\n", mutation_rate);
        printf("Generation: %d\n", generation);

        generation++;
    } while (no_improvement_trials < MAX_WAIT_TRIALS);


    end = clock();


    printf("\n\tExecutado em tempo = %.4f,CFO=%d ", (double) ( (end - start) / (CLOCKS_PER_SEC*num_threads) ), cfo);

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\nFinal max_iter_local_search value: %d\n", max_iter_local_search );

    free( population.reading_individuals->genes );
    free( population.writing_individuals->genes );
}


#endif
