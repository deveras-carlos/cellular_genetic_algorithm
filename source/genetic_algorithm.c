#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

float randgen(float lower_limit, float upper_limit){
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
    int i, j;
    double gene, fit, sum = 0;
    int min, max;
    int quintil = population_size / 5;

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

        // printf( "\tVizinhos de %d: [ %d, %d, %d, %d ]\n", i,
        //     population->reading_individuals[ i ].neighbors[ 0 ], population->reading_individuals[ i ].neighbors[ 1 ],
        //     population->reading_individuals[ i ].neighbors[ 2 ], population->reading_individuals[ i ].neighbors[ 3 ]
        // );


        for ( j = 0; j < 4; j++ )
            population->reading_individuals[ i ].neighbors[ j ] = population->writing_individuals[ i ].neighbors[ j ];

        for ( j = 0; j < individual_size; j++ ){
            gene = (double) randgen( lower_limit, upper_limit );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;
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

int tournament_selection(Population* population, int individual, int tournament_size) {
    int best = -1;
    double best_fitness = INFINITY; // Assuming lower fitness is better
    for (int i = 0; i < tournament_size; i++) {
        int rand_index = rand() % population->population_size;
        double fitness = population->reading_individuals[ rand_index ].fitness;
        if (fitness < best_fitness) {
            best = rand_index;
            best_fitness = fitness;
        }
    }
    return best;
}

void blend_crossover( Population* population, int father, int mother, int son, float alpha ){
    double a, b, r;
    int i;

    a = -alpha;
    b = 1 + alpha;

    for ( i = 0; i < population->individual_size; i++ ){
        r = a + ( ( rand(  ) % 101 ) / 100. ) * ( b - a );
        // printf( "Teste: %lf\n", r );
        // printf( "Son: %d, Father: %d, Mother: %d\n", son, father, mother );
        // printf( "Gene em i: %d\n", i );
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ father ].genes[ i ] +
        r * ( population->reading_individuals[ mother ].genes[ i ] - population->reading_individuals[ father ].genes[ i ] );
        // printf( "Teste 2: %lf\n", population->writing_individuals[ son ].genes[ i ] );
        fix_unfeasible( &( population->writing_individuals[ son ].genes[ i ] ), population->lower_limit, population->upper_limit );
    }
}

void downhill_local_search(
    Chromossome* individual, int individual_size, int step,
    double (*fitness_function)(double*, int), double lower_limit, double upper_limit
){
    int i;
    double step_size = step != 0 ? ( 100. / ( double ) step ) : 1;
    // double step_size = 0.0001;
    double new_fitness, original_fitness = individual->fitness;

    for (i = 0; i < individual_size; i++) {
        double original_gene = individual->genes[i];

        individual->genes[i] = original_gene + step_size;
        fix_unfeasible(&(individual->genes[i]), lower_limit, upper_limit);
        new_fitness = fitness_function(individual->genes, individual_size);

        if (new_fitness < original_fitness) {
            original_fitness = new_fitness;
        } else {
            individual->genes[i] = original_gene - step_size;
            fix_unfeasible(&(individual->genes[i]), lower_limit, upper_limit);
            new_fitness = fitness_function(individual->genes, individual_size);

            if (new_fitness < original_fitness) {
                original_fitness = new_fitness;
            } else {
                individual->genes[i] = original_gene;
            }
        }
    }

    individual->fitness = original_fitness;
}

void genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit ){
    
void genetic_algorithm( unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit ){
    
    clock_t start, end;

    Population population;
    int pa1, pa2;
    double fit, gene;
    int cfo = 0, generation = 0;
    int i, j;
    int no_improvement_trials = 0;
    int num_threads_crossover = (int)(population_size * 0.6);
    int num_threads_local_search = population_size - num_threads_crossover;

    char aux_type;

    int max_iter_local_search = MIN_ITER_LOCAL_SEARCH;

    Chromossome auxiliar_individual;
    double *auxiliar_genes;
    int num_threads = population_size;

    int num_threads = 40;

    srand( ( unsigned ) time( 0 ) );

    start_population( &population, population_size, individual_size, fitness_function, lower_limit, upper_limit );
    population.lower_limit = lower_limit;
    population.upper_limit = upper_limit;

    start = clock();
    do {
        int improvement = 0;

        // Crossover and mutation
        #pragma omp parallel for num_threads( num_threads ) \
        private( pa1, pa2, fit, gene, i, j, max_iter_local_search, generation ) shared( population ) reduction( +:cfo )
        for ( i = 0; i < population.population_size; i++ ){
            if ( population.reading_individuals[ i ].type == LOCAL_SEARCH_TYPE ){
                for ( j = 0; j < max_iter_local_search; j++ )
                    downhill_local_search(
                        &(population.reading_individuals[i]), 
                        population.individual_size,
                        generation,
                        fitness_function, 
                        lower_limit, 
                        upper_limit
                    );
            } else if ( population.reading_individuals[ i ].type == RANDOM_TYPE ){
                for ( j = 0; j < population.individual_size; j++ ){
                    gene = (double) randgen( lower_limit, upper_limit );
                    population.reading_individuals[ i ].genes[ j ] = gene;
                    population.writing_individuals[ i ].genes[ j ] = gene;
                }
                fit = fitness_function( population.reading_individuals[ i ].genes, population.individual_size );
                population.writing_individuals[ i ].fitness = population.reading_individuals[ i ].fitness = fit;
            } else {
                pa1 = population.reading_individuals[ i ].neighbors[ rand() % AMT_NEIGHBORS ];
                pa2 = population.reading_individuals[ i ].neighbors[ rand() % AMT_NEIGHBORS ];
                
                if ( pa1 != pa2 ){
                    blend_crossover( &population, pa1, pa2, i, XALPHA );
                    population.writing_individuals[ i ].fitness = fitness_function( population.writing_individuals[ i ].genes, population.individual_size );
                    if ( round(population.writing_individuals[ i ].fitness * 10000) < round(population.reading_individuals[ i ].fitness * 10000) ){
                        population.reading_individuals[ i ].sels++;
                        improvement = 1;
                    }
                    cfo++;
                }
            }
        }

        // Update reading individuals
        #pragma omp parallel for num_threads( num_threads ) \
        private( auxiliar_genes, i ) shared( population )
        for ( i = 0; i < population.population_size; i++ ){
            if (population.reading_individuals[ i ].sels > population.writing_individuals[ i ].sels ){
                population.reading_individuals[ i ].fitness = population.writing_individuals[ i ].fitness;
                auxiliar_genes = population.reading_individuals[ i ].genes;
                population.reading_individuals[ i ].genes = population.writing_individuals[ i ].genes;
                population.writing_individuals[ i ].genes = auxiliar_genes;
                population.writing_individuals[ i ].sels++;

                if ( population.reading_individuals[ i ].fitness - population.writing_individuals[ i ].fitness > 0.00001 )
                    improvement = 1;
            }
        }

        // // Local search
        // #pragma omp parallel for num_threads(num_threads) private(i, j) shared(population)
        // for (i = 0; i < population.population_size; i++) {
        //     for ( j = 0; j < max_iter_local_search; j++ )
        //         downhill_local_search(
        //             &(population.reading_individuals[i]), 
        //             population.individual_size,
        //             generation,
        //             fitness_function, 
        //             lower_limit, 
        //             upper_limit
        //         );
        // }

        #pragma omp parallel for num_threads(num_threads_local_search) private(i, auxiliar_individual, aux_type)
        for (i = 0; i < population.population_size; i++) {
            int swap_idx = rand() % population.population_size;
            if (swap_idx != i && population.reading_individuals[ i ].type != RANDOM_TYPE && population.reading_individuals[ swap_idx ].type != RANDOM_TYPE) {
                aux_type = population.reading_individuals[ swap_idx ].type;

                auxiliar_individual = population.reading_individuals[i];
                population.reading_individuals[i] = population.reading_individuals[swap_idx];
                population.reading_individuals[swap_idx] = auxiliar_individual;

                population.reading_individuals[ swap_idx ].type = aux_type;
                population.reading_individuals[ i ].type = auxiliar_individual.type;

                aux_type = population.writing_individuals[ swap_idx ].type;

                auxiliar_individual = population.writing_individuals[i];
                population.writing_individuals[i] = population.writing_individuals[swap_idx];
                population.writing_individuals[swap_idx] = auxiliar_individual;

                population.writing_individuals[ swap_idx ].type = aux_type;
                population.writing_individuals[ i ].type = auxiliar_individual.type;
            }
        }

        generation++;

        if (improvement || cfo < MIN_ITERATIONS)
            no_improvement_trials = 0;
        else
            no_improvement_trials++;

        // max_iter_local_search = MIN_ITER_LOCAL_SEARCH + ( int )( ( cfo ) / ( 100 ) );

        printf( "Cfo: %d\n", cfo );

    } while ( no_improvement_trials < MAX_WAIT_TRIALS );

    end = clock();

    printf("\n\tExecutado em tempo = %.4f,CFO=%d ", (double) ( (end - start) / (CLOCKS_PER_SEC*num_threads) ), cfo);

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\nFinal max_iter_local_search value: %d\n", max_iter_local_search );

}


#endif
