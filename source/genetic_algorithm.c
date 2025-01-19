#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

float randgen(float lower_limit, float upper_limit)
{
	float fRandomVal;

	fRandomVal = rand()%101/100.;        // rand entre 0 e 1

	return(	lower_limit + (float) (fRandomVal * (upper_limit - lower_limit)) );

}

double sample_fitness( double* genes, int n ){
    register int i;
    double sum = 100;
    for ( i = 0; i< n; i++ )
        sum += pow( genes[i], 2 );

    return sum;
}

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit
){
    int i, j;
    double gene, fit, sum = 0;

    for ( i = 0; i < population_size; i++ ){
        population->reading_individuals[ i ].genes = ( double* ) malloc( MAX_GENES * sizeof( double ) );
        population->writing_individuals[ i ].genes = ( double* ) malloc( MAX_GENES * sizeof( double ) );

        population->reading_individuals[ i ].neighbors[ 0 ] = i % COLS == 0 ? i + ( COLS - 1 ) : i - 1;         // Left
        population->reading_individuals[ i ].neighbors[ 1 ] = i % COLS == COLS - 1 ? i - ( COLS - 1 ) : i + 1;  // Right

        population->reading_individuals[ i ].neighbors[ 2 ] = ( i / COLS ) > ROWS - 3 ? i - COLS * 2 : i + COLS * 2;         // Top or Bottom
        population->reading_individuals[ i ].neighbors[ 3 ] = ( i / COLS ) == ROWS - 2 || ( i / COLS ) == 0 ? i + COLS : i - COLS * 2;  // Bottom or Top

        // Correction if it's on line 1 or last line
        population->reading_individuals[ i ].neighbors[ 3 ] += ( i / COLS ) == 1 || ( i / COLS ) > ROWS - 2 ? COLS : 0;
        // printf( "%d : %d\n\n", i, ( i / COLS ) );

        for ( j = 0; j < 4; j++ )
            population->writing_individuals[ i ].neighbors[ j ] = population->reading_individuals[ i ].neighbors[ j ];

        for ( j = 0; j < individual_size; j++ ){
            gene = (double) randgen( lower_limit, upper_limit );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;
        }

        fit = fitness_function( population->reading_individuals[ i ].genes, individual_size );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fit;

        if ( i == 0 || i == COLS / 4 || i == COLS / 2 || i == COLS - COLS / 4 ){
            population->reading_individuals[ i ].random = TRUE;
            population->writing_individuals[ i ].random = TRUE;
        } else {
            population->reading_individuals[ i ].random = FALSE;
            population->writing_individuals[ i ].random = FALSE;
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
        r = a + ( rand(  ) % 101/100. ) * ( b - a );
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ father ].genes[ i ] +
        r * ( population->reading_individuals[ mother ].genes[ i ] - population->reading_individuals[ father ].genes[ i ] );
        fix_unfeasible( &( population->writing_individuals[ son ].genes[ i ] ), LOWER_BOUND, UPPER_BOUND );
    }
}

void downhill_local_search(
    Chromossome* individual, int individual_size,
    double (*fitness_function)(double*, int), double lower_limit, double upper_limit
){
    int i;
    double step_size = 0.01; // Step size for tweaking genes
    double new_fitness, original_fitness = individual->fitness;

    for (i = 0; i < individual_size; i++) {
        // Store original gene value
        double original_gene = individual->genes[i];

        // Perturb the gene positively
        individual->genes[i] = original_gene + step_size;
        fix_unfeasible(&(individual->genes[i]), lower_limit, upper_limit);
        new_fitness = fitness_function(individual->genes, individual_size);

        // If improved, keep the change; otherwise, try negative perturbation
        if (new_fitness < original_fitness) {
            original_fitness = new_fitness;
        } else {
            // Restore original value and try negative perturbation
            individual->genes[i] = original_gene - step_size;
            fix_unfeasible(&(individual->genes[i]), lower_limit, upper_limit);
            new_fitness = fitness_function(individual->genes, individual_size);

            if (new_fitness < original_fitness) {
                original_fitness = new_fitness;
            } else {
                // Restore original gene if neither direction improves fitness
                individual->genes[i] = original_gene;
            }
        }
    }

    // Update fitness
    individual->fitness = original_fitness;
}

void genetic_algorithm(  ){
    Population population;
    int pa1, pa2;
    double fit, dvp, med, gene;
    int cfo = 0, generation = 0;
    int i, j;

    double *auxiliar_genes;

    int num_threads = 4;

    srand((unsigned) time(0));

    start_population( &population, MAX_POPULATION, MAX_GENES, sample_fitness, LOWER_BOUND, UPPER_BOUND );

    for ( int i = 0; i < MAX_POPULATION; i++ ){
        printf( "Individual[ %d ] = {\n", i );
        printf( "\tFitness: %lf\n", population.reading_individuals[ i ].fitness );
        printf( "\tGenes : [ " );
        for ( int j = 0; j < MAX_GENES - 1; j++ )
            printf( "%lf, ", population.reading_individuals[ i ].genes[ j ]);
        printf( "%lf ]\n", population.reading_individuals[ i ].genes[ MAX_GENES - 1 ] );
        printf( "\tVizinhos: [ %d, %d, %d, %d ]\n",
            population.reading_individuals[ i ].neighbors[ 0 ], population.reading_individuals[ i ].neighbors[ 1 ],
            population.reading_individuals[ i ].neighbors[ 2], population.reading_individuals[ i ].neighbors[ 3 ]
        );
        printf( "}\n" );
    }

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    do {

        #pragma omp parallel for num_threads( num_threads ) \
        private( pa1, pa2, fit, gene, i, j ) shared( population ) reduction( +:cfo )
        for ( i = 0; i < population.population_size; i++ ){
            if ( population.reading_individuals[ i ].random == TRUE ){
                for ( j = 0; j < population.individual_size; j++ ){
                    gene = (double) randgen( LOWER_BOUND, UPPER_BOUND );
                    population.reading_individuals[ i ].genes[ j ] = gene;
                    population.writing_individuals[ i ].genes[ j ] = gene;
                }
            } else {
                pa1 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];
                pa2 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];

                if ( pa1 != pa2 ){
                    blend_crossover( &population, pa1, pa2, i, XALPHA );
                    population.writing_individuals[ i ].fitness = sample_fitness( population.writing_individuals[ i ].genes, population.individual_size );
                    if ( population.writing_individuals[ i ].fitness < population.reading_individuals[ i ].fitness )
                        population.reading_individuals[ i ].sels++;
                    
                    cfo++;
                }
            }
        }

        #pragma omp parallel for num_threads( num_threads ) \
        private( auxiliar_genes, i ) shared( population )
        for ( i = 0; i < population.population_size; i++ ){
            if (population.reading_individuals[ i ].sels > population.writing_individuals[ i ].sels ){
                population.reading_individuals[ i ].fitness = population.writing_individuals[ i ].fitness;

                auxiliar_genes = population.reading_individuals[ i ].genes;
                population.reading_individuals[ i ].genes = population.writing_individuals[ i ].genes;
                population.writing_individuals[ i ].genes = auxiliar_genes;

                population.writing_individuals[ i ].sels++;
            }
        }

        #pragma omp parallel for num_threads(num_threads) private(i, j) shared(population)
        for (i = 0; i < population.population_size; i++) {
            for ( j = 0; j < MAX_ITER_LOCAL_SEARCH; j++ )
                downhill_local_search(
                    &(population.reading_individuals[i]), 
                    population.individual_size, 
                    sample_fitness, 
                    LOWER_BOUND, 
                    UPPER_BOUND
                );
        }

        generation++;
        printf( "Evaluated individuals %d\n", cfo );
        printf( "Generation number %d\n", generation );
    } while ( cfo < MAX_GEN );

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );
}

#endif