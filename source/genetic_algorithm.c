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
    double sum = 0;
    for ( int i = 0; i < n; i++ )
        sum += genes[ i ];
    return sum;
}

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( double*, int n ), double lower_limit, double upper_limit
){
    int i, j;
    double gene, fit, sum = 0;

    for ( i = 0; i < population_size; i++ ){
        population->reading_individuals[ i ].neighbors[ 0 ] = i % COLS == 0 ? i + ( COLS - 1 ) : i - 1;         // Left
        population->reading_individuals[ i ].neighbors[ 1 ] = i % COLS == COLS - 1 ? i - ( COLS - 1 ) : i + 1;  // Right

        printf( "%d : %d\n\n", i, ( i / COLS ) );
        population->reading_individuals[ i ].neighbors[ 2 ] = ( i / COLS ) > ROWS - 3 ? i - COLS * 2 : i + COLS * 2;         // Top or Bottom
        population->reading_individuals[ i ].neighbors[ 3 ] = ( i / COLS ) == ROWS - 2 || ( i / COLS ) == 0 ? i + COLS : i - COLS * 2;  // Bottom or Top

        // Correction if it's on line 1
        population->reading_individuals[ i ].neighbors[ 3 ] += ( i / COLS ) == 1 ? COLS : 0;

        for ( j = 0; j < 4; j++ )
            population->writing_individuals[ i ].neighbors[ j ] = population->reading_individuals[ i ].neighbors[ j ];

        for ( j = 0; j < individual_size; j++ ){
            gene = (double) randgen( lower_limit, upper_limit );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;
        }

        fit = fitness_function( population->reading_individuals[ i ].genes, individual_size );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fit;

        population->reading_individuals[ i ].RW_type = READING;
        population->writing_individuals[ i ].RW_type = WRITING;

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

void genetic_algorithm(  ){
    Population population;

    start_population( &population, MAX_POPULATION, MAX_GENES, sample_fitness, 0.0F, 10.0F );

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
}

#endif