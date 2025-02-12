#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "genetic_algorithm.h"

#define PI 3.14159265358979F

/* Persistencia */
#define MAXGER 5000
#define NUMCRU 10
#define TAXERR 1.0e-10
#define MAXAVA 10000

// Fitness functions

double griewangk( double x[  ], int n ){
    register int i;
    double sum, pro;
    for ( i = 1, sum = 0, pro = 1; i <= n; i++) {
        sum = sum + pow( x[ i - 1 ], 2) / 4000;
        pro = pro * cos( x[ i - 1 ] / sqrt( i ));
    }

    return ( 1 + sum + pro );
}

double rastringin ( double x[  ], int n ){
    register int i;
    double sum;
    for ( i = 0, sum = 3 * n; i< n; i++ ){
        sum= sum + pow( x[ i ], 2) - 3 * cos( 2 * PI * x[ i ]);
    }

    return sum;
}

double rosenbrock( double x[  ], int n ){
	int i;

	double d = 0;
	for( i = 0; i < n - 1 ; i++ ){
		double p = x[ i + 1 ] - x[ i ] * x[ i ];
		double q = x[ i ] - 1;
		d += 100 * p * p + q * q;
	}
	return d;
}

static const double schw_ans = 418.9828872724336861210758797824382781982;

double schwefel7( double x[  ], int n ){
	int i;
	double d;
	for(i = 0, d = 0; i < n; i++)
		d -= x[i] * sin( sqrt( fabs( x[ i ] ) ) );

	return ( d + schw_ans * ( double )n );
}

double sphere( double* genes, int n ){
    register int i;
    double sum = 100;
    for ( i = 0; i < n; i++ )
        sum += pow( genes[ i ], 2 );

    return sum;
}

// Funções teste
enum { gri, ras, ros, sch, sph } enumFuncs;
static double ( *funccod[  ] )( double[  ], int n ) = {
    griewangk, rastringin, rosenbrock, schwefel7, sphere
};

typedef struct _dados_funcoes_{
    char nome[ 5 ];
    double inf;
    double sup;
} DadosFuncoes;

DadosFuncoes limites_fixos[  ] = {
    { "Gri", -600.0, 600.0 }, { "Ras", -5.12, 5.12 }, { "Ros", -2.048F, 2.048F },
    { "Sch", -500.0F, 500.0F }, { "Sph", -5.12F, 5.12F }
};

int main( int argc, char* argv[  ] ){
    int thread_count;

    int max_vars[ 3 ][ 5 ] = {
        { 10, 10, 10, 10, 10 },
        { 200, 200, 200, 200, 200 },
        { 700, 200, 500, 200, 1000 }
    };

    for ( int j = 0; j < 1; j++ ){
        for ( int i = ras; i <= ras; i++ ){
            // for ( thread_count = 10; thread_count <= 40; thread_count *= 2 ){
            for ( thread_count = 8; thread_count < 9; thread_count++ ){
                printf( "Iniciando execução de %s com %d threads\n", limites_fixos[ i ].nome, thread_count );
                genetic_algorithm( thread_count, max_vars[ j ][ i ], funccod[ i ], limites_fixos[ i ].inf, limites_fixos[ i ].sup );
                printf( "Finalizando execução de %s\n\n", limites_fixos[ i ].nome );
            }
        }
    }



    return 0;
}