/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

struct thread_arguments
{
	int thread_id;
	int start;
	int end;
	double *residuum_array;

	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M; 

	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */

	uint64_t number;         /* Number of threads                              */
	uint64_t method;         /* Gauss Seidel or Jacobi method of iteration     */
	uint64_t interlines;     /* matrix size = interlines*8+9                   */
	uint64_t inf_func;       /* inference function                             */
	uint64_t termination;    /* termination condition                          */
	uint64_t term_iteration; /* terminate if iteration number reached          */
	double   term_precision; /* terminate if precision reached                 */

	
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (void *threadarguments)
{
	struct thread_arguments *targuments = (struct thread_arguments *)threadarguments;
	int thread_id = targuments->thread_id;
	int start = targuments->start;
	int end = targuments->end;
	double *residuum_array = targuments->residuum_array;

	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = targuments->N;
	double const h = targuments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = targuments->term_iteration;

	if (targuments->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	/* initialize m1 and m2 depending on algorithm */
	if (targuments->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = targuments->Matrix[m1];
		double** Matrix_In  = targuments->Matrix[m2];

		maxResiduum = 0;

		/* over all rows */
		for (i = start; i < end; i++)
		{
			double fpisin_i = 0.0;

			if (targuments->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (targuments->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (targuments->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
		residuum_array[thread_id] = maxResiduum;

		targuments->stat_iteration++;
		targuments->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (targuments->termination == TERM_PREC)
		{
			if (maxResiduum < targuments->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (targuments->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	targuments->m = m2;
}

static 
void
which_method(struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	struct thread_arguments threadarguments;

	threadarguments.thread_id = 0;
	threadarguments.start = 0;
	threadarguments.end = arguments->N;
	//threadarguments.residuum_array

	threadarguments.N = arguments->N;
	threadarguments.num_matrices = arguments->num_matrices;
	threadarguments.h = arguments->h;
	threadarguments.Matrix = arguments->Matrix;
	threadarguments.M = arguments->M;

	threadarguments.m = results->m;
	threadarguments.stat_iteration = results->stat_iteration;
	threadarguments.stat_precision = results->stat_precision;

	threadarguments.number = options->number;
	threadarguments.method = options->method; 
	threadarguments.interlines = options->interlines;  
	threadarguments.inf_func = options->inf_func;     
	threadarguments.termination = options->termination;
	threadarguments.term_iteration = options->term_iteration;
	threadarguments.term_precision = options->term_precision;
	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		//variables from stuct
		int start = threadarguments.start;
		int end = threadarguments.end;
		int num_threads = threadarguments.number;

		//new fields
		pthread_t threads[num_threads];
		double residuum_array[num_threads];
		double max_residuum = 0;
		
		//calculated fields 
		int parts = (int) threadarguments.N/num_threads;
		int remain = threadarguments.N%num_threads;

		for (int t = 0; t < num_threads; t++)
		{
			struct thread_arguments targs;
			targs.thread_id = t;
			targs.start = end+1;
			//targs.end = end;
			//threadarguments.residuum_array

			targs.N = arguments->N;
			targs.num_matrices = arguments->num_matrices;
			targs.h = arguments->h;
			targs.Matrix = arguments->Matrix;
			targs.M = arguments->M;

			targs.m = results->m;
			targs.stat_iteration = results->stat_iteration;
			targs.stat_precision = results->stat_precision;

			targs.number = options->number;
			targs.method = options->method; 
			targs.interlines = options->interlines;  
			targs.inf_func = options->inf_func;     
			targs.termination = options->termination;
			targs.term_iteration = options->term_iteration;
			targs.term_precision = options->term_precision;

			end = (t+1) * parts;
			if (t < remain){
				end ++;
			}
			targs.end = end;
			targs.residuum_array = &residuum_array;
			pthread_create(&threads[t], NULL, calculate, (void *)&targs);
			//start = end +1;
		}
		for (int t = 0; t < num_threads; t++)
		{
			pthread_join(threads[t], NULL);
		}
		for (int l = 0; l < num_threads; l++)
		{
			if (residuum_array[l] > max_residuum){
				max_residuum = residuum_array[l];
			}
		}
	}
	else
	{
		calculate(&threadarguments);
	}
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	struct thread_arguments targuments;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	which_method(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
