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
#include <mpi.h>

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

//eigenes Struct für MPI Variablen

struct mpi_parameters
{
	int world_size;
	int world_rank;
	uint64_t start_r;
	uint64_t end_r;
	uint64_t N_r;
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
allocateMatrices (struct calculation_arguments* arguments, struct mpi_parameters* mpi_paras)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
	uint64_t const N_r = mpi_paras->N_r;

	// jeder Prozess allokiert Abschnitt der Matrix (abhängig von eigener Größe N_r)

	arguments->M = allocateMemory(arguments->num_matrices * (N_r + 2) * (N + 1) * sizeof(double));

	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N_r + 2) * sizeof(double*));

		for (j = 0; j <= N_r + 1; j++) 
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N_r + 2)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, struct mpi_parameters* mpi_paras)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	uint64_t const N_r = mpi_paras->N_r;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		// Zeilen angepasst
		for (i = 0; i <= N_r + 1; i++)
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
		// Start des Abschnitts wird für Randinitialisierung benötigt
		uint64_t start_r = mpi_paras->start_r;
		for(i = 0; i < N_r + 2; i++)
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1 + (1 - (h * (i+start_r-1))); // Linke Kante
				// Formel geändert, aber äquivalent soweit wir wissen
				Matrix[j][i][N] = 1 - h * (i+start_r-1); // Rechte Kante
				
				//Matrix[j][i][N] = h * (i+start_r-1); // Rechte Kante
			}
		}

		// Oberer/unterer Rand werden nur vom ersten/letzten Prozess initialisiert

		if (mpi_paras->world_rank == 0)
		{
			for(i = 0; i < N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
				}
			}
		}
		if (mpi_paras->world_rank == mpi_paras->world_size - 1)
		{
			for(i = 0; i < N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][N_r + 1][i] = 1 - (h * i); // Untere Kante
				}
			}
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation for one process                           */
/* ************************************************************************ */

// alte calculate wird nur für einen Prozess verwendet
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}


// calculate für paralleles Jacobi (aus voheriger Übung)
/* ************************************************************************ */
/* calcJacobiMPI: solves the equation for Jacobi with MPI                   */
/* ************************************************************************ */
static
void
calcJacobiMPI (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct mpi_parameters* mpi_paras)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	// Höhe (Anzahl Zeilen)
	int const N_r = mpi_paras->N_r;
	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	// MPI Parameter
	int rank = mpi_paras->world_rank;
	int size = mpi_paras->world_size;
	int start_r = mpi_paras->start_r;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 for Jacobi */
	m1 = 0;
	m2 = 1;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;
		/* over all rows */
		// nur über Abschnitt des Prozesses
		for (i = 1; i < N_r + 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				// i_r ist relativ zur gesamten Matrix
				int i_r = i + start_r - 1;

				fpisin_i = fpisin * sin(pih * (double)i_r);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
		
		// erste/letzte Zeilen werden untereinander ausgetauscht, sodass aktuelle (neu berechnete) Werte
		// in der nächsten Iteration verwendet werden
		// Wir vermuten hier ist ein Fehler.
		
		if (rank == 0)
		{
			MPI_Ssend(&Matrix_Out[N_r][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		}
		else 
		{
			MPI_Recv(&Matrix_Out[0][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank < size - 1)
			{
				MPI_Ssend(&Matrix_Out[N_r][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			}
		}

		if (rank == size - 1)
		{
			MPI_Ssend(&Matrix_Out[1][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}
		else 
		{
			MPI_Recv(&Matrix_Out[N_r + 1][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank > 0)
			{
				MPI_Ssend(&Matrix_Out[1][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			}
		}
		
		// alle Prozesse gleiches Residuum
		double allMaxResiduum;
		MPI_Allreduce(&maxResiduum, &allMaxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		maxResiduum = allMaxResiduum;

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}


/* *************************************************************************************************************** */
/* calcGaussSeidelMPI_iter: solves the equation for Gauß-Seidel with MPI for termination after n iterations        */
/* *************************************************************************************************************** */

// neue calculate wird für Gauß-Seidel verwendet (iterativ und nach Genauigkeit haben eigene Funktion)
static
void
calcGaussSeidelMPI_iter (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct mpi_parameters* mpi_paras)
{
	int i, j;           /* local variables for loops */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	// Höhe (Anzahl Zeilen)
	int const N_r = mpi_paras->N_r;

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	// MPI Parameter
	int rank = mpi_paras->world_rank;
	int size = mpi_paras->world_size;
	int start_r = mpi_paras->start_r;

	int term_iteration = options->term_iteration;

	double** Matrix = arguments->Matrix[0];


	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}


	while (term_iteration > 0)
	{

		// Der erste Prozess kann sofort anfangen, alle anderen müssen die voherige Zeile empfangen
		if (rank > 0)
		{
			MPI_Recv(&Matrix[0][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Sobald wir in der 2. Iteration sind, müssen alle außer der letzte Prozess die erste Zeile vom nachfolgenden Prozess empfangen
		if (results->stat_iteration > 0 && rank < size - 1)
		{
			MPI_Recv(&Matrix[N_r + 1][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		maxResiduum = 0;

		/* over all rows, nur über Teilmatrix */
		for (i = 1; i < N_r + 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				// i_r ist relativ zur gesamten Matrix
				int i_r = i + start_r - 1;
				fpisin_i = fpisin * sin(pih * (double)i_r);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix[i][j] = star;
			}
		}
		
		// alle außer der letzte Prozess senden ihre letzte Zeile an den nächsten
		if (rank < size - 1)
		{
			MPI_Ssend(&Matrix[N_r][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		}

		// solange wir noch vor der letzten Iteration sind,
		// wird die neue erste Zeile an den voherigen Prozess gesendet
		if (rank > 0 && term_iteration > 1) 
		{
			MPI_Ssend(&Matrix[1][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;
		
		if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
	results->m = 0;
}

/* ************************************************************************************************************* */
/* calcGaussSeidelMPI_prec: solves the equation for Gauss-Seidel with termination after certain precision        */
/* ************************************************************************************************************* */

// neue calculate wird für Gauß-Seidel verwendet (iterativ und nach Genauigkeit haben eigene Funktion)
static
void
calcGaussSeidelMPI_prec (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct mpi_parameters* mpi_paras)
{
	int i, j;           /* local variables for loops */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	// Höhe (Anzahl Zeilen)
	int const N_r = mpi_paras->N_r;
	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	// MPI Parameter
	int rank = mpi_paras->world_rank;
	int size = mpi_paras->world_size;
	int start_r = mpi_paras->start_r;
	double maybeMaxResiduum;

	int term_iteration = options->term_iteration;

	double** Matrix = arguments->Matrix[0];

	// Variablen für Terminierung
	int signal = 0;
	int gesendet = 0;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		maxResiduum = 0;
		uint16_t aktuelleIteration = results->stat_iteration;
		//printf("Rang %d ist in Iteration %d\n", rank, aktuelleIteration);
	
		// Der erste Prozess kann sofort anfangen, alle anderen müssen die voherige Zeile empfangen
		// bzw. das maxResiduum vom voherigen Prozess
		if (rank > 0)
		{
			MPI_Recv(&Matrix[0][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// das maxResiduum vom voherigen Prozess wird zunächst zwischengespeichert
			MPI_Recv(&maybeMaxResiduum, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// das Signal, dass den Abbruch übermittelt wird ebenfalls immer übertragen
			MPI_Recv(&signal, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		// Sobald wir so viele Iterationen wie Prozesse haben, kann der erste Prozess ein Signal vom letzten empfangen
		if (rank == 0 && aktuelleIteration >= size - 1)
		{
			MPI_Recv(&signal, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// wenn das übertragene Signal sagt, es soll abgebrochen werden, ist diese Iteration die letzte
		if (signal == 1)
		{
			term_iteration = 0;
		}

		// wir erhalten noch die erste Zeile vom nachfolgenden Prozess, damit die Daten aktuell sind
		if (aktuelleIteration > 0 && rank < size - 1)
		{
			MPI_Recv(&Matrix[N_r + 1][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		//printf("Rang %d berechnet Iteration %d\n", rank, aktuelleIteration);


		/* over all rows */
		for (i = 1; i < N_r + 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				// i_r ist relativ zur gesamten Matrix
				int i_r = i + start_r - 1;
				fpisin_i = fpisin * sin(pih * (double)i_r);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
					
				}

				Matrix[i][j] = star;
			}

		}
		
		// Das Residuum wird mit dem des voherigen Prozess verglichen
		maxResiduum = (maybeMaxResiduum > maxResiduum) ? maybeMaxResiduum : maxResiduum;

		// der letzte Rank hat somit das Residuum der gesamten Matrix und kann feststellen, ob abgebrochen werden muss
		if (rank == size - 1)
		{
			if (maxResiduum < options->term_precision)
			{
				signal = 1;
			}
		}
		
		// alle außer der letzte Prozess senden ihre letzte Zeile an den nächsten, sowie Residuum und das Signal
		if (rank < size - 1)
		{
			MPI_Ssend(&Matrix[N_r][0], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Ssend(&maxResiduum, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Ssend(&signal, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
		}

		// wenn der letzte Prozess noch nicht den Abbruch (signal = 1) gesendet hat, so sendet er immer wieder das signal
		// an den ersten Prozess. Durch die wiederholte Weitergabe vom Signal werden die Iterationen dann schrittweise abgebrochen.
		if (rank == size - 1 && gesendet == 0)
		{
			MPI_Ssend(&signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			// wenn das Signal zum Abbruch gesendet wurde, wird dies gespeichert
			if (signal == 1)
			{
				gesendet = 1;
			}
		}
		
		// der letzte Prozess darf noch nicht abbrechen und tut daher so, als hätte er noch kein Signal erhalten
		signal = 0;

		// wir senden die erste Zeile an den voherigen Prozess, wenn dies nicht die letzte Iteration ist
		if (rank > 0 && term_iteration != 0)
		{
			MPI_Ssend(&Matrix[1][0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}
		
		results->stat_iteration++;
		results->stat_precision = maxResiduum;

	}
	results->m = 0;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_results const* results, struct options const* options, double gesamtspeicher)
{
	//int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	// Der Gesamtspeicher wird vorher berechnet un übergeben
	
	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", gesamtspeicher);
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
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, struct mpi_parameters* mpi_paras)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	uint64_t start_r = mpi_paras->start_r;
	uint64_t end_r = mpi_paras->end_r;


	// Nur einmalig Überschrift
	if (mpi_paras->world_rank == 0)
	{
		printf("Matrix:\n");
		// erste Zeile
		for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[0][x * (interlines + 1)]);
			}
			printf ("\n");
	}

	for (y = 1; y < 8; y++)
	{
		// Zeilen dazwischen, wir berechnen "rückwärts" wo wir in der Gesamtmatrix sind
		uint64_t line = y * (interlines + 1);
		if (line < end_r && line >= start_r)
		{
			for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[line - start_r + 1][x * (interlines + 1)]);
			}
			printf ("\n");
		}
	}

	if (mpi_paras->world_rank == mpi_paras->world_size - 1)
	{
		// letzte Zeile
		for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[mpi_paras->N_r + 1][x * (interlines + 1)]);
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

	/*
	einige Erklärungen:
	Die Matrix wird sowohl für Gauss-Seidel als auch für Jacobi in Zeilenabschnitte
	aufgeteilt, sodass jeder Prozess in etwa gleich viele Zeilen enthält.

	Das parallele Jacobi-Verfahren wurde von der voherigen Übung übernommen, lediglich
	die Berechnung von "fpisin_i" wurde korrigiert.

	Für die Berechnung mit einem Prozess wurde die ursprüngliche calculate Funktion verwendet.

	Die parallele Berechnung von Gauss-Seidel wurde aufgeteilt in
	zwei Funktionen für die zwei Abbruchkriterien.

	Das Verfahren, dass nach Iterationszahl abbricht, sollte einen deutlichen Laufzeitvorteil haben.

	Das Verfahren nach Genauigkeit hat einen relativ schwachen speedup, es ist jedoch schneller.
	Hierfür haben wir jedoch getestet, dass das Verhältnis von interlines und nprocs passend sein sollte
	sei interlines ≈ x * nprocs, so sollte für x in etwa gelten: 3 < x < 9
	Sonst braucht das parallele Programm aufgrund der Kommunikation oft mehr Zeit.
	*/

	// MPI initialisieren
	MPI_Init(NULL, NULL);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	struct mpi_parameters mpi_paras;

	// get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_paras.world_size);
    // get rank of process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_paras.world_rank);

	// Signal für displayMatrix Reihenfolge
	int signal = 0;

	// Beschreibung soll nur einmal ausgegeben werden
	if (mpi_paras.world_rank == 0)
	{
	printf("============================================================\n");
	printf("Program for calculation of partial differential equations.  \n");
	printf("============================================================\n");
	printf("(c) Dr. Thomas Ludwig, TU München.\n");
	printf("    Thomas A. Zochler, TU München.\n");
	printf("    Andreas C. Schmidt, TU München.\n");
	printf("============================================================\n");
	printf("\n");
	}

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	// Position und Länge in der Gesamtmatrix
	mpi_paras.start_r = (mpi_paras.world_rank*(arguments.N - 1))/mpi_paras.world_size + 1;
	mpi_paras.end_r = ((mpi_paras.world_rank + 1)*(arguments.N - 1))/mpi_paras.world_size + 1;
	mpi_paras.N_r = mpi_paras.end_r - mpi_paras.start_r;

	if (mpi_paras.world_size == 1)
	{
		// einzelner Prozess
		mpi_paras.start_r = 1;
		mpi_paras.end_r = arguments.N;
		mpi_paras.N_r = mpi_paras.end_r - mpi_paras.start_r;

		allocateMatrices(&arguments, &mpi_paras);
		initMatrices(&arguments, &options, &mpi_paras);
		gettimeofday(&start_time, NULL);
		calculate(&arguments, &results, &options);
		gettimeofday(&comp_time, NULL);

		
		double gesamtspeicher = (arguments.N + 1) * (mpi_paras.N_r + 2) * sizeof(double) * arguments.num_matrices / 1024.0 / 1024.0;
		displayStatistics(&results, &options, gesamtspeicher);
		displayMatrix(&arguments, &results, &options, &mpi_paras);
		MPI_Barrier(MPI_COMM_WORLD);
		freeMatrices(&arguments);
	}
	else
	{
		// Matrizen initialisieren
		allocateMatrices(&arguments, &mpi_paras);
		initMatrices(&arguments, &options, &mpi_paras);

		gettimeofday(&start_time, NULL);
		if (options.method == METH_GAUSS_SEIDEL)
		{
			// parallel mit Gauss-Seidel
			// Barrieren für korrekte Zeitmessung
			MPI_Barrier(MPI_COMM_WORLD);
			// Funktion je nach Verfahren
			if (options.termination == TERM_ITER)
			{
				calcGaussSeidelMPI_iter(&arguments, &results, &options, &mpi_paras);
			}
			else
			{
				calcGaussSeidelMPI_prec(&arguments, &results, &options, &mpi_paras);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		else
		{
			// parallel mit Jacobi
			// Barrieren für korrekte Zeitmessung
			MPI_Barrier(MPI_COMM_WORLD);

			calcJacobiMPI(&arguments, &results, &options, &mpi_paras);

			MPI_Barrier(MPI_COMM_WORLD);
		}
		gettimeofday(&comp_time, NULL);

		// Gesamtspeicher ist Summe der einzelnen Speicher
		double speicher_r = (arguments.N + 1) * (mpi_paras.N_r + 2) * sizeof(double) * arguments.num_matrices / 1024.0 / 1024.0;
		double gesamtspeicher;
		MPI_Allreduce(&speicher_r, &gesamtspeicher, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Residuum ist maximale Fehlernorm aller Prozesse
		double fehler_r = results.stat_precision;
		double fehler_all = 0;
		MPI_Allreduce(&fehler_r, &fehler_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		results.stat_precision = fehler_all;

		// Matrix in richtiger Reihenfolge ausgeben
		if (mpi_paras.world_rank == 0)
		{
			displayStatistics(&results, &options, gesamtspeicher);
			displayMatrix(&arguments, &results, &options, &mpi_paras);
			MPI_Ssend(&signal, 1, MPI_INT, mpi_paras.world_rank + 1, 0, MPI_COMM_WORLD);
		}
		if (mpi_paras.world_rank < mpi_paras.world_size && mpi_paras.world_rank > 0)
		{
			//printf("%d: vor rec", mpi_paras.world_rank);
			MPI_Recv(&signal, 1, MPI_INT, mpi_paras.world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//printf("%d: zwischen rec und disp", mpi_paras.world_rank);
			displayMatrix(&arguments, &results, &options, &mpi_paras);
			if (mpi_paras.world_rank < mpi_paras.world_size - 1)
			{
				MPI_Ssend(&signal, 1, MPI_INT, mpi_paras.world_rank + 1, 0, MPI_COMM_WORLD);
			}
		}
			MPI_Barrier(MPI_COMM_WORLD);
			freeMatrices(&arguments);
	}
	
	// finalize MPI Environment
	MPI_Finalize();
	return 0;
}
