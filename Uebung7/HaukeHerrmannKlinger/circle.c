#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int*
init (int N, int rank)
{	
	//Größe + ein Puffer allokieren
	int* buf = (int*)malloc(sizeof(int) * (N + 1));

	//jeder Prozess eigenen Zufallsseed
	srand(time(NULL) + rank);

	for (int i = 0; i < N; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	return buf;
}

int
circle (int* buf, int size, int rank, int *N_lokal, int number_to_check)
{
	int iterationen = -1;
	int fertig = 0;
	int new_N_lokal;
	//Puffer für Arrayabschnittinhalt
	int puffer[*N_lokal];

	while (fertig == 0)
	{
		//Größe der einzelnen Array Abschnitte übermitteln
		if (rank == 0)
		{
			MPI_Ssend(N_lokal, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&new_N_lokal, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Recv(&new_N_lokal, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank < size - 1)
			{
				MPI_Ssend(N_lokal, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Ssend(N_lokal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}

		//Inhalte der Arrays übergeben
		if (rank == 0)
		{
			MPI_Ssend(buf, *N_lokal, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(buf, new_N_lokal, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			//alten Array Inhalt zwischenspeichern
			for (int i = 0; i<*N_lokal; i++)
			{
				puffer[i]=buf[i];
			}
			MPI_Recv(buf, new_N_lokal, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank < size - 1)
			{
				MPI_Ssend(puffer, *N_lokal, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Ssend(puffer, *N_lokal, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}

		//Arraygröße aktualisieren
		*N_lokal = new_N_lokal;

		//Überprüfen, ob Iteration abgebrochen werden muss
		if (rank == size - 1)
		{
			if (buf[0] == number_to_check)
			{
				fertig = 1;
			}
			for (int i = 0; i < size - 1; i++)
			{
				MPI_Ssend(&fertig, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
		}
		else
		{
			MPI_Recv(&fertig, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		iterationen++;
	}
	return iterationen;
}

int
main (int argc, char** argv)
{
	int N;
	int rank;
	int* buf;
    int ret;
	int N_lokal;
	int start_lokal;
	int end_lokal;
	int number_to_check;
	int size;
	int signal = 0;

    // init MPI Environment
    MPI_Init(NULL, NULL);

    // get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get rank of process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}
	
	// gesamtes Array length
	N = atoi(argv[1]);

	// Position im gesamten Array
	start_lokal = (N*rank)/size;
	end_lokal = (N*(rank + 1))/size;
	N_lokal = end_lokal - start_lokal;
	printf("%d: start ist %d, ende ist %d\n", rank, start_lokal, end_lokal);

	// Teilarray initialisieren
	buf = init(N_lokal, rank);

	//Falls nur ein Prozess
	if (size == 1)
	{
		printf("\nBEFORE\n");
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
		printf("\nIterationen:%d, Abbruchwert:%d\n\n",0, buf[0]);
		printf("\nAFTER\n");
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
		return EXIT_SUCCESS;
	}

	if (rank == 0)
	{
		number_to_check = buf[0];
		MPI_Ssend(&number_to_check, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
		//printf("%d\n", number_to_check);
	}
	if (rank == size - 1)
	{
		MPI_Recv(&number_to_check, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("%d\n", number_to_check);
	}


	//Erster Prozess printet zuerst
	if (rank == 0)
	{
		printf("\nBEFORE\n");
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
	}
	else 
	{
		//Prozesse printen ihren Teil nach dem Vorgänger
		MPI_Recv(&signal, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
	}

	//Alle Prozesse geben ein Signal an den Nachfolger
	if (rank < size - 1)
	{
		MPI_Ssend(&signal, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	ret = circle(buf, size, rank, &N_lokal, number_to_check);

	MPI_Barrier(MPI_COMM_WORLD);
	sleep(0.1);

	//Ausgabe wie am Anfang nur mit "AFTER"
	if (rank == 0)
	{
		printf("\nIterationen:%d, Abbruchwert:%d\n\n",ret, number_to_check);
		printf("\nAFTER\n");
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
	}
	else 
	{
		MPI_Recv(&signal, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = 0; i < N_lokal; i++)
		{
			printf("rank %d: %d\n", rank, buf[i]);
		}
	}
	if (rank < size - 1)
	{
		MPI_Ssend(&signal, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return EXIT_SUCCESS;
}
