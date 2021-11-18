#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    // init MPI Environment
    MPI_Init(&argc, &argv);

    // get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    
    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];
    // declare message feedback
    int message = 0;

    gettimeofday(&tv, NULL);
    gethostname(hostname, 30);

    time = tv.tv_sec;
    micro_sec = tv.tv_usec;
    
    strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
    // testen senden oder empfangen
    if (world_rank < world_size - 1)
    {   
        // output setzen
        snprintf(output, 80, "[%d] %s : %s.%d\n",world_rank, hostname, time_string, (int)micro_sec);
        // output senden
        MPI_Send(&output, 80, MPI_CHAR, world_size - 1, 0, MPI_COMM_WORLD);
        MPI_Send(&micro_sec, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD);
        // message Feedback empfangen
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Recv(&message, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // wenn message == 1 dann erfolgreich
    }   
    else
    {
        // Fall nur ein Prozess, also einfach nur printen einmal
        if(world_size == 1)
        {
            printf("[%d] %s : %s.%d\n",world_rank, hostname, time_string, (int)micro_sec);
            printf("[%d] Kleinseter MS_Anteil: %d \n",world_rank, micro_sec);
            printf("[%d] Größte Differenz: %d \n", world_rank, 0);
        }
        else
        {
            int buffer;
            int min = 999999;
            int max = 0;
            for (int i = 0; i < world_size-1; i++)
            {
                MPI_Recv(&output, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&buffer, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(buffer < min)
                {
                    min = buffer;
                }
                if(buffer > max)
                {
                    max = buffer;
                }
                printf("%s", output);
            }
            printf("[%d] Kleinseter MS_Anteil: %d \n",world_rank, min);
            printf("[%d] Größte Differenz: %d \n", world_rank, max-min);


            //Obwohl wir mehrere Barriers haben, kann die Reihenfolge der Prints nur so eingehalten werden
            //Das ist nicht schön.
            sleep(1);


            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < world_size-1; i++)
            {
                message = 1;
                MPI_Send(&message, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        message = 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // finalize MPI Environment
    if(message == 1)
    {
        printf("[%d] beendet jetzt!\n", world_rank);
    }
    MPI_Finalize();

    // printf("%s\n", output);
    // printf("%d\n", (int)micro_sec);

    return 0;
}
