#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>

#define min(a,b) (((a)<(b))?(a):(b))

// int min(int a, int b)
// {
//     if (b > a) 
//         b = a;
//     return b;
// }

int main(int argc, char *argv[])
{

    int manager_rank = 0;
    int nrows = 10;
    int ncols = 10;
    // printf("Matrix %dx%d\n",nrows,ncols);

    int b[ncols];
    
    int my_rank;
    int term_tag = nrows;
    int sender;
    int next_row  = 0;

    MPI_Init(&argc,&argv);
    int nprocs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int maxIdx = min(nprocs - 1, nrows);
    // printf("nprocs: %d\n",nprocs);
    // printf("maxIdx: %d\n",maxIdx);

    MPI_Status status;
    
    if (my_rank == manager_rank)
    {

        int row_index;
        int ans;

        int A[nrows][ncols];
        int c[nrows];

        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; i < ncols; j++)
            {
                A[i][j] = ncols*i+j;
            }
            
        }
        for (int i = 0; i < ncols; i++)
        {
            b[i] = i;
        }
        
        MPI_Bcast(b, ncols, MPI_INT, manager_rank, MPI_COMM_WORLD);

        // pelos ranks dos processos
        for (int i = 1; i <=maxIdx; i++)
        {
            MPI_Send(A[next_row], ncols, MPI_INT, i, next_row, MPI_COMM_WORLD);
            next_row++;
        }


        for (int i = 0; i < nrows; i++)
        {
            MPI_Recv(&ans, 1 , MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row_index = status.MPI_TAG;
            // printf("row_idx: %d\n",row_index);
            sender = status.MPI_SOURCE;
            c[row_index] = ans;

            if (next_row < nrows)
            {
                // printf("ola %d\n", next_row);
                MPI_Send(A[next_row], ncols, MPI_INT, sender, next_row, MPI_COMM_WORLD);
                next_row++;
            }
            else
            {
                MPI_Send(MPI_BOTTOM, 0, MPI_INT, sender, term_tag, MPI_COMM_WORLD);
                // printf("ending: %d\n",sender);
            }
        }

        // printf("Array: c=[ ");
        // for (int i = 0; i < nrows; i++)
        // {
        //     printf("%d, ",c[i]);
        // }
        // printf("\b\b]\n");

        // printf("nrows: %d\n",nrows);
        printf("Array c = \n");
		for (int i=0; i<nrows; i++) 
        {
			printf("%d \n",c[i]);
		}
		printf("\n");
              
    } 
    else 
    {
        MPI_Bcast(b, ncols, MPI_INT, manager_rank, MPI_COMM_WORLD);

        int row[ncols];
        int row_index;
        int ans;
        
        while (1)
        {
            
            MPI_Recv(row, ncols, MPI_INT, manager_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row_index = status.MPI_TAG;

            if (row_index == term_tag)
                break;

            ans = 0;

            for (int i = 0; i < ncols; i++)
            {
                ans += row[i]*b[i];
                
            }
            // printf("ans(%d)=%d\n",my_rank,ans);

            MPI_Send(&ans, 1, MPI_INT, manager_rank, row_index, MPI_COMM_WORLD);  
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == manager_rank)
        printf("All processes have finished.\n");

    MPI_Finalize();

    return 0;
}
