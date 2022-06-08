// Aula 8 (dia 8/6)
// Alexandre Rodrigues 92993

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
// #include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))

#define NXMAX 500
#define TOL 1e-7
#define MAXIT 5e5
#define L 1

int main(int argc, char *argv[])
{
    int nprocs;
    int myid; 
    int nx, ny;

    int manager_rank = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if (myid == manager_rank)
    {
        printf("Introduza numero de pontos {max %d, 0 para sair}: ",NXMAX);
        scanf(" %d", &nx);
    }
    MPI_Bcast(&nx, 1, MPI_INT, myid, MPI_COMM_WORLD);
    ny = nx;

    if (nx == 0)
    {
        MPI_Finalize();
        return 0;
    }

    if (nx > NXMAX)
    {
        MPI_Finalize();
        return 1;
    }

    int ndims = 1;
    int dims[1] = {nprocs};
    int periodic[1] = {0};
    MPI_Comm comm1D;
    int newid;
    int nbrbottom, nbrtop;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm1D);

    MPI_Comm_rank(comm1D, &newid);

    MPI_Cart_shift(comm1D, 0, 1, &nbrbottom , &nbrtop);

    // printf("myid=%d, newid=%d, bot=%d, top=%d\n", myid, newid, nbrbottom, nbrtop);   

    int firstrow;
    int myrows;

    if (newid == manager_rank)
    {   
        int listfirstrow[nprocs];
        int listmyrows[nprocs];


        int nrows = (int)((double)(ny-2)/(double)+0.5f);
        
        for (int i = 0; i < nprocs; i++)
        {
            listfirstrow[i] = 1 + i *  nrows;
            listmyrows[i] = nrows;
        }

        // altera o numero de linhas do ultimo
        listmyrows[nprocs-1] = ny - 2 - (nprocs - 1) * nrows;

        MPI_Scatter(listfirstrow, 1, MPI_INT, &firstrow, 1, MPI_INT, newid, comm1D);

        MPI_Scatter(listmyrows, 1, MPI_INT, &myrows, 1, MPI_INT, newid, comm1D);
        printf("\n");

    }
    else
    {
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstrow, 1, MPI_INT, manager_rank, comm1D);

        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &myrows, 1, MPI_INT, manager_rank, comm1D);
    }

    printf("newid=%d, firstrow=%d, lastrow=%d\n", newid, firstrow, firstrow+myrows-1);

    MPI_Finalize();

    return 0;
}