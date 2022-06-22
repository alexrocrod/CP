// Aula 8 (dia 8/6)
// Alexandre Rodrigues 92993

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))

#define NXMAX 500
#define TOL 1e-7
#define MAXIT 5e5
#define L 1

double f(double x, double y)
{
    return x*y;
}

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
        printf("Introduza numero de pontos {max %d, 0 para sair}: \n",NXMAX);
        scanf(" %d", &nx);
    }
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

    int ndims = 2;
    int dims[2] = {nprocs};
    int periodic[2] = {0};
    MPI_Comm comm1D;
    int newid;
    int nbrbottom, nbrtop;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm1D);

    MPI_Comm_rank(comm1D, &newid);

    MPI_Cart_shift(comm1D, 0, 1, &nbrbottom , &nbrtop);

    printf("myid=%d, newid=%d, bot=%d, top=%d\n", myid, newid, nbrbottom, nbrtop);   

    int firstrow;
    int myrows;

    if (newid == manager_rank)
    {   
        int listfirstrow[nprocs];
        int listmyrows[nprocs];


        int nrows = (int)((double)(ny-2)/(double)nprocs + 0.5f);
        
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

    MPI_Barrier(comm1D);
    printf("newid=%d, firstrow=%d, lastrow=%d\n", newid, firstrow, firstrow+myrows-1);


    double (*Vold)[nx], (*Vnew)[nx], (*myf)[nx];
    Vold = calloc(myrows + 2, sizeof(*Vold));
    Vnew = calloc(myrows + 2, sizeof(*Vnew));
    myf = calloc(myrows + 2, sizeof(*myf));

    double h = ((double)2 * L) / ((double) nx - 1);
    for (int j = 1; j < nx -1 ; j++)
    {
        for (int i = 1; i < myrows + 1; i++)
        {
            myf[i][j] = f(-L + j * h, -L + (firstrow + i - 1) * h);
        }
        
    }

    // initialize to zeros (but calloc already does it)
    if (newid == manager_rank){
        for (int j = 0; j < nx; j++)
        {
            Vnew[0][j] = 0.;
        }
    }

    // initialize to zeros (but calloc already does it)
    if (newid == nprocs - 1){
        for (int j = 0; j < nx; j++)
        {
            Vnew[myrows+1][j] = 0.;
        }
    }

    for (int i = 1; i <= myrows; i++)
    {
        Vnew[i][0] = 0.;
        Vnew[i][nx-1] = 0.;
    }

    double tm1 = MPI_Wtime();

    for (int iter = 0; iter < MAXIT; iter++)
    {
        double sums[2] = {0.0,0.0};
        double global_sums[2];

        for (int j = 1; j < nx -1 ; j++)
        {
            for (int i = 1; i < myrows + 1; i++)
            {
                Vnew[i][j] = (Vold[i+1][j] + Vold[i-1][j] + Vold[i][j+1] + Vold[i][j-1]  + h * h  * myf[i][j]) / 4.0;
                sums[0] += (Vnew[i][j] - Vold[i][j]) * (Vnew[i][j] - Vold[i][j]);
                sums[1] += Vnew[i][j] * Vnew[i][j];
            }
            
        }

        MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm1D);
        
        if (sqrt(global_sums[0]/global_sums[1]) < TOL)
        {
            if (newid == manager_rank)
            {
                printf("calculo durou %f segundos\n", MPI_Wtime() - tm1);
                double (*V)[nx];
                V = calloc(ny, sizeof(*V));

                for (int i = 0; i < myrows +1; i++)
                {
                    for (int j = 0; j < nx; j++)
                    {
                        V[i][j] = Vnew[i][j];
                    }
                    
                }

                for (int i = myrows + 1; i < ny; i++)
                {
                   MPI_Recv(V[i], nx, MPI_DOUBLE, MPI_ANY_SOURCE, i, comm1D, MPI_STATUS_IGNORE);
                    
                }

                tm1 = MPI_Wtime();

                FILE *pfile;
                pfile = fopen("results_2022.dat","w");

                for (int i = 0; i < ny; i++)
                {
                    for (int j = 0; i < nx; j++)
                    {
                        fprintf(pfile, "%.5f ", V[i][j]);
                    }
                    fprintf(pfile,"\n");
                }
                

                printf("escrita durou %f segundos\n", MPI_Wtime() - tm1);

                free(V);
            }
            else
            {
                for (int i = 1; i < myrows + 1; i++)
                {
                   MPI_Send(Vnew[i], nx, MPI_DOUBLE, 0, firstrow + i - 1, comm1D); 
                } 
                if (newid == nprocs -1)
                {
                    MPI_Send(Vnew[myrows+1], nx, MPI_DOUBLE, 0, ny-1, comm1D);
                }
            }

            break;
        }

        // comunicações sentido ascendente
        MPI_Sendrecv(Vnew[myrows], nx, MPI_DOUBLE, nbrtop, 0, Vnew[0] , nx, MPI_DOUBLE, nbrbottom, 0, comm1D, MPI_STATUS_IGNORE);

        // comunicações sentido descendente
        MPI_Sendrecv(Vnew[1], nx, MPI_DOUBLE, nbrbottom, 1, Vnew[myrows+1] , nx, MPI_DOUBLE, nbrtop, 1, comm1D, MPI_STATUS_IGNORE);
        
        for (int i = 0; i < myrows +2; i++)
        {
            for (int j = 0; j < nx; j++)
            {
                Vold[i][j] = Vnew[i][j];
            }
            
        }
        
    }

    free(Vold);
    free(Vnew);
    free(myf);

    


    MPI_Finalize();

    return 0;
}