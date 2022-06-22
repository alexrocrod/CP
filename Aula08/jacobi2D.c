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

    int nprocs_col = (int) nprocs/2;

    int ndims = 2;
    int dims[2] = {nprocs_col, 2};
    int periodic[2] = {0,0};
    MPI_Comm comm2D;
    int newid;
    int nbrbottom, nbrtop, nbrleft, nbrright;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm2D);

    if (comm2D == MPI_COMM_NULL) // processo extra é terminado (se nprocs é impar)
    {
        MPI_Finalize();
        return 0;
    }

    nprocs = nprocs_col * dims[1];

    MPI_Comm_rank(comm2D, &newid);

    MPI_Cart_shift(comm2D, 0, 1, &nbrbottom , &nbrtop);
    MPI_Cart_shift(comm2D, 1, 1, &nbrleft , &nbrright);

    printf("myid=%d, newid=%d, bot=%d, top=%d, left=%d, right=%d\n", myid, newid, nbrbottom, nbrtop, nbrleft, nbrright);   

    int firstrow, firstcol;
    int myrows, mycols;

    if (newid == manager_rank)
    {   
        int listfirstrow[nprocs];
        int listmyrows[nprocs];

        int listfirstcol[nprocs];
        int listmycols[nprocs];

        int nrows = (int)((double)(ny-2)/(double)nprocs_col + 0.5);
        
        // Linhas
        for (int i = 0; i < nprocs_col; i++)
        {
            listfirstrow[2*i] = 1 + i *  nrows;
            listmyrows[2*i] = nrows;
            listfirstrow[2*i+1] = 1 + i *  nrows;
            listmyrows[2*i+1] = nrows;
        }
        // altera o numero de linhas do penultimo e do ultimo
        listmyrows[nprocs-2] = ny - 2 - (nprocs_col - 1) * nrows;
        listmyrows[nprocs-1] = ny - 2 - (nprocs_col - 1) * nrows;

        // Colunas
        int ncols_temp = (int)((nx-2)/2);
        for (int i = 0; i < nprocs_col; i++)
        {
            listfirstcol[2*i] = 1;
            listmycols[2*i] = ncols_temp;
            listfirstcol[2*i+1] = ncols_temp + 1;
            listmycols[2*i+1] = nx - 2 - ncols_temp;
        }

        MPI_Scatter(listfirstrow, 1, MPI_INT, &firstrow, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listmyrows, 1, MPI_INT, &myrows, 1, MPI_INT, newid, comm2D);

        MPI_Scatter(listfirstcol, 1, MPI_INT, &firstcol, 1, MPI_INT, newid, comm2D);
        MPI_Scatter(listmycols, 1, MPI_INT, &mycols, 1, MPI_INT, newid, comm2D);

        printf("\n");
    }
    else
    {
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstrow, 1, MPI_INT, manager_rank, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &myrows, 1, MPI_INT, manager_rank, comm2D);

        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &firstcol, 1, MPI_INT, manager_rank, comm2D);
        MPI_Scatter(MPI_BOTTOM, 1, MPI_INT, &mycols, 1, MPI_INT, manager_rank, comm2D);
    }

    MPI_Barrier(comm2D);
    printf("newid=%d, firstrow=%d, lastrow=%d, firstcol=%d, lastcol=%d\n", newid, firstrow, firstrow+myrows-1, firstcol, firstcol+mycols-1);


    double (*Vold)[mycols+2], (*Vnew)[mycols+2], (*myf)[mycols+2];
    Vold = calloc(myrows + 2, sizeof(*Vold));
    Vnew = calloc(myrows + 2, sizeof(*Vnew));
    myf = calloc(myrows + 2, sizeof(*myf));

    double h = ((double)2 * L) / ((double) nx - 1);

    for (int j = 1; j < mycols + 1 ; j++)
    {
        for (int i = 1; i < myrows + 1; i++)
        {
            myf[i][j] = f(-L + (firstcol + j - 1) * h, -L + (firstrow + i - 1) * h);
        }
        
    }

    // initialize to zeros (but calloc already does it)
    if (newid == manager_rank || newid == 1){
        for (int j = 0; j < mycols+2; j++)
        {
            Vnew[0][j] = 0.;
            Vold[0][j] = 0.;
        }
    }

    // initialize to zeros (but calloc already does it)
    if (newid == nprocs - 1 || newid == nprocs - 1){
        for (int j = 0; j < mycols + 2; j++)
        {
            Vnew[myrows+1][j] = 0.;
            Vold[myrows+1][j] = 0.;
        }
    }

    if (newid % 2 == 0)
    {
        for (int i = 1; i < myrows + 1; i++)
        {
            Vnew[i][0] = 0.;
            Vold[i][0] = 0.;
        }
    }
    else
    {
       for (int i = 1; i < myrows + 1; i++)
        {
            Vnew[i][mycols+1] = 0.;
            Vold[i][mycols+1] = 0.;
        } 
    }

    MPI_Datatype column;
    MPI_Type_vector(myrows + 2, 1, mycols + 2 , MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    double tm1 = MPI_Wtime();

    for (int iter = 0; iter < MAXIT; iter++)
    {
        double sums[2] = {0.0,0.0};
        double global_sums[2];

        for (int j = 1; j < mycols -1 ; j++)
        {
            for (int i = 1; i < myrows + 1; i++)
            {
                Vnew[i][j] = (Vold[i+1][j] + Vold[i-1][j] + Vold[i][j+1] + Vold[i][j-1]  + h * h  * myf[i][j]) / 4.0;
                sums[0] += (Vnew[i][j] - Vold[i][j]) * (Vnew[i][j] - Vold[i][j]);
                sums[1] += Vnew[i][j] * Vnew[i][j];
            }
            
        }

        MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, comm2D);
        
        if (sqrt(global_sums[0]/global_sums[1]) < TOL)
        {
            if (newid == manager_rank)
            {
                printf("calculo durou %f segundos\n", MPI_Wtime() - tm1);
                printf("%d iteracoes\n", iter);

                // double (*V)[nx];
                // V = calloc(ny, sizeof(*V));

                // for (int i = 0; i < myrows +1; i++)
                // {
                //     for (int j = 0; j < nx; j++)
                //     {
                //         V[i][j] = Vnew[i][j];
                //     }
                    
                // }

                // for (int i = myrows + 1; i < ny; i++)
                // {
                //    MPI_Recv(V[i], nx, MPI_DOUBLE, MPI_ANY_SOURCE, i, comm1D, MPI_STATUS_IGNORE);
                    
                // }

                // tm1 = MPI_Wtime();

                // FILE *pfile;
                // pfile = fopen("results_2022.dat","w");

                // for (int i = 0; i < ny; i++)
                // {
                //     for (int j = 0; i < nx; j++)
                //     {
                //         fprintf(pfile, "%.5f ", V[i][j]);
                //     }
                //     fprintf(pfile,"\n");
                // }
                

                // printf("escrita durou %f segundos\n", MPI_Wtime() - tm1);

                // free(V);
            }
            else
            {
            //     for (int i = 1; i < myrows + 1; i++)
            //     {
            //        MPI_Send(Vnew[i], nx, MPI_DOUBLE, 0, firstrow + i - 1, comm1D); 
            //     } 
            //     if (newid == nprocs -1)
            //     {
            //         MPI_Send(Vnew[myrows+1], nx, MPI_DOUBLE, 0, ny-1, comm1D);
            //     }
            }

            break;
        }

        // comunicações sentido ascendente
        MPI_Sendrecv(Vnew[myrows], mycols+2, MPI_DOUBLE, nbrtop, 0, Vnew[0] , mycols+2, MPI_DOUBLE, nbrbottom, 0, comm2D, MPI_STATUS_IGNORE);

        // comunicações sentido descendente
        MPI_Sendrecv(Vnew[1], mycols+2, MPI_DOUBLE, nbrbottom, 1, Vnew[myrows+1] , mycols+2, MPI_DOUBLE, nbrtop, 1, comm2D, MPI_STATUS_IGNORE);

        // comunicações sentido para direita
        MPI_Sendrecv(&(Vnew[0][mycols]), 1, column, nbrright, 2, &(Vnew[0][0]), 1, column, nbrleft, 2, comm2D, MPI_STATUS_IGNORE);

        // comunicações sentido para esquerda
        MPI_Sendrecv(&(Vnew[0][1]), 1, column, nbrleft, 3, &(Vnew[0][mycols+1]), 1, column, nbrright, 3, comm2D, MPI_STATUS_IGNORE);
        
        for (int i = 0; i < myrows + 2; i++)
        {
            for (int j = 0; j < mycols + 2; j++)
            {
                Vold[i][j] = Vnew[i][j];
            }
            
        }
        
    }

    MPI_Type_free(&column);

    free(Vold);
    free(Vnew);
    free(myf);

    MPI_Finalize();

    return 0;
}