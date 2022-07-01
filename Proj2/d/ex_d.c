//
//  Jacobi_1d_v1.c
//  
//
//  Created by Rui Costa on 04/06/2021.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXITER 500000
#define NXMAX 500
#define L 1.0
#define TOL 1e-6

double f(double x, double y){
    //return 2.0 - x*x - 10.0*y + 50.0*x*y;
    return 7*sin(2*M_PI*x)*cos(3*M_PI*x)*sin(2*M_PI*y)*cos(3*M_PI*y);
}

int main(int argc, char *argv[]){

    int i, j;
    int myid, numprocs;
    int nx, ny;
    int ndims = 2;
    int dims[2];
    int periodic[2] = {1,1};
    int newid, nbrtop, nbrbottom, nbrleft, nbrright;
    int mytoprow, myrows, myleftcol, mycols;
    int Jaciter;
    double h;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if (myid == 0) {
        printf("Introduza numero de pontos (0 para sair, maximo %d pontos): \n", NXMAX);
        scanf(" %d",&nx);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if  (nx <= 0 || nx > NXMAX) {
        MPI_Finalize();
        return 1;
    }
    
    ny = nx;
    
    dims[0] = (int) numprocs / 2;
    dims[1] = 2;
    
    numprocs = dims[0]*dims[1];
    
    MPI_Comm comm2d;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &comm2d);
    
    if (comm2d == MPI_COMM_NULL) {
        MPI_Finalize();
        return 0;
    }
    
    MPI_Comm_rank(comm2d, &newid);
    MPI_Cart_shift(comm2d, 0, 1, &nbrtop, &nbrbottom);
    MPI_Cart_shift(comm2d, 1, 1, &nbrleft, &nbrright);
    
    printf("myid %d, newid %d, nbrtop %d, nbrbottom %d, nbrleft %d, nbrright %d\n", myid, newid, nbrtop, nbrbottom, nbrleft, nbrright);
    MPI_Barrier(comm2d);
  
    // calcular a decomposicao do dominio
    if (newid == 0) {
        int numrows = (int)(((double)2*(nx-2)/numprocs) + 0.5);
        int remaining = nx;
        int listtoprows[numprocs], listnumrows[numprocs];
        for (i=0; i<numprocs-2; i+=2) {
            listnumrows[i] = numrows + 1;
            listnumrows[i+1] = numrows + 1;
            listtoprows[i] = nx - remaining;
            listtoprows[i+1] = nx - remaining;
            remaining += -numrows;
            
        }
        listnumrows[numprocs-2] = remaining-1;
        listnumrows[numprocs-1] = remaining-1;
        listtoprows[numprocs-2] = nx - remaining + 1;
        listtoprows[numprocs-1] = nx - remaining + 1;
        
        int listleftcols[numprocs], listnumcols[numprocs];
        for (i=0; i<numprocs; i+=2) {
            listleftcols[i] = 0;
            listnumcols[i] = (ny-2)/2 + 1;
            listleftcols[i+1] = (ny-2)/2 + 1;
            listnumcols[i+1] = ny - listnumcols[i];
           // printf("%d %d mylc %d myrc %d \n",numprocs,i , listleftcols[i], listleftcols[i]+listnumcols[i]-1);
           // printf("%d %d mylc %d myrc %d \n",numprocs,i+1 , listleftcols[i+1], listleftcols[i+1]+listnumcols[i+1]-1);
            
        }
        
        MPI_Scatter(listnumrows, 1, MPI_INT, &myrows, 1, MPI_INT, newid, comm2d);
        MPI_Scatter(listtoprows, 1, MPI_INT, &mytoprow, 1, MPI_INT, newid, comm2d);
        
        MPI_Scatter(listnumcols, 1, MPI_INT, &mycols, 1, MPI_INT, newid, comm2d);
        MPI_Scatter(listleftcols, 1, MPI_INT, &myleftcol, 1, MPI_INT, newid, comm2d);
        
    } else {
        MPI_Scatter(NULL, 1, MPI_INT, &myrows, 1, MPI_INT, 0, comm2d);
        MPI_Scatter(NULL, 1, MPI_INT, &mytoprow, 1, MPI_INT, 0, comm2d);
        
        MPI_Scatter(NULL, 1, MPI_INT, &mycols, 1, MPI_INT, 0, comm2d);
        MPI_Scatter(NULL, 1, MPI_INT, &myleftcol, 1, MPI_INT, 0, comm2d);
    }
    
    printf("%d %d   mytr %d  mybr %d  mylc %d myrc %d,  mycols %d myrows %d \n", myid, newid, mytoprow, mytoprow+myrows-1, myleftcol, myleftcol + mycols -1, mycols, myrows);
    
    // Alocar matrizes em memoria contigua myVold, myVnew e myf, de dimensoes (myrows+2) X ny;
    double (*myVold)[mycols+2], (*myVnew)[mycols+2], (*myf)[mycols+2];
    myVold = calloc(myrows+2, sizeof(*myVold));
    myVnew = calloc(myrows+2, sizeof(*myVnew));
    myf = calloc(myrows+2, sizeof(*myf));
    
    h = 2.0 * L / (nx);
    
    // Calcular myf:
    for (i=1; i<myrows+1; i++){
        for (j=1; j<mycols+1; j++){
            myf[i][j] = f(-L + h*(myleftcol+j-1), L - h*(mytoprow+i-1));
        }
    }
    
    
    // Criar novo datatype para passar colunas
    MPI_Datatype column;
    MPI_Type_vector(myrows+2, 1, mycols+2, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    double tm1 = MPI_Wtime();
    
    // Iterar Jacobi:
    for (Jaciter = 0; Jaciter < MAXITER; Jaciter++) {


        // Calcular Vnew e obter somas das diferenças (Vnew-Vold)^2 e de Vnew^2
        double sums[2];
        sums[0] = 0.0;
        sums[1] = 0.0;

        // Calcular pares
        for (i=1; i<myrows+1; i++) {
            for (j=1; j<mycols+1; j+=2) {
                    myVnew[i][j] = 0.25 * (myVnew[i-1][j] + myVnew[i][j-1] + myVnew[i][j+1] + myVnew[i+1][j] - h * h * myf[i][j]);
                    sums[0] += (myVnew[i][j]-myVold[i][j])*(myVnew[i][j]-myVold[i][j]);
                    sums[1] += myVnew[i][j]*myVnew[i][j];
                    // printf("i %d, j %d, myVnew[i][j] %f \n", i, j, myVnew[i][j]);
            }
        }

        // Comunicar aos vizinhos
        MPI_Sendrecv(&myVnew[1][1], mycols, MPI_DOUBLE, nbrtop, 0, &myVnew[myrows+1][1], mycols, MPI_DOUBLE, nbrbottom, 0, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[myrows][1], mycols, MPI_DOUBLE, nbrbottom, 1, &myVnew[0][1], mycols, MPI_DOUBLE, nbrtop, 1, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[1][1], 1, column, nbrleft, 2, &myVnew[1][mycols+1], 1, column, nbrright, 2, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[1][mycols], 1, column, nbrright, 3, &myVnew[1][0], 1, column, nbrleft, 3, comm2d, MPI_STATUS_IGNORE);

       
        // Calcular ímpares
        for (i=1; i<myrows+1; i++) {
            for (j=2; j<mycols+1; j+=2) {
                    myVnew[i][j] = 0.25 * (myVnew[i-1][j] + myVnew[i][j-1] + myVnew[i][j+1] + myVnew[i+1][j] - h * h * myf[i][j]);
                    sums[0] += (myVnew[i][j]-myVold[i][j])*(myVnew[i][j]-myVold[i][j]);
                    sums[1] += myVnew[i][j]*myVnew[i][j];
            }
        }

        // Comunicar aos vizinhos 
        MPI_Sendrecv(&myVnew[1][1], mycols, MPI_DOUBLE, nbrtop, 0, &myVnew[myrows+1][1], mycols, MPI_DOUBLE, nbrbottom, 0, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[myrows][1], mycols, MPI_DOUBLE, nbrbottom, 1, &myVnew[0][1], mycols, MPI_DOUBLE, nbrtop, 1, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[1][1], 1, column, nbrleft, 2, &myVnew[1][mycols+1], 1, column, nbrright, 2, comm2d, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&myVnew[1][mycols], 1, column, nbrright, 3, &myVnew[1][0], 1, column, nbrleft, 3, comm2d, MPI_STATUS_IGNORE);


        // Agregar somas e testar criterio de paragem
        double sums_global[2];
        MPI_Allreduce(&sums, &sums_global, 2, MPI_DOUBLE, MPI_SUM, comm2d);
        
        printf("(%d) iter %d  err = %e \n", newid, Jaciter, sums_global[0]/sums_global[1]);
        
        if (sqrt(sums_global[0]/sums_global[1]) < TOL) {

            
            if (newid == 0) {
                printf("Numero de iterações: %d\n", Jaciter);
                printf("%f seg para calcular.\n", MPI_Wtime()-tm1);
                tm1 = MPI_Wtime();
            }
            
            // Criar novo datatype para definir a file view de cada processo (cada processo vê apenas a parte pela qual é responsável na matriz global)
            int globalsize[] = {nx, ny};
            int localsize[] = {myrows, mycols};
            int start_inds[] = {mytoprow, myleftcol};
    
            
            //printf("myid %d localsize[0] %d localsize[1] %d start_inds[0] %d start_inds[1] %d \n",myid,localsize[0],localsize[1],start_inds[0],start_inds[1]);
            MPI_Datatype submatrix;
            MPI_Type_create_subarray(2, globalsize, localsize, start_inds, MPI_ORDER_C, MPI_DOUBLE, &submatrix);
            MPI_Type_commit(&submatrix);
            
            MPI_File pf;
            MPI_File_open(comm2d, "results_d_C.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &pf);
            MPI_File_set_view(pf, 0, MPI_DOUBLE, submatrix, "native", MPI_INFO_NULL);
            
            
            // Criar novo datatype para enviar apenas pontos locais pelos quais é responsável (não queremos enviar pontos fantasma)
            int matrixsize[] = {myrows+2,mycols+2};
            start_inds[0] = 1;
            start_inds[1] = newid % 2;
            if (newid == 0 || newid == 1) {
                start_inds[0]--;
            }

            // printf("myid %d, matrixsize[0] %d, matrixsize[1] %d, start_inds[0] %d, start_inds[1] %d \n", myid, matrixsize[0], matrixsize[1], start_inds[0], start_inds[1]);

            
            MPI_Datatype localmatrix;
            MPI_Type_create_subarray(2, matrixsize, localsize, start_inds, MPI_ORDER_C, MPI_DOUBLE, &localmatrix);
            MPI_Type_commit(&localmatrix);
            
            // Escrever ficheiro binário
            MPI_File_write_all(pf, myVnew, 1, localmatrix, MPI_STATUS_IGNORE);
            MPI_File_close(&pf);
            
            
            MPI_Type_free(&submatrix);
            MPI_Type_free(&localmatrix);
            
            if (newid == 0) {
                printf("%f seg para escrever ficheiro.\n", MPI_Wtime()-tm1);
            }
    
            break;
        }

        // Na direcção 'vertical':
        // Este sendrecv envia para cima e recebe de baixo
        MPI_Sendrecv(myVnew[1], mycols+2, MPI_DOUBLE, nbrtop, 0, myVnew[myrows+1], mycols+2, MPI_DOUBLE, nbrbottom, 0, comm2d, MPI_STATUS_IGNORE);
        // Este sendrecv envia para baixo e recebe de cima
        MPI_Sendrecv(myVnew[myrows], mycols+2, MPI_DOUBLE, nbrbottom, 1, myVnew[0], mycols+2, MPI_DOUBLE, nbrtop, 1, comm2d, MPI_STATUS_IGNORE);

        // Na direcção 'horizontal':
        // Este sendrecv envia para a esquerda e recebe da direita
        MPI_Sendrecv(&(myVnew[0][1]), 1, column, nbrleft, 2, &(myVnew[0][mycols+1]), 1, column, nbrright, 2, comm2d, MPI_STATUS_IGNORE);
        // Este sendrecv envia para a direita e recebe da esquerda
        MPI_Sendrecv(&(myVnew[0][mycols]), 1, column, nbrright, 3, &(myVnew[0][0]), 1, column, nbrleft, 3, comm2d, MPI_STATUS_IGNORE);
        
        // Vold = Vnew
        for (i=0; i<myrows+2; i++) {
            for (j=0; j<mycols+2; j++) {
                myVold[i][j] = myVnew[i][j];
            }
        }
        
    }
    
    free(myf);
    free(myVnew);
    free(myVold);
    MPI_Type_free(&column);
    
    
    MPI_Finalize();
    return 0;
}