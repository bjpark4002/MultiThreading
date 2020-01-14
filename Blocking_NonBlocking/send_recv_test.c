
/*      Cpt S 411, Introduction to Parallel Computing
 *      School of Electrical Engineering and Computer Science
 *      
 *      Example code
 *      Send receive test:
 *      rank 1 sends to rank 0 (all other ranks sit idle)
 *      For timing use of C gettimeofday() is recommended.
 * */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>

#define BufSize 1


int main(int argc,char *argv[]){

   int rank,p;
   struct timeval t1,t2;


   int number = 1;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&p);

   MPI_Request request, request2;
   MPI_Status status;


   printf("my rank=%d\n",rank);
   printf("Rank=%d: number of processes =%d\n",rank,p);

   assert(p>=2);


   if(rank==1) {
                int x = 10;
                int dest = 0;
                int i = 1;

        //-------------Blocking Send--------------------
                for(i=1 ; i < 2097152 ; i*=2)
                {
                        char message[i];
                        int j = 0;
                        for (j; j<i; j++)
                        {
                                message[j] = 'a';
                        }
                        gettimeofday(&t1,NULL);
                        MPI_Send(&message,i,MPI_BYTE,dest,0,MPI_COMM_WORLD);
                        gettimeofday(&t2,NULL);
                        int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                        printf(" [Blocking Send] %d microsecond, %d bytes \n",tSend, i);
                }


        //--------------Non-Blocking Send---------------


                for(i=1 ; i < 2097152 ; i*=2)
                {
                        char message[i];
                        int j = 0;
                        for (j; j<i; j++)
                        {
                                message[j] = 'a';
                        }
                        gettimeofday(&t1,NULL);
                        MPI_Isend(&message,i,MPI_BYTE,dest,0,MPI_COMM_WORLD,&request);
                        MPI_Wait(&request, &status);
                        gettimeofday(&t2,NULL);
                        int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                        printf(" [Non-Blocking Send] %d microsecond, %d bytes \n",tSend, i);
                }





   }
   else if (rank==0) {
                int y=0;
                int i = 1;


//---------------------Blocking Recv------------------------
                for(i=1; i< 2097152 ; i*=2)
                {
                        char message[i];

                        int j = 0;
                        for (j; j<i; j++)
                        {
                                message[j] = 'a';
                        }
                        gettimeofday(&t1,NULL);
                        MPI_Recv(&message,i,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                        gettimeofday(&t2,NULL);
                        int tRecv = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                        printf("[Blocking Recv] %d microsecond, %d bytes \n",tRecv, i);
}

//--------------------Non-Blocking Recv------------------------


                for(i=1; i< 2097152 ; i*=2)
                {
                        char message[i];

                        int j = 0;
                        for (j; j<i; j++)
                        {
                                message[j] = 'a';
                        }
                        gettimeofday(&t1,NULL);
                        MPI_Recv(&message,i,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
//                        MPI_Wait(&request,&status);

                        gettimeofday(&t2,NULL);
                        int tRecv = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                        printf("[Non-Blocking Recv] %d microsecond, %d bytes \n",tRecv, i);
                }





   }
   MPI_Finalize();
}







