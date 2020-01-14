#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>


/*
program runs with command below
1. mpicc -o p2 main.c  // to compile
2. mpirun -np <number of process> p2 <size of input>  // to check result



To run and see the result from each function, uncomment out function calls in the main.
For Sum
1.myReduceSum() 
2.naiveReduceSum() 
3.MPILibraryReduceSum() 
For Max 
1.myReduceMax() 
2.naiveReduceMax() 
3.MPILibraryReduceMax()


*/



void naiveReduceSum( int*, int, int, int);
void naiveReduceMax( int*, int, int, int);
void reduceWithSameNP(int *, int, int,int);
void myReduceSum(int *, int, int, int); //for sum
void myReduceMax(int *, int, int, int); // for max
void MPILibraryReduceSum(int *,int,int,int);
void MPILibraryReduceMax(int *,int,int,int);//for max

int *generateArray(int);
void dispArray( int *, int);


int lg2(int ); 
int fromBinaryToDecimal( int );
void decToBinary(int, int*,int,int );
int togglingTth( int , int ,int,int);

int main(int argc, char *argv[]){
	int rank,p,n;
	struct timeval t1,t2;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

//	printf("my rank=%d\n",rank);
//	printf("Rank=%d: number of processes =%d\n",rank,p);

	if(argc<1){
		printf("Usage : p2 {array size} \n");
		exit(1);
	}
	n = atoll(argv[1]);
//	printf("Debug : array size = %d\n",n);
	int *array = generateArray(n);
//	dispArray(array,n);
//	reduceWithSameNP(array,rank,n,p );

//	printf("-------------------sum------------------\n");

//sum part
	myReduceSum(array, rank,n ,p);
	MPILibraryReduceSum(array, rank, n, p);
	naiveReduceSum(array, rank, n, p);
//	myReduceSum(array, rank, n, p);

//	myReduceSum(array, rank, n, p);
//	naiveReduceSum(array, rank, n, p);

//	MPILibraryReduceSum(array, rank, n, p);



//	printf("----------------------------------------\n");
//	printf("\n");

//max part
//	myReduceMax(array, rank, n, p);
//	naiveReduceMax(array, rank, n, p);
//	MPILibraryReduceMax(array, rank, n, p);
// 	myReduceMax(array, rank, n, p);
//	naiveReduceMax(array, rank, n, p);




	MPI_Finalize();
}

void MPILibraryReduceSum(int *array, int rank, int n, int p)
{
	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	int localSum =0;
	int portion =n/p;
	int i =0;
	int globalSum=0;
//	printf(" portion == %d\n",portion);

	for(i = rank*portion; i< (rank*portion)+portion; i++)
	{
		localSum += array[i];
	}
	//printf("localSum = %d\n",localSum);
	MPI_Allreduce(&localSum, &globalSum,1 , MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if( rank == (p-1))
	{
		gettimeofday(&t2,NULL);

		int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
		printf(" MPI_Allreduce Sum == %d,  %d microsecond\n",globalSum, tSend);

	}
}


void MPILibraryReduceMax(int *array, int rank, int n, int p)
{

	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
        int localMax =0;
        int portion =n/p;
        int i =0;
        int globalMax=0;
        for(i = rank*portion; i< (rank*portion)+portion; i++)
        {

		if(localMax < array[i]){
			localMax = array[i];
		}
        }
        MPI_Allreduce(&localMax, &globalMax,1 , MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        if( rank == (p-1))
        {

		gettimeofday(&t2,NULL);
                int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                printf(" MPI_Allreduce Max == %d, %d microsecond\n",globalMax,tSend);
        }
}




	// for sum.
void myReduceSum( int * array, int rank , int n , int p)
{

	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
	//printf(" my Reduce called ****\n");
        int localSum = 0;
        int t = 0,i=0 ;
        int partner;
        int bit = lg2(p);
	int portion = n/p;
	int recievedSum =0;
	MPI_Status status;




	for(i = rank*portion ; i< (rank*portion)+portion; i++)
	{
//		printf("??\n");
		localSum += array[i];
	}
	for(t = 0 ;  t<bit; t++)
	{
                partner = togglingTth(t,bit,rank,p);
//		printf(" rank = %d ___ partner = %d\n",rank,partner);
		if(partner < rank){
			MPI_Recv(&recievedSum,1,MPI_INT,partner, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			localSum = localSum + recievedSum;
		}
		else{
			MPI_Send(&localSum,1, MPI_INT, partner, 0, MPI_COMM_WORLD);
		}
        }
	if(rank ==(p-1))
	{
		gettimeofday(&t2,NULL);
                int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
		printf("myReduce Sum == %d, %d microsecond\n",localSum,tSend);
	}
}

void myReduceMax( int * array, int rank , int n , int p)
{

	struct timeval t1,t2;
	gettimeofday(&t1,NULL);
        int localMax = 0;
        int t = 0,i=0 ;
        int partner;
        int bit = lg2(p);
        int portion = n/p;
        int recievedMax =0;
        MPI_Status status;
        for(i = rank*portion ; i< (rank*portion)+portion; i++)
        {
		if (array[i] > localMax){
			localMax = array[i];
        	}
	}
        for( t = 0; t < bit;t++)
        {
			partner = togglingTth(t,bit,rank,p);
            if(partner < rank){
                MPI_Recv(&recievedMax,1,MPI_INT,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            	if(localMax<recievedMax){	
					localMax = recievedMax;
				}
			}
			else{
                        MPI_Send(&localMax,1, MPI_INT, partner, 0, MPI_COMM_WORLD);
                }
               // printf("partner rank == %d || my rank == %d\n\n", partner,rank);
        }
        if(rank ==(p-1))
        {
		gettimeofday(&t2,NULL);
                int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                printf("myReduceMax Max == %d, %d microsecond\n",localMax,tSend);
        }
}

int*  generateArray(int n){
	int *a = (int *)malloc(sizeof(int)*n);
	int i = 0;
	for(i = 0 ; i < n ; i++){
		a[i] = rand() % 10 +1 ;
	}
	return a;
}

void dispArray(int *array, int n){
	int i;
	printf("Debug : --------- array contents ------------\n");
 	for(i=0; i<n; i++){
		printf("                   %d      \n",array[i]);
	}
	printf("---------------------------------------------\n");
}


void reduceWithSameNP( int * array, int rank, int n, int p)
{

	assert(n == p);
	int localSum=array[rank];
	int recievedSum = 0;
	MPI_Status status;

	if(rank > 0 )
	{
		MPI_Recv(&recievedSum, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		localSum += recievedSum;
		if(rank < (p-1) )
		{
			MPI_Send(&localSum, 1 , MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		}
	}

	if(rank == 0)
	{
		MPI_Send(&localSum,1,MPI_INT,rank+1,0,MPI_COMM_WORLD );

	}

	if(rank == p-1)
	{
		printf(" reduceWithSameNP()  sum === %d\n",localSum);
	}
}

void naiveReduceSum( int * array, int rank , int n , int p)
{
	struct timeval t1,t2;



	gettimeofday(&t1,NULL);
	assert(n%p == 0);
//	printf(" p = %d\n",p);
//	printf(" log_2( p ) = %d\n",lg2(8));
	int localSum = 0;
	int recievedSum = 0;
	int portion = n/p;
	int i = 0;
	MPI_Status status;


        for(i = rank*portion ; i < (rank*portion)+portion ; i++){
		localSum += array[i];
        }


	if(rank > 0 ){
		MPI_Recv(&recievedSum, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	//	printf(" rank %d 's recieved Sum = %d\n",rank, recievedSum);
		localSum += recievedSum;
//		for(i = rank*portion ; i < (rank*portion)+portion ; i++){
//			localSum += array[i];
//		}
//		localSum += recievedSum;
		if(rank < (p-1))
		{
//			for(i = rank*portion ; i < (rank*portion)+portion ; i++){
  //                      localSum += array[i];
//                	}
			MPI_Send(&localSum, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		}
	}
	if(rank == 0){
//		for(i = rank*portion ; i < (rank*portion)+portion ; i++){
//			localSum += array[i];
//		}
	//	printf(" rank 0's localSum to send = %d\n",localSum);
		MPI_Send(&localSum, 1 , MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}
	if(rank == p-1){
		gettimeofday(&t2,NULL);
                int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
		printf("naiveReduce  Sum == %d, %d microsecond\n",localSum, tSend);
	}
}

void naiveReduceMax( int * array, int rank , int n , int p)
{
	struct timeval t1, t2;
	gettimeofday(&t1,NULL);
        assert(n%p == 0);
        int localMax = 0;
        int recievedMax = 0;
        int portion = n/p;
        int i = 0;
        MPI_Status status;
        if(rank > 0 ){
                MPI_Recv(&recievedMax, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        //      printf(" rank %d 's recieved Sum = %d\n",rank, recievedSum);
		if( localMax <recievedMax){
			localMax = recievedMax;
		}
                for(i = rank*portion ; i < (rank*portion)+portion ; i++){
			if( array[i] >= localMax){
				localMax = array[i];
			}
                }
                if(rank < (p-1))
                {
                        MPI_Send(&localMax, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
                }
        }
        if(rank == 0){
                for(i = rank*portion ; i < (rank*portion)+portion ; i++){
                    if(array[i]> localMax){
			localMax = array[i];
			}
                }
        //      printf(" rank 0's localSum to send = %d\n",localSum);
                MPI_Send(&localMax, 1 , MPI_INT, rank+1, 0, MPI_COMM_WORLD);

        }
        if(rank == p-1){

		gettimeofday(&t2,NULL);
                int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                printf("naiveReduceMax  Max == %d, %d microSecond\n",localMax, tSend);
        }
}

int togglingTth( int t, int bit, int rank,int p)
{
	int *binaryNum = (int *)malloc(sizeof(int)*bit);
	int init=0;
	for(init = 0; init< bit; init++){
		binaryNum[init] = 0;
	}
//	decToBinary(rank,binaryNum,bit,p);


        int k = bit-1;
        int sum=0;
        int tem=0;

        for (k ; k>=0 ; k--)
        {
                tem = sum+ipow(2,k);
                if(tem <=rank  )
                {
                        binaryNum[k] = 1;
                        sum += ipow(2,k);
                }
                else{
                        binaryNum[k] = 0;
                }
        }

	//flip

	if( binaryNum[t] == 0){
		binaryNum[t] =1;
	}
	else{
		binaryNum[t] = 0;
	}

	sum = 0;
	for(k =0; k< bit; k++)
	{
		if(binaryNum[k]==1){
			sum += ipow(2,k);
		}
	}


	free(binaryNum);
	return sum;

}


int lg2(int input)
{
	int value=0;
	while(input !=1){
		input /=2;
		value++;
	}
	return value;
}

int ipow(int base, int exp)
{

	if(exp == 0)
	{
		return 1;
	}
	int result = 1;
	for(;;)
	{
		if(exp &1)
			result *= base;
		exp >>= 1;
		if(!exp)
			break;
		base *= base;

	}
	return result;

}

int fromBinaryToDecimal( int input)
{
	int decimal,i,a;
	for(i = 0 ; input!=0; ++i)
	{
		a = input%10;
		decimal=(a)*(ipow(2,i))+decimal;
		input=input/10;
	}

	return decimal;
}

void decToBinary(int n , int* binaryNum ,int bit, int p )
{

//	printf("in Dec rank =%d \n",n);

	int i = bit-1;
	int sum=0;
	int tem=0;

	for (i ; i>=0 ; i--)
	{
		tem = sum+ipow(2,i);
		if(tem <=n  )
		{
			binaryNum[i] = 1;
			sum += ipow(2,i);
		}
		else{
			binaryNum[i] = 0;
		}
	}


}