#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>


//matrix functions
void createMatrix(int power, int matrix[2][2]);
void printMatrix(int matrix[2][2], int rank);
void matrixMultiplication(int first[2][2], int second[2][2]);
void matrixSum(int first[2][2], int second[2][2]);
void copyMatrix(int first[2][2], int second[2][2]);
void p_element_pp(int matrix[2][2], int p, int rank, int output[2][2]);
void parallel_serial_matrix(int x_zero, int mOff[2][2], int portion, int rank);

//togle functions
int ipow(int base, int exp);
int togglingTth( int , int ,int,int);
int lg2(int);

//testing
void serial_baseline(int n, int x_zero);
int serial_matrix_helper(int first[2], int second[2][2]);
void serial_matrix(int n, int x_zero);
// these are prime numbers
int A = 523;
int B = 883;
int P = 11801;

int main(int argc, char *argv[]){

	if(argc<3){
                printf(" Usage : p3 <seed> <n> \n");
                exit(1);
        }
	struct timeval t1, t2;
	int rank,p;
	int x_zero = atoll(argv[1]);
	int n = atoll(argv[2]);
	// Step1,  initializing M, M^0, Array Xlocal of size n/p part.
	int M[2][2];
	int M_zero[2][2];
	int M_global[2][2];
	createMatrix(1,M);
	createMatrix(0,M_zero);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// step 2. at each process Pj

	// to exclude dummy time
	int emptyGarbageTime =0 ;
        int garbage=0;
	for(emptyGarbageTime; emptyGarbageTime< 99;emptyGarbageTime++)
	{
		garbage++;
        }



	gettimeofday(&t1,NULL);
	int Mlocal[2][2];
	createMatrix(0, Mlocal); // Mlocal = M^0

	int portion = n/p;
       	int startIndex = rank*portion;
	int localMatrixArray[portion][2][2];
	int offset[2][2] = {{0,0}, {0,0}  };
	int mOff[2][2];
	int i;
	if( rank != 0){  //  0's process's local matrix should be M^0
		for( i = 0 ; i <portion; i++){ // initialize local Matrix array.
			createMatrix(1,localMatrixArray[i]);
			matrixMultiplication(Mlocal, localMatrixArray[i]);
		}
	}
	//Step 4
	p_element_pp(Mlocal, p, rank, mOff);
	parallel_serial_matrix(x_zero, mOff, portion,rank); // this prints out the values in series. each process outputs n/p values.
	gettimeofday(&t2,NULL);

	if(rank ==0){
	        float secs = (t2.tv_sec - t1.tv_sec) * 1000.0f + (t2.tv_usec - t1.tv_usec) / 1000.0f;
	        printf("-parallel-           p = %d  n = %d  %f millisecond\n",p,n,secs);

		gettimeofday(&t1,NULL);
		serial_baseline(n,x_zero);     // test
		gettimeofday(&t2,NULL);
		secs = (t2.tv_sec - t1.tv_sec) * 1000.0f + (t2.tv_usec - t1.tv_usec) / 1000.0f;
                printf("-serial_baseline-    p = %d  n = %d  %f millisecond\n",p,n,secs);

		gettimeofday(&t1,NULL);
		serial_matrix(n, x_zero);	// test
		gettimeofday(&t2,NULL);
                secs = (t2.tv_sec - t1.tv_sec) * 1000.0f + (t2.tv_usec - t1.tv_usec) / 1000.0f;
                printf("-serial_matrix-      p = %d  n = %d  %f millisecond\n",p,n,secs);

	}
	MPI_Finalize();
}
void parallel_serial_matrix(int x_zero, int mOff[2][2], int portion, int rank){
	int i = 0;
	int base[2] = {x_zero, 1};
	int M[2][2] = {{A,0},{B,1}};
	int Mnext[2][2];

	copyMatrix(Mnext,mOff); //initialize Mnext to M Offset

	int output;
//	printf("[process %d]   ",rank);
	for( i =0; i < portion; i++){

		output = serial_matrix_helper(base, Mnext);
		matrixMultiplication(Mnext,M);
//		printf("%d ",output);
	}
//	printf("\n");
}

void p_element_pp(int x[2][2], int p, int rank, int output[2][2]){

	int l[2][2];
	int g[2][2];
	int g_remote[2][2];
	int bit = lg2(p);
	int t = 0;
	int partner;
	MPI_Status status;
	MPI_Request request;

        copyMatrix(l,x);  // initialize local & global matrix.
        copyMatrix(g,l);

	for(t=0; t< bit; t++){
		partner = togglingTth(t,bit,rank,p);
		//send g to partner  & recv g_remote from partner and initialize g_remote
		MPI_Sendrecv(&g[0][0],4, MPI_INT, partner, 0, &g_remote[0][0], 4, MPI_INT, partner,0,  MPI_COMM_WORLD, &status  );

		if(partner < rank){
			matrixMultiplication(l,g_remote);
		}
//		matrixSum(g,g_remote);
		matrixMultiplication(g,g_remote);
	}
	copyMatrix(output,l);
}
void matrixMultiplication(int first[2][2], int second[2][2]){

	first[0][0] =  ( (first[0][0]*second[0][0]) + (first[0][1] * second[1][0]) ) % P;
	first[0][1] =  ( (first[0][0]*second[0][1]) + (first[0][1] * second[1][1]) ) % P;
	first[1][0] =  ( (first[1][0]*second[0][0]) + (first[1][1] * second[1][0]) ) % P;
	first[1][1] =  ( (first[1][0]*second[0][1]) + (first[1][1] * second[1][1]) ) % P;


}
void matrixSum(int first[2][2], int second[2][2]){

	int i=2, j =2;
	for(i = 0 ; i < 2 ; i ++){
		for(j = 0 ; j <2 ; j ++){
			first[i][j] = first[i][j] + second[i][j];
		}
	}

}
void copyMatrix(int first[2][2], int second[2][2]){

        int i=2, j =2;
        for(i = 0 ; i < 2 ; i ++){
                for(j = 0 ; j <2 ; j ++){
                        first[i][j] =  second[i][j];
                }
        }

}
void serial_baseline(int n, int x_zero){
	int * series;
	int i = 1, tem;
	series = (int*)malloc(n* sizeof(int));
	series[0] = x_zero;
//	printf("%d ",series[0]);
	for(i ; i< n; i++){
		tem = (A*series[i-1]) + B;
		series[i] = tem%P;
//		printf("%d ", series[i]);
	}
//	printf("\n");
}

int serial_matrix_helper(int first[2], int second[2][2]){
	int x_i;
	x_i =  ( (first[0] * second[0][0] ) + ( second[1][0] )  )  % P  ;

	return x_i;


}
void serial_matrix(int n, int x_zero){
	int * series;
	int i = 1, tem;
	int M[2][2] = {{A,0} , {B,1}  };
	int Mnext[2][2] = {{A,0}, {B,1}};
	int base[2] = {x_zero,1};
	series = (int*)malloc(n* sizeof(int));
	series[0] = x_zero;
//	printf("%d ",series[0]);
	for(i; i< n; i++){
		series[i] = serial_matrix_helper(base, Mnext  );
		matrixMultiplication(Mnext,M);
//		printf("%d ",series[i]);
	}
//	printf("\n");

}


void createMatrix(int power,int matrix[2][2]){ // this function is for generating M and M0.
	if(power == 0){ // if power ==0, create M^0
		matrix[0][0] = 1;
		matrix[0][1] = 0;
		matrix[1][0] = 0;
		matrix[1][1] = 1;
	}
	else{ // create M
		matrix[0][0] = A;
		matrix[0][1] = 0;
		matrix[1][0] = B;
		matrix[1][1] = 1;
	}
}
void printMatrix(int matrix[2][2], int rank) {
	int i,j;
//	printf("-----------------[%d]Matrix[%d]------------\n",rank, rank);
	for( i = 0; i < 2; i++){
		printf(" %d iter.  rank[%d] ",i,rank);
		for( j = 0 ; j < 2 ; j++){
			printf(" %d ", matrix[i][j]);
		}
		printf("\n");
	}
}


//////// blow code is from P2//////////

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



int togglingTth( int t, int bit, int rank,int p)
{
	int *binaryNum = (int *)malloc(sizeof(int)*bit);
	int init=0;
	for(init = 0; init< bit; init++){
		binaryNum[init] = 0;
	}

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

