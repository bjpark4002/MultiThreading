#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <assert.h>

void piEst(int n, int p, struct drand48_data * drand)
{
	int i = 0;
	double sum=0;
	#pragma omp parallel for private(i) shared(n,drand) reduction(+:sum) schedule(static)
	for (i = 0 ; i < n ; i++){
		double x,y;		
		drand48_r(drand,&x);
		drand48_r(drand,&y);
		if( (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5) <= 0.25){
			sum+=1.0;
		} 
	}
	printf("sum = %f success rate = %f  pi =%f \n",sum, sum/(double)n, 4.0*sum/(double)n);

}




int main(int argc, char *argv[])
{
	int n, p;
	double x,y;
	struct drand48_data drandVar;
	

	int zz,z = 0; 

	// loop {number of iterations} [number of threads]

	if(argc<2) {
		printf("Usage: p4 n p\n");
		exit(1);
	}
	
	n = atoll(argv[1]);
	printf("n = %d\n",n);
	if(argc==3) {
		p = atoi(argv[2]);
		assert(p>=1);
		printf("p = %d\n",p);
	}
	omp_set_num_threads(p);
	#pragma omp parallel
	{
		assert(p==omp_get_num_threads());
		int rank = omp_get_thread_num();
		printf("Rank=%d: my world has %d threads\n",rank,p);
	}  // end of my omp parallel region

	int i,a=0,sum =0;
	double time = omp_get_wtime();
	srand48_r(n,&drandVar);
	piEst(n,p,&drandVar);
	
	time = omp_get_wtime() - time;
	printf("\n %f seconds \n ", time);

	return 0;
}