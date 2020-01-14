/*

How To Run?

gcc-9 -fopenmp main.c -o main // compile.  in gcc-9, 9 means the version of gcc.
./main k d                    // running code.  d is double format.




*/







#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define QUEUESIZE 5 // top 5 page ranks.

typedef struct Node{        // Node structure.
    int vertexNum;          // this will store the vertex index.
    int visit;
    double pr;
    struct Node *next;      // this will point to Next now that is linked to itself.

}Node;  

typedef struct Queue{
    int count;
    Node *front;
    Node *rear;

}Queue;

void initializeQueue(Queue *q){
    q->count = 0;
    q->front = NULL;
    q->rear = NULL;
}

int isEmpty(Queue *q){
    return (q->rear ==NULL);
}
void dequeue(Queue *q){
    Node *tmp;
    tmp = q->front;
    q->front = q->front->next;
    q->count--;
    free(tmp);
}

void enqueue(Queue *q, Node *value){
    // if(q->count <QUEUESIZE){

    Node *tmp;
    tmp = malloc(sizeof(Node));
    tmp->vertexNum = value->vertexNum;
    tmp->visit = value->visit;
    tmp->pr = value->pr;
    tmp->next = NULL;

    
    if(!isEmpty(q)){

        if(q->count < 5){
            if(q->rear->pr< tmp->pr){
                q->rear->next = tmp;
                q->rear = tmp;
                q->count++;
            }
        }   
        if(q->rear->pr < tmp->pr  && q->count ==5){  // if new node has bigger rank, push.
            dequeue(q); // dequeue the lowest pr vertex because we already have 5 entries.
            q->rear->next = tmp;
            q->rear = tmp;
            q->count++;
        }
        
    }else{
        q->front = tmp;
        q->rear = tmp;
        q->count++;
    }

}

void display(Node * head,int index){

    if(head == NULL){
        printf("\n");
    }else{
        printf("[%d]Vertex:%d\nPageRank:%lf\n\n",index, head->vertexNum, head->pr);
        index--;
        display(head->next,index);
    }


}

typedef struct List{  // List structure will be used to store Head Node of a list.
    Node *head;
}List;



int getEdgeNum(List *arr, int index){    // this function returns the number of edges in vertices.
    int i = 0;
    Node *ptr;
    ptr = arr[index].head;
    if(ptr == NULL){
        return -1;
    }
    while(ptr){
        i++;
        ptr = ptr->next;
    }
    return i-1;
}

int getRandomVertex(List*arr, int size, int seed){
    Node * ptr;
    //do random. 
    while(1){
        int randomNum = rand_r(&seed) % size; // get randome number which is 0~ last index
        ptr = arr[randomNum].head;  // put in ptr
        if(ptr!= NULL){ // if it is not NULL, then return
            return randomNum; 
        }
    }
}

int getNeightborVetex(List* arr, int vertexNum, int seed ,int size){
    
    // printf(" in with seed = %d\n",seed);
    Node * tmp = arr[vertexNum].head;
    int count = 2;
    int i = 1;
    
    int marker = -1; // to make sure it doesn't get an empty neighbor as a result.
    int totalEdge = getEdgeNum(arr, vertexNum);

    int jumpTo = (seed+count) % totalEdge;
    jumpTo++; //prevent 0;
    int validResult = 0;

    for(i = 0 ; i < totalEdge; i++){ // make sure not all neighbors are empty vertex.
        if(getEdgeNum(arr,tmp->vertexNum) != -1){
            marker = 1;
            validResult = tmp->vertexNum; // save the last valid neighbor that is not empty.
            if ( jumpTo == i){
                return tmp->vertexNum;          // if selected neighbor is not empty,
            }
            
        }
        tmp = tmp->next;
    }
    if(marker == -1){ // if all the neighbors are empty, then just pick randome node out of all nodes
        int sendResult = getRandomVertex(arr, size ,seed);
        return sendResult;
    }
    return validResult;

}



void insertVertex(int a, int b, List *arr){
    Node *dest, *tmp, *src;
    if(arr[a].head == NULL){
        src = (Node*)malloc(sizeof(Node));
        src->vertexNum = a;
        src->visit = 0 ;
        src->pr = 0.0;
        src->next = NULL;
        arr[a].head = src;
    }
    

    dest = (Node*)malloc(sizeof(Node));
    dest->vertexNum = b;
    dest->visit = 0 ;
    dest->pr = 0.0;
    dest->next = NULL;

    tmp = arr[a].head;
    while(tmp->next != NULL){
        tmp = tmp->next;
    }
    tmp->next = dest;
    

    
}


void printMatrix(List *arr, int size){

    int i = 0 ;
    // int count = 0 ;
    for( i = 0 ; i <= size ; i++){
        Node * ptr = arr[i].head;
        // if(arr[i].head->vertexNum == 916427){
        //     printf("i = %d\n",i);
        //     return;
        // }
        while(ptr != NULL){
            // count++;
            printf(" %d[%d,%.2lf] ",ptr->vertexNum, ptr->visit,ptr->pr*100000);
            ptr = ptr->next;
        }
        printf("\n");
        // printf(" total = %d, edges = %d \n",count,count-1 );
    }
}

void printLine(List*arr, int index){
    int i = 0 ;
    int count = 0;
    Node *ptr = arr[index].head;
    while(ptr !=NULL){
        printf(" %d[%d] \n",ptr->vertexNum, ptr->visit);
        count++;
        ptr = ptr->next;
    }
    printf(" total = %d, edges = %d \n",count,count-1 );
    return;
}

void updateRank(List *arr, int indexForMatrix, int k, int last_index){
    int i = 0 ;

    double temVisit;

    Queue *q;       //for top 5 output.
    q = malloc(sizeof(Queue));
    initializeQueue(q);


    for( i = 0 ; i <= last_index ; i++){
        Node * ptr = arr[i].head;
        temVisit = 0.0;

        if(ptr != NULL){
            temVisit += ptr->visit;
            ptr->pr = temVisit/indexForMatrix; 
            ptr->pr = ptr->pr / k;
            enqueue(q,ptr);  // putting in queue to output top 5.
        }

    }
    int outputRank = 5;
    printf("queue size = %d\n",q->count);
    display(q->front,outputRank);
    


}



void pageEstimate(char * filename, int k, double d){
    char line[256];
    int first_index = -1;
    //double matrix[v][e];
    int i=0,j =0;
    int v=-1,e=-1;
    printf("Your file : %s\n",filename);
    FILE *fp;

    if ( (fp=fopen(filename,"r")) == NULL){
        printf("Failed to open file\n");
        return ;
    }
    int tem  = 0, lineNum =0 ;
    char * token;

    int fileFormat = 0;
    
    int incr = 0;
    while( fgets(line, 256, fp)!= NULL){ // get the number of vertices and edges.  
        if(line[0] == '#' && incr == 0){
            fileFormat = 1;
        }
        if( line[0] == '#'){   // this part detects the lines starting with # in the input file.
            token =strtok(line," ");
            while(token != NULL){   // keep tokenizing until it sees the Vertex # and Edge #
                if( strcmp(token,"Nodes:") == 0){
                    token = strtok(NULL," ");
                    v = atoll(token);
                }
                if( strcmp(token,"Edges:") == 0){
                    token = strtok(NULL," ");
                    e = atoll(token);
                }else{
                    token = strtok(NULL," ");
                }
            }
            tem++;
        }else{
            lineNum++;

            token = strtok(line,"\t");
            if( first_index == -1){
                // printf(" %s \n",line);
                first_index = atoi(token);
            }
        }
        incr ++;
    
        
    }
    
    int last_index = atoi(token);

    // printf("first index = %d, last index = %d, total lines = %d\n",first_index,last_index, lineNum);
    // return;
    v= v <0? lineNum: v; // v could be -1 depending upon file format then we use lineNum.

    List *adjacencyList = (List*)malloc(last_index * sizeof(List));  // allocate memory for the list of linked list heads.
    for( i = 0; i < last_index ; i ++){
        // printf("lst = %d i= %d v = %d \n",adjacencyList[i],i,v);
        adjacencyList[i].head = NULL;
    }

    
    // printf("size of adjacencyList = %d\n",v*sizeof(List));
    rewind(fp);
    tem = 0;
    int firstEle = -1;
    int secondEle = -1;
    int indexForMatrix =0;

    indexForMatrix = last_index- first_index;


    if( fileFormat == 0){
        // printf(" format = 0 \n");
        while(fgets(line,256,fp)!=NULL){
        
            if(line[0]!= '#'){

                if( (token = strtok(line," ")) != NULL){   //first element of pair
                    firstEle = atoi(token);
                }
                if( (token = strtok(NULL, " ")) != NULL){     //seoncd element of pair.
                    secondEle = atoi(token);
                }

                insertVertex( firstEle, secondEle, adjacencyList);
                // insertVertex(secondEle, secondEle, firstEle, adjacencyList);
            }
        }
    }

    if( fileFormat == 1){
        // printf(" format = 1 \n");
        while(fgets(line,256,fp)!=NULL){
        
            if(line[0]!= '#'){

                if( (token = strtok(line,"\t")) != NULL){   //first element of pair
                    firstEle = atoi(token);
                }
                if( (token = strtok(NULL, "\t")) != NULL){     //seoncd element of pair.
                    secondEle = atoi(token);
                }
                // if(firstEle == tem ){           // detects edges with different vertax
                //     // printf("add vertex to the same head\n");
                //     // insertVertex(firstEle, firstEle, secondEle, adjacencyList);
                // }
                // if(firstEle != tem){
                //     // printf("add vertex to the different head\n");
                //     // indexForMatrix++;
                //     // insertVertex(firstEle, firstEle, secondEle, adjacencyList);
                //     tem = firstEle;
                // }
                // printf(" insert2 before\n");
                insertVertex( firstEle, secondEle, adjacencyList);
                // insertVertex(secondEle, secondEle, firstEle, adjacencyList);
            }
        }
    }

    fclose(fp);
  
    // printMatrix(adjacencyList,indexForMatrix);
    

    
    double runTime = omp_get_wtime();
    int totalWalk= 0;

    #pragma omp parallel for schedule(static) shared(adjacencyList,indexForMatrix,totalWalk)
    for(i =0 ; i <= indexForMatrix ; i++ ){  //  index from 0 ~ # of vertexes
        // printf(" my rank = %d \n",omp_get_thread_num());    

        // printf("check ");
        if(adjacencyList[i].head != NULL){ // this if statement filters empty heads.
            Node * temNode ;
            temNode = adjacencyList[i].head;        //update the very first Node.
            #pragma omp atomic
                temNode->visit +=1; // update the visit count.
            totalWalk++;
            
            // #pragma omp parallel for schedule(static) shared(k,i) private(adjacencyList)
            for(j = 1 ; j <= k ; j++){
                int z = 0;


                int myseed = omp_get_thread_num()*(i+1)*(j+1)*last_index;
                int probability = d*100; // this is my damping ratio.
                int randomNum = rand_r(&myseed) % 101; // get toss / tail probability 


                // 1-d  == tail
                // d == head.
                if(probability >= randomNum  ){                         //    if d, jump to a random Node find a Node that is not NUll (while()) then jump ->                  
                    int randomIndex= getRandomVertex(adjacencyList,last_index,myseed);
                    temNode = adjacencyList[randomIndex].head;
                }   
                else{                  //    if 1-d, jump to a connected Node When you jump to another Node, the probability of each Node = 1/(d(i) +1 ) then jump -> 
                    int jumpToNeighbor = getNeightborVetex(adjacencyList,temNode->vertexNum, myseed,last_index);
                    temNode = adjacencyList[jumpToNeighbor].head;
                }
                // omp_set_lock(&my_lock);
                #pragma omp atomic
                    temNode->visit += 1;//update visit
                // omp_unset_lock(&my_lock);
                totalWalk++;
            }

        }

    
    }


    runTime = omp_get_wtime() - runTime;
    printf("Total time For Random Walk = %f seconds \n", runTime);

    runTime = omp_get_wtime();
    updateRank(adjacencyList,indexForMatrix, k ,last_index );
    // printMatrix(adjacencyList,last_index); 
    runTime = omp_get_wtime() - runTime;
    printf("Total time for Calculating Rank = %f seconds \n", runTime);

    // printMatrix(adjacencyList,last_index); 

}


int main(int argc, char *argv[])
{
        int V= 10, E = 10;

        
        int k=0;
        float d=0;// random walk, damping ratio
        int n = 10; // number of Node. 10 is just a random number to use temperarilly.
        if(argc<2) {
                printf("Input: {K: length of random walk} {D: damping ratio}\n");
                exit(1);
        }

        //loops = atoll(argv[1]);
        k = atoll(argv[1]);
        printf("k = %d\n",k);
        if(argc==3) {

                d = atof(argv[2]);
                if ( d> 0 && d < 1 ){
                    printf("d = %.2lf\n",d);
                }else{
                    printf( "check d value [0 < d < 1] ") ;
                    return 0;
                }
        }
        int p = 8;

        char filename[32];
        printf("Enter file name: (ex. input.txt)\n:");
        scanf("%32s",filename); //read no more than 32 characters. // prevent buffer overflow

        omp_set_num_threads(p);
        #pragma omp parallel
        {
            assert(p==omp_get_num_threads());
            //printf("Debug: number of threads set = %d\n",omp_get_num_threads());

            int rank = omp_get_thread_num();
            printf("Info: Rank=%d: my world has %d threads\n",rank,p);
        }  // end of my omp parallel region
        pageEstimate(filename,k,d);    // expected d = 90%


        
        return 0;
}