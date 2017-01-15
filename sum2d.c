#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);
  }
  return rand_nums;
}

int calc(int rank,float *suba,float *subb,float *subc,int sub_len){
for(int i=0;i<sub_len;i++){
    subc[i]=suba[i]+subb[i];
}
return rank;
}

int main(int argc, char** argv){
if (argc != 3) {
    fprintf(stderr, "Usage: m n of arrays\n");
    exit(1);
  }
int m,n;
m=atoi(argv[1]);
n=atoi(argv[2]);
srand(time(NULL));

MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  float *a=NULL,*b=NULL,*c=NULL;

  if(world_rank==0){
    printf("Creating a from process id=%d\n",world_rank);
    a=create_rand_nums(m*n);
    printf("Creating b from process id=%d\n",world_rank);
    b=create_rand_nums(m*n);
    c=(float *)malloc(sizeof(float)*m*n);
  }

    int sub_len=0;
  if(world_rank==0){
    sub_len=(m*n)/world_size;
      if((m*n)%world_size==0){
      ;
      }else{
      sub_len++;
      }
  }

  MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&sub_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    int *sub_len_arr=NULL,*offset=NULL;

    if(world_rank==0){
        sub_len_arr=malloc(sizeof(int)*world_size);
        offset=malloc(sizeof(int)*world_size);
        sub_len_arr[0]=sub_len;
        offset[0]=0;
        for(int i=1;i<world_size;i++){
        offset[i]=i*sub_len;
        if(i==world_size-1)
        sub_len_arr[i]=m*n-offset[i];
        else
        sub_len_arr[i]=sub_len;
        }
    }
    if((m*n)%world_size!=0){
        if(world_rank==0&&(m*n)%world_size!=0){
        int x=sub_len_arr[world_size-1];
        MPI_Send(&x, 1, MPI_INT, world_size-1, 0, MPI_COMM_WORLD);
        }else if(world_rank==world_size-1){
        MPI_Recv(&sub_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    float *suba=(float *)malloc(sizeof(float)*sub_len);
    float *subb=(float *)malloc(sizeof(float)*sub_len);
    float *subc=(float *)malloc(sizeof(float)*sub_len);

    MPI_Scatterv(a, sub_len_arr, offset, MPI_FLOAT, suba,
              sub_len, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b, sub_len_arr, offset, MPI_FLOAT, subb,
              sub_len, MPI_FLOAT, 0, MPI_COMM_WORLD);
    printf("Calculation Complete in process id=%d\n",calc(world_rank,suba,subb,subc,sub_len));

    MPI_Gatherv(subc, sub_len, MPI_FLOAT, c, sub_len_arr, offset, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank==0){

    printf("Array 1:\n");
    for(int i=0;i<m*n;i++){
    printf("%f ",a[i]);
    if((i+1)%n==0)
    printf("\n");
    }

    printf("Array 2:\n");
    for(int i=0;i<m*n;i++){
    printf("%f ",b[i]);
    if((i+1)%n==0)
    printf("\n");
    }

    printf("Array 3:\n");
    for(int i=0;i<m*n;i++){
    printf("%f ",c[i]);
    if((i+1)%n==0)
    printf("\n");
    }

    free(a);
    free(b);
    free(c);

    free(sub_len_arr);
    free(offset);
    }

    free(suba);
    free(subb);
    free(subc);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
