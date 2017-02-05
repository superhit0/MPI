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

float *calc_row(float *a,float *b[],int n,int o){
  float *ans=malloc(sizeof(float)*o);
  for(int i=0;i<o;i++){
    float sum=0;
    for(int j=0;j<n;j++){
      sum+=(a[j]**(b[i]+j));
    }
    ans[i]=sum;
  }
  return ans;
}

int main(int argc, char** argv){
  if (argc != 4) {
    fprintf(stderr, "Usage: m n o of arrays\n");
    exit(1);
  }

  int m,n,o;
  m=atoi(argv[1]);
  n=atoi(argv[2]);
  o=atoi(argv[3]);
  srand(time(NULL));

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int lim,*nop=NULL,count;

  float *a[m],*b[o],*c[m];
  if(world_rank==0){
    for(int i=0;i<m;i++){
      a[i]=create_rand_nums(n);
    }
    for(int i=0;i<o;i++){
      b[i]=create_rand_nums(n);
    }
    for(int i=0;i<m;i++){
      c[i]=malloc(sizeof(float)*o);
    }
    if(world_size>=m){
      nop=malloc(sizeof(int)*1);
      count=1;
      lim=m-1;
    }else{
      lim=world_size-1;
      count=world_size;
      nop=malloc(sizeof(int)*world_size);
      for(int i=0;i<world_size;i++){
        nop[i]=0;
      }
      for(int i=0;i<m;i++){
        nop[i%world_size]++;
      }
    }
  }else{
    for(int i=0;i<m;i++){
      a[i]=malloc(sizeof(float)*n);
    }
    for(int i=0;i<o;i++){
      b[i]=malloc(sizeof(float)*n);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&lim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank!=0)
  nop=malloc(sizeof(int)*count);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(nop, count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<m;i++){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(a[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    for(int i=0;i<o;i++){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(b[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(world_rank<=lim){
    int offset=0;
    for(int i=0;i<world_rank;i++){
      offset+=count==1?1:nop[i];
    }
    //int end=count==1?offset+1:offset+nop[world_rank];
    int no_rows=count==1?1:nop[world_rank];
    float *sub_ans[no_rows];
    for(int i=0,j=offset;i<no_rows;i++,j++){
      sub_ans[i]=calc_row(a[j],b,n,o);
      if(world_rank==0){
        c[i]=sub_ans[i];
      }else{
        MPI_Send(sub_ans[i], o, MPI_FLOAT, 0, j, MPI_COMM_WORLD);
      }
    }
  }

  if(world_rank==0){
    int row2proc[m];
    int start=0;
    if(count==1){
      start=1;
      for(int i=0;i<m;i++){
        row2proc[i]=i;
      }
    }else{
      for(int i=0,x=0;x<m&&i<world_size;i++){
        if(i==1)start=x;
        for(int j=0;x<m&&j<nop[i];j++,x++){
          row2proc[x]=i;
        }
      }
    }
    for(int i=start;i<m;i++){
      MPI_Recv(c[i],o,MPI_FLOAT,row2proc[i],i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    printf("A=\n");
    for(int i=0;i<m;i++){
      for(int j=0;j<n;j++){
        printf("%f ",a[i][j]);
      }
      printf("\n");
    }
    printf("B=\n");
    for(int i=0;i<n;i++){
      for(int j=0;j<o;j++){
        printf("%f ",b[i][j]);
      }
      printf("\n");
    }
    printf("C=\n");
    for(int i=0;i<m;i++){
      for(int j=0;j<o;j++){
        printf("%f ",c[i][j]);
      }
      printf("\n");
    }
  }
  printf("\n");
  MPI_Finalize();
}
