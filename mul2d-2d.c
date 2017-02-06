#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);
  }
  return rand_nums;
}

float calc_ele(float *a,float *b,int n){
  float ans=0;
  for(int i=0;i<n;i++){
    ans+=(a[i]*b[i]);
  }
  //float *p=&ans;
  return ans;
}

int main(int argc, char** argv){
  if (argc != 2) {
    fprintf(stderr, "Usage: n of arrays[nXn matrix]\n");
    exit(1);
  }

  int n;
  n=atoi(argv[1]);
  srand(time(NULL));

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int sq=(int)sqrt(world_size);
  float *a[n],*b[n],c[n][n];
  int *row,*rowoff,proc[n][n];
  row=malloc(sizeof(int)*sq);
  rowoff=malloc(sizeof(int)*sq);
  if(world_rank==0){
    for(int i=0;i<n;i++){
      a[i]=create_rand_nums(n);
    }
    for(int i=0;i<n;i++){
      b[i]=create_rand_nums(n);
    }
    for(int i=0;i<sq;i++){
      row[i]=n/sq;
    }
    for(int i=0;i<n%sq;i++){
      row[i]++;
    }
    rowoff[0]=0;
    for(int i=1;i<sq;i++){
      rowoff[i]=rowoff[i-1]+row[i-1];
    }
    for(int i=0;i<sq;i++){
      for(int j=0;j<sq;j++){
        for(int k=0;k<row[i];k++){
          for(int l=0;l<row[j];l++){
            proc[rowoff[i]+k][rowoff[j]+l]=(i*sq)+j;
          }
        }
      }
    }
  }else{
    for(int i=0;i<n;i++){
      a[i]=malloc(sizeof(float)*n);
    }
    for(int i=0;i<n;i++){
      b[i]=malloc(sizeof(float)*n);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(row, sq, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rowoff, sq, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<n;i++){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(a[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    for(int i=0;i<n;i++){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(b[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(world_rank<sq*sq){
    for(int i=rowoff[world_rank/sq];i<rowoff[world_rank/sq]+row[world_rank/sq];i++){
      for(int j=rowoff[world_rank%sq];j<rowoff[world_rank%sq]+row[world_rank%sq];j++){
        float val=calc_ele(a[i],b[j],n);
        if(world_rank==0)
        c[i][j]=val;
        else
        MPI_Send(&val, 1, MPI_FLOAT, 0, (i*n)+j, MPI_COMM_WORLD);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank==0){
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        if(proc[i][j]!=0){
          MPI_Recv(&c[i][j],1,MPI_FLOAT,proc[i][j],(i*n)+j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }
    }
    printf("A=\n");
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        printf("%f ",a[i][j]);
      }
      printf("\n");
    }
    printf("B=\n");
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        printf("%f ",b[i][j]);
      }
      printf("\n");
    }
    printf("C=\n");
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        printf("%f ",c[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Finalize();
}
