#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){
  if (argc != 2) {
      fprintf(stderr, "Usage: n of partitions\n");
      exit(1);
    }
  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  long n;
  n=atol(argv[1]);
  double base=(double)1/(double)n;

  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Bcast(&base, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  long div=n/world_size;
  double sum=0;
  long end=world_rank==world_size-1?n-1:(world_rank+1)*div;
  for(long i=world_rank*div;i<end;i++){
    sum+=((base*4)/(1+(i*base)*(i*base)));
  }
  double allsum=0;
  MPI_Reduce(&sum,&allsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank==0){
    printf("PI= %.17g \n",allsum);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}
