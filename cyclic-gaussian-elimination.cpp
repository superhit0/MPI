// Parallel Gaussian Elimination (cyclic-striped mapping)
// Programmed by Shawn Bratcher and Jeff Howe
// CSE 160 Term Project (Fall 96)

#include <iostream>
#include <stdlib.h>

#include "mpi.h"
using namespace std;
int size,rank; 

void serial_gaussian( double *A, double *b, double *y, int n )
{
  int i, j, k;

  cout << "In serial algorithm" << endl;

  for( k=0; k<n; k++ ) {		// k = current row
    for( j=k+1; j<n; j++ ) {		// in division step
      if( A[k*n+k] != 0)
        A[k*n+j] = A[k*n+j] / A[k*n+k];
      else
        A[k*n+j] = 0;
    }

    if( A[k*n+k] != 0 )			// calculates new value
      y[k] = b[k] / A[k*n+k];		// for equation solution
    else
      y[k] = 0.0;

    A[k*n+k] = 1.0;			// sets UTM diagonal value

    for( i=k+1; i<n; i++ ) {		// Guassian elimination occurs
      for( j=k+1; j<n; j++ )		// in all remaining rows
        A[i*n+j] -= A[i*n+k] * A[k*n+j];

      b[i] -= A[i*n+k] * y[k];
      A[i*n+k] = 0.0;
    }
  }
}

// distribute_matrix
// distributes an x by y matrix among the processors using the cyclic-mapping
// algorithim.
void distribute_matrix( double *M, double *LM, int x, int y )
{
  MPI_Status status;
  int i, j, p;

  if( rank == 0 ) {
    for( p=size-1; p>=0; p-- ) {
      for( i=p; i<y; i+=size )
        for( j=0; j<x; j++ ) {
          LM[ (i/size)*x+j ] = M[ i*x+j ];
        }
      if( p != 0 ) 
        MPI_Send( LM, (y/size)*x, MPI_DOUBLE, p, 10, MPI_COMM_WORLD );
    }
  }
  else
    MPI_Recv( LM, (y/size)*x, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status );
}

// gather_matrix
// gathers an x by y matrix distributed among the processors using the cyclic-mapping
// algorithim.
void gather_matrix( double *M, double *LM, int x, int y )
{
  MPI_Status status;
  int i, j, p;

  if( rank == 0 ) {
    for( p=0; p<size; p++ ) {
      if( p != 0 ) 
        MPI_Recv( LM, (y/size)*x, MPI_DOUBLE, p, 10, MPI_COMM_WORLD, &status );
      for( i=p; i<y; i+=size )
        for( j=0; j<x; j++ )
          M[ i*x+j ] = LM[ (i/size)*x+j ];
    }
  }
  else
    MPI_Send( LM, (y/size)*x, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
}

// print_matrix
// For debugging purposes, prints an x by y matrix.
void print_matrix( double *M, int x, int y ) 
{
  int i, j;

  for( i=0; i<y; i++ ) {
    for( j=0; j<x; j++ )  
      cout << M[ i*x+j ] << "  ";
    cout << endl;
  }
}

// parrallel_gaussian
// This is the parallel Gaussian elimination algorithim using cyclic-striped
// mapping as well as pipelining.
void parallel_gaussian( double *A, double *b, double *y, int n )
{
  MPI_Status status;
  MPI_Request request;
  MPI_Comm com = MPI_COMM_WORLD;
  double *LA = new double[ n*n/size ];		// declare local storage for
  double *Lb = new double[ n/size ];		// matricies
  double *Ly = new double[ n/size ];
  double *cur_row = new double[n+1];
  double stime, ftime, tottime=0;
  int i, j, k, np = (int) n/size;
  int pred, succ, start, flag=1;
  int ksn, in, ink;

  pred = ((rank-1)) % size;
  succ = (rank+1) % size;

  distribute_matrix( A, LA, n, n );		// Distribute matrix A and b
  distribute_matrix( b, Lb, 1, n );
  
  for( k=0; k<n; k++ ) {
    if( (k%size) == rank ) {
      ksn = (k/size)*n;
      for( j=k+1; j<n; j++ )					// Peform division step
        LA[ ksn+j ] = LA[ ksn+j ] / LA[ ksn+k ];
      Ly[ k/size ] = Lb[ k/size ] / LA[ ksn+k ];
      LA[ ksn+k ] = 1.0;

      for( j=0; j<n; j++ )
        cur_row[j] = LA[ (k/size)*n+j ];
      cur_row[n] = Ly[ k/size ];

      MPI_Send( cur_row, n+1, MPI_DOUBLE, succ, 20, com );	// Send row to successor
    }
    else {			// Receive row from predecessor
      MPI_Recv( cur_row, n+1, MPI_DOUBLE, pred, 20, com, &status );
      if( succ != k%size ) {	// Forward row to successor
        MPI_Send( cur_row, n+1, MPI_DOUBLE, succ, 20, com );
      }
    }
// This line used to test cyclic-mapping without pipelining.
//    MPI_Bcast( cur_row, n+1, MPI_DOUBLE, k%size, com );

    start = (rank <= k%size) ? (int) k/size+1 : (int) k/size;	// compute starting point

    for( i=start; i<np; i++ ) {			// Perform elimination step
      in = i*n; ink = in+k;
      for( j=k+1; j<n; j++ ) {
        LA[ in+j ] -= LA[ ink ] * cur_row[j];
      }
      Lb[i] -= LA[ ink ] * cur_row[n];
      LA[ ink ] = 0.0;
    }
  } 
  gather_matrix( A, LA, n, n );		// Gather matricies
  gather_matrix( y, Ly, 1, n );
} 

// print_equations
// takes a matrix of size n*n and an array of size n and prints
// their results onto stdout.  Because A is a UTM at the end
// stage of the program, zero values are excluded from the printout
void print_equations( double *A, double *y, int n )
{
  int i, j;

  for( i=0; i<n; i++ ) {
    for( j=0; j<n; j++ ) {
      if( A[i*n+j] != 0 ) {
        cout << A[i*n+j] << "x" << j;
        if( j<n-1 ) cout << " + ";
      }
      else
        cout << "      ";
    }
    cout << " = " << y[i] << endl;
  }
}

// print_solution
// For testing purposes, performs back-substition and prints out
// solutions for all unkowns.
void print_solution( double *A, double *y, int n )
{
  int i, j;
  double *x = new double[n];

  for( i=n-1; i>=0; i-- ) {
    x[i] = y[i];
    for( j=n-1; j>i; j-- )
      x[i] -= x[j] * A[i*n+j];
  }
  for( i=0; i<n; i+=4 ) {
    for( j=i; j<i+4 && j<n; j++ )
      cout << "x[" << j << "] = " << x[j] << "  ";
    cout << endl;
  }
}

main( int argc, char *argv[] )
{
  double *A, *b, *y;
  double stime, ftime;
  int i, j, n, r;
  
  if( argc < 2 ) {
    cout << "Usage\n";
    cout << "  Arg1 = number of equations / unkowns\n";
    return -1;
  }

  n = atoi(argv[1]);

  if( rank==0 ) {
    A = new double[n*n];	// Processor 0 contains the full matrix
    b = new double[n];
    y = new double[n];
  }

  for( i=0; i<n; i++ ) {	// Compute coefficient values using random
    b[i] = 0.0;			// numbers to prevent scalar of equations from
    for( j=0; j<n; j++ ) {	// forming.  Set up solutions to be 
      r = rand();		// x0=0,x1=1,...,xn=n
      A[i*n+j] = r;
      b[i] += j*r;
    }
  }

  MPI_Init (&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  stime = MPI_Wtime();			// Start Timing

  if (size == 1)
        serial_gaussian ( A, b, y, n );
  else
  {
    if ( ( n % size ) != 0 ) {
        cout << "Unknowns must be multiple of processors." << endl;
        return -1;
    }
    if ( rank==0 )
        cout << "In parallel algorithm: " << n/size << endl;

    parallel_gaussian( A, b, y, n );
  }
       
  ftime = MPI_Wtime();			// Finish Timing

  MPI_Barrier( MPI_COMM_WORLD );

  if ( rank==0 )
  {
     print_equations( A, y, n );
     print_solution( A, y, n );
     cout << "n = " << n << "  p = " << size;
     cout << "  x[" << n-1 << "] = " << y[n-1] << endl;
     cout << "  Total Time = " << ftime-stime << endl;
  }
}
