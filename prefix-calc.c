#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

int getMid(int s, int e) {  return s + (e -s)/2;  }

int min(int a, int b) { return (a<=b)?a:b; }

int *create_rand_nums(int num_elements) {
  int *rand_nums = (int *)malloc(sizeof(int) * num_elements);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() % 20);
    if(rand_nums[i]==-1)
    i--;
  }
  return rand_nums;
}

int constructSTUtil(int arr[], int ss, int se, int *st, int si)
{
    if (ss == se)
    {
        st[si] = ss;
        return st[si];
    }
    int mid = getMid(ss, se);
    st[si] =  min(constructSTUtil(arr, ss, mid, st, si*2+1) ,
              constructSTUtil(arr, mid+1, se, st, si*2+2));
    return st[si];
}
int constructSTUtil2(int arr[], int ss, int se, int *st, int si)
{
    if (ss == se)
    {
        st[si] = arr[ss];
        return arr[ss];
    }
    int mid = getMid(ss, se);
    st[si] =  -1;
    constructSTUtil2(arr, ss, mid, st, si*2+1);
    constructSTUtil2(arr, mid+1, se, st, si*2+2);
    return st[si];
}

int *constructST(int arr[], int n,int al)
{
    int x = (int)(ceil(log2(n)));
    int max_size = 2*(int)pow(2, x) - 1;
    int *st = malloc(sizeof(int)*max_size);
    if(al==0)
    constructSTUtil(arr, 0, n-1, st, 0);
    else
    constructSTUtil2(arr, 0, n-1, st, 0);
    return st;
}

void bottomUpUpdate(int tree[],int st[],int ss,int se,int si,int world_rank){
  int pre=(si&1)?(si/2):(si/2)-1;
  if(ss == se){
    if(world_rank!=0&&st[si]==world_rank&&st[pre]!=st[si]){
      MPI_Send(&tree[si], 1, MPI_INT, st[pre], 0, MPI_COMM_WORLD);
    }
    return;
  }
  int mid=getMid(ss,se);
  bottomUpUpdate(tree,st,ss,mid,si*2+1,world_rank);
  bottomUpUpdate(tree,st,mid+1,se,si*2+2,world_rank);
  if(world_rank!=st[si])
  return;
  int recv;
  MPI_Recv(&recv,1,MPI_INT,st[si*2+2],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  tree[si*2+2]=recv;
  tree[si] = tree[si*2+1] + recv;
  if(world_rank!=0&&st[pre]!=st[si]){
    MPI_Send(&tree[si], 1, MPI_INT, st[pre], 0, MPI_COMM_WORLD);
    return;
  }
}

void topDownUpdate(int tree[],int st[],int ss,int se,int si,int world_rank,int *ans){
  int left=si*2+1;
  int right=si*2+2;
  int pre=(si&1)?(si/2):(si/2)-1;
  if(world_rank==st[si]){
    if(ss==se){
      if(ss!=0&&st[pre]!=st[si]){
        MPI_Recv(&tree[si],1,MPI_INT,st[pre],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
      if(ss==0){
        ans[0]=tree[si];
      }else{
        MPI_Send(&tree[si], 1, MPI_INT, 0, ss, MPI_COMM_WORLD);
      }
      return;
    }
    if(si!=0&&st[pre]!=st[si]){
      MPI_Recv(&tree[si],1,MPI_INT,st[pre],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    tree[left]=tree[si]-tree[right];
    tree[right]=tree[si];
    MPI_Send(&tree[right], 1, MPI_INT, st[right], 0, MPI_COMM_WORLD);
  }else{
    if(ss==se)
    return;
  }
  int mid=getMid(ss,se);
  topDownUpdate(tree,st,ss,mid,si*2+1,world_rank,ans);
  topDownUpdate(tree,st,mid+1,se,si*2+2,world_rank,ans);
}

int main(int argc, char** argv){
  if (argc != 2) {
    fprintf(stderr, "Usage: n size of arrays\n");
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
  int *st = NULL;
  int *arr=NULL;
  int *tree=NULL;
  int *ans=NULL;
  int size=0;
  if(world_rank == 0){
    int x = (int)(ceil(log2(n)));
    int max_size = 2*(int)pow(2, x) - 1;
    arr=create_rand_nums(n);
    st=constructST(arr, n, 0);
    tree=constructST(arr, n, 1);
    ans=malloc(sizeof(int)*n);
    size=max_size;
  }else{
    int x = (int)(ceil(log2(n)));
    int max_size = 2*(int)pow(2, x) - 1;
    size=max_size;
    st = malloc(sizeof(int)*max_size);
    tree = malloc(sizeof(int)*max_size);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(tree, size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(st, size, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  bottomUpUpdate(tree,st,0,n-1,0,world_rank);
  topDownUpdate(tree,st,0,n-1,0,world_rank,ans);
  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank==0){
    printf("Array: \n");
    for(int i=0;i<n;i++){
      printf("%d ",arr[i]);
    }
    printf("\nSum: \n");
    for(int i=1;i<n;i++){
      MPI_Recv(&ans[i],1,MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    for(int i=0;i<n;i++){
      printf("%d ",ans[i]);
    }
    printf("\n");
  }

  MPI_Finalize();
}
