# MPI
Example Programs of MPICH2

#Giude to Install MPICH2 on a system
http://mpitutorial.com/tutorials/installing-mpich2/

#Compile
mpicc -o [executable-file] [c-source-program-file].c

#Run
mpirun -n [no. of processes] --hosts [hosts list(IPs can be used as well)] ./[executable-file]
