export DYLD_LIBRARY_PATH=$LIBRARY_PATH

DIR="./"  # PATH to the folder where input files exist. Output will be generated here as well.
WAVECAR="./WAVECAR"  # path to the wavecar is located
NPROC=10
EXEC="../../vaspbaum.mpi"  # path to the "vaspbaum" program
mpirun -np $NPROC $EXEC -path $DIR -set_unfold 
