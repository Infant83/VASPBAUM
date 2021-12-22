export DYLD_LIBRARY_PATH=$LIBRARY_PATH

DIR="/Users/Infant/iCloud/VASP/TMDC_TBFIT/WSe2-data/4.WSe2-G_hetero/0.0.00/3.OPT_SLIDING/2.case_B/RELAX_FULL/SOC/ENCUT400/BAND_UNFOLD/unfold_graphene/1.path"
#DIR="/Users/Infant/iCloud/code/GitHub/VASPBAUM/example/1.WSe2-twisted_bilayer/"
WAVECAR="$DIR/WAVECAR"
NPROC=8
make
#mpirun -np $NPROC vaspbaum.mpi -path $DIR -nosoc -unfold -sigma 0.02  -nediv 2000 -norm T -ef  -0.0
#./vaspberry.mpi -path $DIR -nosoc -unfold -sigma 0.015 -ien 0.0 -fen 4.0 -nediv 2000 -norm T
#./vaspberry.mpi -path $DIR -nosoc -unfold -sigma 0.015 -ien 0.0 -fen 4.0 -nediv 2000 -norm T
#mpirun -np $NPROC ./vaspberry.mpi -wf $WAVECAR -nosoc  -cd 2 -sigma 0.015 -ien 0.0 -fen 4.0 -nediv 2000 -out $out_folder
#./vaspberry.mpi -wf $WAVECAR -nosoc  -cd 2 -sigma 0.015 -ien 0.0 -fen 4.0 -nediv 2000 -out $out_folder
#python cd_plot_line.py

DIR="/Users/Infant/iCloud/VASP/TMDC_TBFIT/WSe2-data/5.WSe2-twist/1.m1n2_21.79/1.WSe2_only_layer/BAND_UNFOLD/nseg_18/vaspberry"
# 1. prepare unfold
#mpirun -np $NPROC vaspbaum.mpi -path $DIR  -nosoc  -set_unfold 

# 2. unfold
#mpirun -np $NPROC vaspbaum.mpi -path $DIR  -nosoc  -unfold -sigma 0.02  -nediv 2000 -norm T -ef  -0.0

# 3. CD
#mpirun -np $NPROC vaspbaum.mpi -path $DIR  -nosoc -sigma 0.02  -nediv 2000 -norm T -ef  -0.0 -cd 2 -ien 0.0 -fen 4.0

# 4. unfold & CD
#mpirun -np $NPROC vaspbaum.mpi -path $DIR  -nosoc -unfold -sigma 0.02  -nediv 2000 -norm T -ef  -0.0 -cd 2 -ien 0.0 -fen 4.0

DIR="/Users/Infant/iCloud/VASP/TMDC_TBFIT/WSe2-data/4.WSe2-G_hetero/0.0.00/3.OPT_SLIDING/2.case_B/RELAX_FULL/SOC/ENCUT400/BAND_UNFOLD/unfold_graphene/1.path"
# 2. unfold
 mpirun -np $NPROC vaspbaum.mpi -path $DIR  -soc -unfold -sigma 0.02  -nediv 2000 -norm T -ef  -0.0
# 4. unfold & CD
 mpirun -np $NPROC vaspbaum.mpi -path $DIR  -soc -unfold -sigma 0.02  -nediv 2000 -norm T -ef  -0.0 -cd 2 -ien 0.0 -fen 4.0

