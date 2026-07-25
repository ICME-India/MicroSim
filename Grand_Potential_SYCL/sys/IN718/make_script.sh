#bin/bash
cd ..
make clean
make DEBUG=1
cd IN718

../KKS_FD_solver Input_tdb_new.in Filling.in output
