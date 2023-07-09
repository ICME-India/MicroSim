module load spack/0.17
. /home-ext/apps/spack/share/spack/setup-env.sh

#spack load cuda /v75gxwp
#spack load gsl  /h5q52xd

#spack load cuda /pbeecsj
#spack load gsl  /a7p4vaz
spack load intel-mpi@2019.10.317 /6icwzn3
spack load intel-oneapi-compilers@2021.4.0
spack load intel-oneapi-mkl /gsap3zd
spack load hdf5  /ssgjscn

spack load cuda /rbt77ll
spack load gsl  /h5q52xd

make -f Makefile.pravega 
 

