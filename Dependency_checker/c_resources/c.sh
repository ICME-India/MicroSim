#!/bin/sh


gslchkr ()
{
cd c_resources
echo "******************************************************************"
gcc ./gsl_checker.c -o ./gsl_chkr.out
gsl_status=$?
chmod 777 *
./gsl_chkr.out
#echo $gsl_status
if [ "$gsl_status" -eq 0 ]; then
echo "GSL(GNU Scientific library) is correctly installed and configured: \e[1;32m STATUS : OK \e[0m"
else
echo "GSL (GNU Scientific Library) is not installed or not properly configured. \e[1;31m STATUS : NOT FOUND \e[0m " 
fi
cd ..

echo ""
}



fftw3chkr ()
{
cd c_resources
echo "******************************************************************"
gcc ./fftw3checker.c -o ./fftw3_chkr.out
fftw3_status=$?
chmod 777 *
./fftw3_chkr.out
#echo $fftw3_status
if [ "$fftw3_status" -eq 0 ]; then
echo "fftw3 library is correctly installed and configured: \e[1;32m STATUS : OK \e[0m"
else
echo "fftw3 library is not installed or not properly configured. \e[1;31m STATUS : NOT FOUND \e[0m " 
fi
cd ..
echo ""
}

CUDA_Aware_MPI ()
{
cd c_resources
echo "******************************************************************"
gcc ./cudaAwareMPI_check.c -o ./cudaAwareMPI_chkr.out
chmod 777 *
./cudaAwareMPI_chkr
cudaMPI_status=$?
#echo $fftw3_status
if [ "$cudaMPI_status" -eq 0 ]; then
echo "CUDA aware OpenMpi is propperly installed and configured: \e[1;32m STATUS : OK \e[0m"
else
echo "CUDA aware OpenMpi is not installed or not properly configured. \e[1;31m STATUS : NOT FOUND \e[0m " 
fi
cd ..
echo ""
}
