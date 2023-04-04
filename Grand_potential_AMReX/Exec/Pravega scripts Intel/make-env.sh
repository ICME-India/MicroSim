#!/bin/sh

module load spack/0.17
. /home-ext/apps/spack/share/spack/setup-env.sh

module load compiler/intel/oneapi/compiler/latest

source /home/ext/pub/compiler/intel/oneapi/setvars.sh

spack load gcc /yqde46k
module load gnu11
spack load gsl /h5q52xd

export LD_LIBRARY_PATH=/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/lib:$LD_LIBRARY_PATH
