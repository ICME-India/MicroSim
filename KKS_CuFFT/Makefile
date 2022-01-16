SHELL = /bin/sh
NVCC = nvcc   
COMPOPS = -Xptxas -O3
#COMPOPS = -Xptxas -g  
ARCH = -arch=sm_70 
LIBS =  -lcufft -lcuda -lcudart
#INCS =  -I/usr/local/cuda-11.2/samples/common/inc
HEADERS = functions/global_vars.h

SRCDIR = .
BINDIR = DATA
SRCFILES = $(SRCDIR)/KKS_cuFFT.cu
all: build

build: build_dir $(HEADERS)
	$(NVCC) $(COMPOPS) $(ARCH) $(SRCFILES) -o microsim_kks_cufft $(LIBS)

build_dir:
	mkdir -p $(BINDIR)

clear:
	-\rm $(BINDIR)/*.*

clean:
	-\rm *.out
	-\rm -r $(BINDIR)
#	@\rm -f job.*.*
## End of the makefile
