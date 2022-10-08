
# Problem size.
size ?= 840
# Number of processes when using the "qsub-*" targets.
ngpus ?= 1

####################################################################################################

CUDA_INSTALL_PATH = /usr/local/cuda-11.6
CC = g++
OPTFLAG = -O2 -fomit-frame-pointer -ftree-vectorize -ftree-vectorizer-verbose=0 -funroll-loops
NVCC = ${CUDA_INSTALL_PATH}/bin/nvcc
INCDIR = -I./src/cuda/common/inc/
FLAGS = ${OPTFLAG} -I${CUDA_INSTALL_PATH}/include -Wall -g ${INCDIR}
NVFLAGS = -O2 -I${CUDA_INSTALL_PATH}/include --compiler-options -fno-strict-aliasing --ptxas-options=-v -g ${INCDIR}
BITS = $(shell getconf LONG_BIT)
ifeq (${BITS},64)
	LIBSUFFIX := 64
endif
IFLAGS = -I../common/include/ -I./include/
LFLAGS = -L${CUDA_INSTALL_PATH}/lib${LIBSUFFIX} -lm -lstdc++ -lcudart

####################################################################################################

# <shortcut targets
default: compile
x: qsub
# /shortcut targets>

.prepare-compilation:
	@mkdir -p ./argo/bin
compile: .prepare-compilation
	${CC} ${FLAGS} ${IFLAGS} -c ./src/main.cpp \
		-o ./argo/bin/main.o
	${CC} ${FLAGS} ${IFLAGS} -c ../common/src/check_solution.c \
		-o ./argo/bin/check_solution.o
	${NVCC} ${NVFLAGS} ${IFLAGS} -DUNIX -c ./src/jacobi_gpu.cu \
		-o ./argo/bin/jacobi_gpu.o
	${CC} ${LFLAGS} ./argo/bin/main.o ./argo/bin/check_solution.o ./argo/bin/jacobi_gpu.o \
		-o ./argo/bin/jacobi_gpu.x

.prepare-qsub:
	@mkdir -p ./argo/inputs ./argo/outputs &&\
	cp ../inputs/input_template.txt ./argo/inputs/input_${size}.txt &&\
	sed -i 's/%SIZE%/${size}/g' ./argo/inputs/input_${size}.txt
qsub: .prepare-qsub
	@cat ./PBS_template.sh |\
		sed -e 's/%NGPUS%/${ngpus}/' -e 's/%SIZE%/${size}/' |\
		tee ./argo/inputs/temp.txt |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id &&\
		../common/scripts/rename_pbs_config.sh ./argo/inputs/temp.txt

clean-bin:
	rm -f ./argo/bin/*
clean-inputs:
	rm -f ./argo/inputs/*
clean-outputs:
	rm -f .latest_qsub_job_id ./argo/outputs/*
clean-all:
	rm -rf .latest_qsub_job_id ./argo/

