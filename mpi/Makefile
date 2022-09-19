
# Example executions:
# 1) $ make run-release
#	Will run the release executable with default problem size 840.
# 2) $ make run-release size=1680 n=8
#	Will run the release executable with supported problem size 1680 and # of processes 8.

# Receives argument "size". Default value: 840.
size ?= 840
n ?= 4

THIS_SOURCES := $(shell find ./src/ -name "*.c")

default:
	mkdir -p ./bin

compile-release: default
	mpicc -O3 ${THIS_SOURCES} -lm -o ./bin/mpi_release.x

# MPI profile
compile-debug: default
	mpicc -O3 -g ${THIS_SOURCES} -L/opt/mpiP-3.5/build/lib -lmpiP -lbfd -lunwind \
	-lm -o./bin/mpi_debug.x

compile-all: compile-debug compile-release

run-debug:
	cat ../inputs/input_template.txt | sed -e 's/%SIZE%/${size}/g' | mpirun -np ${n} ./bin/mpi_debug.x
run-release:
	cat ../inputs/input_template.txt | sed -e 's/%SIZE%/${size}/g' | mpirun -np ${n} ./bin/mpi_release.x

clean:
	find ./bin/ -name "*.x" -exec rm {} \;