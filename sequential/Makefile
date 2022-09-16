
# Example executions:
# 1) $ make run-release
#	Will run the release executable with default problem size 840.
# 2) $ make run-release size=1680
#	Will run the release executable with supported problem size 1680.

# Receives argument "size". Default value: 840.
size ?= 840

COMMON_SOURCES := $(shell find ../common/src/ -name "*.c")
THIS_SOURCES := $(shell find ./src/ -name "*.c")

default:
	mkdir -p ./bin

compile-debug: default
	mpicc -g ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-lm -o ./bin/sequential_debug.x
compile-release: default
	mpicc -O3 ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-lm -o ./bin/sequential_release.x
compile-all: compile-debug compile-release

run-debug:
	cat ../inputs/input_template.txt | sed -e 's/%SIZE%/${size}/g' | ./bin/sequential_debug.x
run-release:
	cat ../inputs/input_template.txt | sed -e 's/%SIZE%/${size}/g' | ./bin/sequential_release.x

clean:
	find ./bin/ -name "*.x" -exec rm {} \;
