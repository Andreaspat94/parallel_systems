
# Example executions:
# 1) $ make run-release
#	Will run the release executable with the default problem size 840.
# 2) $ make run-release size=1680
#	Will run the release executable with problem size 1680.

size ?= 840

COMMON_SOURCES := $(shell find ../common/src/ -name "*.c")
THIS_SOURCES := $(shell find ./src/ -name "*.c")

.prepare-compilation:
	@mkdir -p ./argo/bin
compile-debug: .prepare-compilation
	mpicc -g ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-lm -o ./argo/bin/debug.x
compile-release: .prepare-compilation
	mpicc -O3 ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-lm -o ./argo/bin/release.x
compile-all: compile-debug compile-release

.prepare-run:
	@mkdir -p ./argo/inputs &&\
	cp ../inputs/input_template.txt ./argo/inputs/input_${size}.txt &&\
	sed -i 's/%SIZE%/${size}/g' ./argo/inputs/input_${size}.txt
run-debug: .prepare-run
	@./argo/bin/debug.x < ./argo/inputs/input_${size}.txt
run-release: .prepare-run
	@./argo/bin/release.x < ./argo/inputs/input_${size}.txt

.prepare-qsub: .prepare-run
	@mkdir -p ./argo/outputs
qsub-debug: .prepare-qsub
	@cat ./PBS_template.sh |\
		sed -e 's/%TARGET%/debug/' -e 's/%SIZE%/${size}/' |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id
qsub-release: .prepare-qsub
	@cat ./PBS_template.sh |\
		sed -e 's/%TARGET%/release/' -e 's/%SIZE%/${size}/' |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id

show-latest-qsub-outputs:
	@../common/scripts/show_qsub_job_outputs.sh || true

clean-bin:
	rm -f ./argo/bin/*
clean-inputs:
	rm -f ./argo/inputs/*
clean-outputs:
	rm -f .latest_qsub_job_id ./argo/outputs/*
clean-all:
	rm -rf .latest_qsub_job_id ./argo/

