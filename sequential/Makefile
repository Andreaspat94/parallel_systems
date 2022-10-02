
# Example executions:
# 1) $ make run-release
#	Will run the release executable with the default problem size 840.
# 2) $ make run-release size=1680
#	Will run the release executable with problem size 1680.

# Problem size.
size ?= 840
# Number of processes when using the "qsub-*" targets.
np ?= 1
# Optimization to use. Valid values: {0,1,1x,2,2x,3,3x}.
opt ?= 0

COMMON_SOURCES := $(shell find ../common/src/ -name "*.c")
THIS_SOURCES := $(shell find ./src/ -name "*.c")

# <shortcut targets
default: compile-all
x: qsub-release
s: show-latest-qsub-info
# /shortcut targets>

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
	@./argo/bin/debug.x < ./argo/inputs/input_${size}.txt -O${opt}
run-release: .prepare-run
	@./argo/bin/release.x < ./argo/inputs/input_${size}.txt -O${opt}

.prepare-qsub: .prepare-run
	@mkdir -p ./argo/outputs
qsub-debug: .prepare-qsub
	@cat ./PBS_template.sh |\
		sed -e 's/%NP%/${np}/g' -e 's/%TARGET%/debug/' -e 's/%SIZE%/${size}/' -e 's/%OPT_LEVEL%/${opt}/' |\
		tee ./argo/inputs/temp.txt |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id &&\
		../common/scripts/rename_pbs_config.sh ./argo/inputs/temp.txt
qsub-release: .prepare-qsub
	@cat ./PBS_template.sh |\
		sed -e 's/%NP%/${np}/g' -e 's/%TARGET%/release/' -e 's/%SIZE%/${size}/' -e 's/%OPT_LEVEL%/${opt}/' |\
		tee ./argo/inputs/temp.txt |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id &&\
		../common/scripts/rename_pbs_config.sh ./argo/inputs/temp.txt

show-latest-qsub-info:
	@../common/scripts/show_qsub_job_info.sh || true

report-optimizations-performance: compile-all
	@mkdir -p ./argo/reports &&\
	./scripts/report_optimizations_performance.sh

clean-bin:
	rm -f ./argo/bin/*
clean-inputs:
	rm -f ./argo/inputs/*
clean-outputs:
	rm -f .latest_qsub_job_id ./argo/outputs/*
clean-all:
	rm -rf .latest_qsub_job_id ./argo/

