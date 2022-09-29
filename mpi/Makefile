
# Example executions:
# 1) $ make run-release
#	Will run the release executable with the default problem size 840 and # of processes 4.
# 2) $ make run-release size=1680 np=8
#	Will run the release executable with problem size 1680 and # of processes 8.

size ?= 840	# Problem size.
np ?= 4		# Number of processes when using the "qsub-*" targets.

COMMON_SOURCES := $(shell find ../common/src/ -name "*.c")
THIS_SOURCES := $(shell find ./src/ -name "*.c")

# <shortcut targets
default: compile-all
x: qsub-release
s: show-latest-qsub-info
# /shortcut targets>

.prepare-compilation:
	@mkdir -p ./argo/bin
compile-profiled: .prepare-compilation
	mpicc -O3 -g ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-L/opt/mpiP-3.5/build/lib -lmpiP -lbfd -lunwind \
		-lm -o ./argo/bin/profiled.x
compile-release: .prepare-compilation
	mpicc -O3 ${COMMON_SOURCES} ${THIS_SOURCES} \
		-I../common/include/ -I./include/ \
		-lm -o ./argo/bin/release.x
compile-all: compile-profiled compile-release

.prepare-qsub:
	@mkdir -p ./argo/inputs &&\
	cp ../inputs/input_template.txt ./argo/inputs/input_${size}.txt &&\
	sed -i 's/%SIZE%/${size}/g' ./argo/inputs/input_${size}.txt &&\
	mkdir -p ./argo/outputs
qsub-profiled: .prepare-qsub
	@cat ./PBS_template.sh |\
		./scripts/replace_spread_variables.sh ${np} |\
		sed -e 's/%TARGET%/profiled/' -e 's/%SIZE%/${size}/' -e 's/%NP%/${np}/' |\
		tee ./argo/inputs/temp.txt |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id &&\
		../common/scripts/rename_pbs_config.sh ./argo/inputs/temp.txt
qsub-release: .prepare-qsub
	@cat ./PBS_template.sh |\
		./scripts/replace_spread_variables.sh ${np} |\
		sed -e 's/%TARGET%/release/' -e 's/%SIZE%/${size}/' -e 's/%NP%/${np}/' |\
		tee ./argo/inputs/temp.txt |\
		qsub -o ./argo/outputs/ -e ./argo/outputs/ |\
		tee .latest_qsub_job_id &&\
		../common/scripts/rename_pbs_config.sh ./argo/inputs/temp.txt

show-latest-qsub-info:
	@../common/scripts/show_qsub_job_info.sh || true

clean-bin:
	rm -f ./argo/bin/*
clean-inputs:
	rm -f ./argo/inputs/*
clean-outputs:
	rm -f .latest_qsub_job_id ./argo/outputs/*
clean-all:
	rm -rf .latest_qsub_job_id ./argo/

