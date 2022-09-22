#!/bin/bash

np="$1"

nodes=`python -c "from math import ceil; print(ceil($np/8))"`
mpiprocs=`python -c "from math import ceil; print(ceil($np/$nodes))"`

echo "`< /dev/stdin`" | sed -e "s/%SELECT%/$nodes/" -e "s/%MPIPROCS%/$mpiprocs/"


