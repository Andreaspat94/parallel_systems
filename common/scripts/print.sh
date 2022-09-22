#!/bin/bash

GREEN='\033[32m'
RED='\033[31m'
NC='\033[0m'

function print_red {
    echo -en "${RED}"
    echo "$@"
    echo -en "${NC}"
}

function print_green {
    echo -en "${GREEN}"
    echo "$@"
    echo -en "${NC}"
}

function print_error {
    >&2 print_red $@
}


