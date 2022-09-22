#!/bin/bash

RESET='\033[39m'
LIGHT_RED='\033[91m'
LIGHT_GREEN='\033[92m'
LIGHT_YELLOW='\033[93m'
LIGHT_CYAN='\033[96m'
#RESET='\033[0m'

function print_red {
    echo -en "${LIGHT_RED}"
    echo "$@"
    echo -en "${RESET}"
}

function print_green {
    echo -en "${LIGHT_GREEN}"
    echo "$@"
    echo -en "${RESET}"
}

function print_yellow {
    echo -en "${LIGHT_YELLOW}"
    echo "$@"
    echo -en "${RESET}"
}

function print_cyan {
    echo -en "${LIGHT_CYAN}"
    echo "$@"
    echo -en "${RESET}"
}

function print_error {
    >&2 print_red $@
}


