#!/bin/bash

this_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

export PYTHONPATH=${this_directory}
export SVJ_PROCESSING_ROOT=${this_directory}

