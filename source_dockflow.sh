#!/bin/bash
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
export DOCKFLOW_ROOT=$SCRIPT_DIR
alias dockflow="$DOCKFLOW_ROOT/full_pipeline/dockflow"
