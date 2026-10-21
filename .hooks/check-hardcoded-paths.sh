#!/bin/sh

# This script is to be used as a Git pre-commit hook (see
# .pre-commit-config.yaml). It raises an error if it can find hard-coded paths
# to directories used in NeuroSpin.

! grep -E -- '\W/(neurospin|drf|volatile)' "$@"
