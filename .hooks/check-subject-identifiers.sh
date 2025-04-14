#!/bin/sh

# This script is to be used as a Git pre-commit hook (see
# .pre-commit-config.yaml). It raises an error if it can find internal subject
# identifiers in the committed scripts.

! grep -E -- '2[1-5]F[0-9]{3,}|LA[0-9]{5,}' "$@"
