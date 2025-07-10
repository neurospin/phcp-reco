#!/bin/sh

# This script is to be used as a Git pre-commit hook (see
# .pre-commit-config.yaml). It raises an error if it can find internal subject
# identifiers in the committed scripts.

success=true

for filename in "$@"; do
    if printf '%s' "$filename" | grep -E '2[1-5]F[0-9]{3,}|LA[0-9]{5,}' >/dev/null; then
        echo "Subject identifier found in a filename: '${filename}'"
        success=false
    fi
done

if grep -E -- '2[1-5]F[0-9]{3,}|LA[0-9]{5,}' "$@"; then
    success=false
fi

$success
