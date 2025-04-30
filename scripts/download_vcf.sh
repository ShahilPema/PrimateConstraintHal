#!/bin/bash

# Attempt to download the file using bs command
bs  download file -i $1 -o .

find . -type f -name "*.gz" -exec sh -c '
    mv "$1" .
    dir=$(dirname "$1")
    if [ "$dir" != "." ] && [ -z "$(ls -A "$dir")" ]; then
        rmdir "$dir"
    fi
' sh {} \;

rm *.json
