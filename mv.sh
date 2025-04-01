#!/bin/bash

# Iterate over directories matching the pattern
for dir in $(ls -d work/*/*/*final.bed); do
    # Determine target directory (parent of the *final.bed directory)
    target_dir=$(dirname "$dir")

    # Move all contents from current *final.bed directory up one level
    mv "$dir"/* "$target_dir"/

    echo "Moved contents of $dir to $target_dir"
done
