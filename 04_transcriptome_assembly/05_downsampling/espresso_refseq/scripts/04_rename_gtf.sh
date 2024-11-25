#!/bin/bash

for filepath in data/espresso_q/*/samples_N2_R0_updated.gtf; do
    # Extract the directory name (the part matching '*')
    dirname=$(basename $(dirname "$filepath"))
    
    # Construct the new filename
    new_filepath="data/espresso_q/$dirname/${dirname}.gtf"
    
    # Rename the file
    mv "$filepath" "$new_filepath"
done
