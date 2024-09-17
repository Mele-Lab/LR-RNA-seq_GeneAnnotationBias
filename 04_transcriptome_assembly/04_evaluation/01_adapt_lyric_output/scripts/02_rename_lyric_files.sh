#!/bin/bash

#!/bin/bash

# Loop through the files provided as arguments or if no argument, process all matching files
for file in "$@"; do
  # Extract the relevant parts using a more precise regex
  if [[ $file =~ ont_HpreCap_0\+\_([0-9]+)([A-Z]+[0-9]+)([A-Z]+[0-9]+)\..*gff_updated.gff.gz ]]; then
    # Construct the new filename with "_updated"
    new_name="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}_updated.gff.gz"
    # Rename the file
    mv "$file" "$new_name"
  elif [[ $file =~ ont_HpreCap_0\+\_([0-9]+)([A-Z]+[0-9]+)([A-Z]+[0-9]+)\..*gff.gz ]]; then
    # Construct the new filename without "_updated"
    new_name="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}.gff.gz"
    # Rename the file
    mv "$file" "$new_name"
  fi
done

