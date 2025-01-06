#!/bin/zsh

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input-file.tex>"
  exit 1
fi

# Input and output file names
input_file="$1"
output_file="${input_file%.tex}.txt"

# Apply sed commands
sed -e 's/\\begin{align}/\\f[/g' \
    -e 's/\\end{align}/\\f]/g' \
    -e 's/\\v/\\vec/g' \
    -e 's/\\diffp{/\\frac{\\partial /g' \
    -e 's/\\//g' \
    -e 's/\nonumber//g' \
    -e 's/&//g' \
    -e 's/\$/\\f\$/g' \
    "$input_file" > "$output_file"

echo "Conversion complete. Output saved to $output_file"
