#!/bin/bash

# Check if a directory is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 directory [output_file]"
  exit 1
fi

# Get the absolute path of the directory
directory=$(cd "$1" && pwd)

# Define the output file (default is "file_paths.txt" if not provided)
output_file=${2:-file_paths.txt}

# Find all files and print their absolute paths to the output file
find "$directory" -type f -print > "$output_file"

echo "File paths have been saved to $output_file"

