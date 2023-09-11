#!/bin/bash
# Get the directory where the script resides
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")")"

# Change the working directory to the script's directory
cd "$script_dir"

image_name="noetic_srrg"

# docker rmi $image_name
docker build -t $image_name .