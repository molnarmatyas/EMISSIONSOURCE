#!/bin/bash

# Step 1: Rename files to a temporary name
for i in {1..100}; do
  old_file="z-prueba_EOS2_${i}.root"
  temp_file="temp_${i}.root"
  mv "$old_file" "$temp_file"
done

# Step 2: Rename files from the temporary name to the final desired name
for i in {1..100}; do
  temp_file="temp_${i}.root"
  new_file="z-prueba_EOS2_$((i + 80)).root"
  mv "$temp_file" "$new_file"
done
