#!/bin/bash

# Redirect both stdout and stderr to a file
exec > extracted_params.h 2>&1

# Input: fitall.log
LOGFILE="fitall.log"

# Check if file exists
if [[ ! -f "$LOGFILE" ]]; then
  echo "Error: $LOGFILE not found!"
  exit 1
fi

# Function to extract parameters and format as C++ array for a specific energy
extract_param() {
  local energy="$1"     # Energy to filter
  local param="$2"      # Parameter to search (e.g., "par2")
  local label="$3"      # Label for output
  local asym_label="$4" # Asymmetric error label
  local values=()
  local err_low=()
  local err_high=()
  local processing_energy=0

  while IFS= read -r line; do
    # Check for the start of the relevant energy block
    if [[ $line =~ Running\ processing\ and\ fit\ for\ energy:\ $energy ]]; then
      processing_energy=1
      continue
    fi

    # Stop processing when another energy block begins
    if [[ $processing_energy -eq 1 && $line =~ Running\ processing\ and\ fit\ for\ energy:\  ]]; then
      break
    fi

    # Extract parameter values
    if [[ $processing_energy -eq 1 && $line =~ $param=([^+]+) ]]; then
      values+=("${BASH_REMATCH[1]}")
    fi

    # Extract asymmetric errors
    if [[ $processing_energy -eq 1 && $line =~ $asym_label:\ \+([^[:space:]]+)\ --([^[:space:]]+) ]]; then
      err_high+=("${BASH_REMATCH[1]}")
      err_low+=("${BASH_REMATCH[2]}")
    fi
  done < "$LOGFILE"

  # Join values with commas for C++ arrays
  local values_str=$(IFS=,; echo "${values[*]}")
  local err_low_str=$(IFS=,; echo "${err_low[*]}")
  local err_high_str=$(IFS=,; echo "${err_high[*]}")

  # Print values as C++ array
  echo "double ${label}Values[10] = {$values_str};" # 5 or any other number of K_T bins
  echo "double ${label}ErrLow[10] = {$err_low_str};"
  echo "double ${label}ErrHigh[10] = {$err_high_str};"
  echo
}

# Energy to extract
ENERGY="$1"

if [[ -z "$ENERGY" ]]; then
  echo "Usage: $0 <energy>"
  echo "Example: $0 9.2GeV"
  exit 1
fi

# Extract and format data for alpha, N, and R
echo "// C++ Arrays for ROOT script for energy: $ENERGY"
extract_param "$ENERGY" "par2" "alpha" "Asym. err alpha"
extract_param "$ENERGY" "par0" "N" "Asym. err N"
extract_param "$ENERGY" "par1" "R" "Asym. err R"