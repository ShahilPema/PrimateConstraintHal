#!/usr/bin/env python3

import sys

def validate_maf(filename):
    try:
        with open(filename, 'r') as file:
            in_block = False
            line_num = 0
            for line in file:
                line_num += 1
                line = line.strip()
                if not line:
                    in_block = False
                    continue
                if line.startswith('a'):
                    if in_block:
                        print(f"Error: Unexpected 'a' line at {line_num}")
                        return False
                    in_block = True
                elif line.startswith('s'):
                    if not in_block:
                        print(f"Error: 's' line outside of block at {line_num}")
                        return False
                    parts = line.split()
                    if len(parts) != 7:
                        print(f"Error: Incorrect number of fields in 's' line at {line_num}")
                        return False
                    try:
                        int(parts[2])
                        int(parts[3])
                        int(parts[5])
                    except ValueError:
                        print(f"Error: Non-integer values in 's' line at {line_num}")
                        return False
        print("MAF file appears to be valid.")
        return True
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate_maf.py <maf_file>")
        sys.exit(1)
    else:
        result = validate_maf(sys.argv[1])
        if result:
            sys.exit(0)  # Successful validation
        else:
            sys.exit(1)  # Validation failed

