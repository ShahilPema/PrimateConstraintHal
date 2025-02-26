import os
import sys

def split_bed_file(input_file, max_lines=100000):
    basename = os.path.splitext(os.path.basename(input_file))[0]
    with open(input_file, 'r') as bed:
        file_number = 0
        line_count = 0
        current_outfile = open(f"{basename}_{file_number}.bed", 'w')

        for line in bed:
            if line_count == max_lines:
                current_outfile.close()
                file_number += 1
                current_outfile = open(f"{basename}_{file_number}.bed", 'w')
                line_count = 0

            current_outfile.write(line)
            line_count += 1

        current_outfile.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_bed_file.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    split_bed_file(input_file)
