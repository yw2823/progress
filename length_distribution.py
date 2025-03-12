import sys
from collections import defaultdict


def main(argv):
    """Process input sequences and count read lengths."""
    
    if len(argv) < 4:
        print(f'Usage: python {argv[0]} inputfilename outfilename [-f | -q]')
        print('\tUse - for stdin instead of an input filename')
        sys.exit(1)

    input_filename = argv[1]
    output_filename = argv[2]
    seq_type = argv[3]

    # Open output file
    with open(output_filename, 'w') as outfile:
        # Handle input stream (file or stdin)
        input_stream = sys.stdin if input_filename == '-' else open(input_filename)

        read_count = 0
        length_counts = defaultdict(int)  # Dictionary to store length counts

        # Process FASTA files
        if seq_type == '-f':
            seq = ""
            for line in input_stream:
                if read_count % 1_000_000 == 0 and read_count > 0:
                    print(f'{read_count} reads processed')

                if line.startswith(">"):  # Header line
                    if seq:  # Process previous sequence
                        length_counts[len(seq)] += 1
                    seq = ""  # Reset sequence
                else:
                    seq += line.strip()
                read_count += 1

            # Process last sequence
            if seq:
                length_counts[len(seq)] += 1

        # Process FASTQ files
        elif seq_type == '-q':
            skip_next = False
            seq_next = False

            for line in input_stream:
                if line.startswith('+'):  
                    skip_next = True  # Skip quality score header
                    continue
                if skip_next:
                    skip_next = False  # Skip quality score line
                    continue
                if line.startswith('@'):  
                    seq_next = True  # Sequence follows
                    continue
                if seq_next:
                    length_counts[len(line.strip())] += 1  # Count sequence length
                    seq_next = False
                    read_count += 1
                    if read_count % 1_000_000 == 0:
                        print(f'{read_count} reads processed')

        # Write length distribution to output file
        for length in sorted(length_counts.keys()):
            outfile.write(f'{length}\t{length_counts[length]}\n')

if __name__ == '__main__':
    main(sys.argv)