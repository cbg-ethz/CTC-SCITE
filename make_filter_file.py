import csv

def parse_and_write_tsv(input_file_path, output_file_path):
    # Open the input TSV file for reading
    with open(input_file_path, 'r', newline='', encoding='utf-8') as input_tsv:
        # Create a CSV reader for the input TSV file
        reader = csv.reader(input_tsv, delimiter='\t')

        # Skip the header
        header = next(reader, None)

        # Find the index of the last column
        last_column_index = len(header) - 1

        # Collect rows with "CTC" in the last column
        selected_rows = [row for row in reader]# if (row[last_column_index] == "CTC" or row[last_column_index] == "WBC" or row[last_column_index] == "tumor_cell")]

    # Open the output TSV file for writing
    with open(output_file_path, 'w', newline='', encoding='utf-8') as output_tsv:
        # Create a CSV writer for the output TSV file
        writer = csv.writer(output_tsv, delimiter='\t')

        # Write the header to the output file
        writer.writerow(header[:1])  # Write only the first column as header

        # Write selected first column values to the output file
        for row in selected_rows:
            writer.writerow([row[0]])

# Replace 'input_file.tsv' and 'dna_samples_nofilter.tsv' with the actual paths to your input and output TSV files
input_file_path = 'gDNA_sample_annotation_final.tsv'
output_file_path = 'dna_samples_nofilter.tsv'
parse_and_write_tsv(input_file_path, output_file_path)

