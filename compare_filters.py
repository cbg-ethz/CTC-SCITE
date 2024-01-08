def compare_files(file1_path, file2_path):
    # Read the content of the first file and split it into words
    with open(file1_path, 'r', encoding='utf-8') as file1:
        words_file1 = set(file1.read().split())

    # Read the content of the second file and split it into words
    with open(file2_path, 'r', encoding='utf-8') as file2:
        words_file2 = set(file2.read().split())

    # Find the words that are in file1 but not in file2
    unique_to_file1 = words_file1 - words_file2

    # Find the words that are in file2 but not in file1
    unique_to_file2 = words_file2 - words_file1

    return unique_to_file1, unique_to_file2

# Replace 'file1.txt' and 'file2.txt' with the actual paths to your files
file1_path = 'dna_samples_nofilter.txt'
file2_path = 'Br16_notAllMissingValues.txt'

unique_words_file1, unique_words_file2 = compare_files(file1_path, file2_path)

print("Words unique to file1:", unique_words_file1)
print("Words unique to file2:", unique_words_file2)

