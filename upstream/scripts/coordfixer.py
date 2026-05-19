import pandas as pd
import sys
import os



# Input chromosome (e.g., "chr1") from command line argument

if len(sys.argv) != 3:
    print("Usage: python script.py <chromosome>")
    sys.exit(1)

chromosome = sys.argv[1]
OUTPUT_FOLDER = sys.argv[2]
# Define the directory paths and filenames based on the input chromosome
directory_path = OUTPUT_FOLDER  + '/modkit_splitbam_referenced_bychrandhap_filtered/'
directory_path_output = OUTPUT_FOLDER  + '/modkit_referenced_splitbam_mergedbychr/'
print(directory_path)
print(chromosome)
file_H1_name = os.path.join(directory_path, f'modkit_extract_{chromosome}_H1_filtered.tsv')
file_H2_name = os.path.join(directory_path, f'modkit_extract_{chromosome}_H2_filtered.tsv')
file_noH_name = os.path.join(directory_path, f'modkit_extract_{chromosome}_noH_filtered.tsv')
output_file_name = os.path.join(directory_path_output, f'modkit_extract_{chromosome}_merged_filtered_coordmod.tsv')

# Ensure the output directory exists
os.makedirs(directory_path_output, exist_ok=True)

# Load the files with selected columns only (1, 3, 4, 6, 11)
cols_to_keep = [0, 2, 3, 5, 10]  # Note: 0-indexed, corresponding to columns 1, 3, 4, 6 in pandas
file_H1 = pd.read_csv(file_H1_name, sep='\t', header=None, usecols=cols_to_keep)
file_H2 = pd.read_csv(file_H2_name, sep='\t', header=None, usecols=cols_to_keep)
file_noH = pd.read_csv(file_noH_name, sep='\t', header=None, usecols=cols_to_keep)

# Adjust coordinates for column 3 (index 2) if column 6 (index 5) is "-"
file_H1[2] = file_H1.apply(lambda row: row[2] - 1 if row[5] == "-" else row[2], axis=1)
file_H2[2] = file_H2.apply(lambda row: row[2] - 1 if row[5] == "-" else row[2], axis=1)
file_noH[2] = file_noH.apply(lambda row: row[2] - 1 if row[5] == "-" else row[2], axis=1)

# Add the Haplotag column for each file
file_H1['Haplotag'] = 'H1'
file_H2['Haplotag'] = 'H2'
file_noH['Haplotag'] = 'noH'

# Keep only columns 1, 3, 4, and the Haplotag column (final structure: 1, 3, 4, Haplotag)
file_H1 = file_H1[[0, 2, 3, 10, 'Haplotag']]
file_H2 = file_H2[[0, 2, 3, 10, 'Haplotag']]
file_noH = file_noH[[0, 2, 3, 10,  'Haplotag']]

# Merge the three files into one
merged_file = pd.concat([file_H1, file_H2, file_noH])

# Sort the merged file by column 1 (alphabetically) and column 4 (numerically)
#merged_file_sorted = merged_file.sort_values(by=[0, 2])

# Save the final sorted merged file
merged_file.to_csv(output_file_name, sep='\t', header=False, index=False)

# Function to find min and max coordinates
def find_min_max(df):
    return df[2].min(), df[2].max()

# Find min and max coordinates for the final merged file
min_max_coords = find_min_max(merged_file)

# Print the min and max coordinates for the final merged file
print(f"{chromosome} Merged Min and Max:", min_max_coords)
