import pandas as pd
import glob
import sys
import os

# Specify the directory where your CSV files are located
input_dir = sys.argv[1]
input_dir.replace('\\', '/')
if input_dir[-1] != '/':
    input_dir += '/'

# Specify the pattern of CSV files you want to combine (e.g., all CSV files with .csv extension)
file_pattern = '*.csv'

# Get a list of CSV files in the directory that match the pattern
csv_files = glob.glob(input_dir + file_pattern)

# Initialize an empty DataFrame to store the combined data
combined_data = pd.DataFrame()

# Loop through each CSV file, read it, and concatenate it horizontally to the combined_data DataFrame
for csv_file in csv_files:
    data = pd.read_csv(csv_file, header=None)  # Assuming there are no headers in the input files
    combined_data = pd.concat([combined_data, data], axis=1)

# Save the combined data to a new CSV file
combined_data.to_csv('combined_output_tmp.csv', index=False, header=False)

# Prepend delimiter
with open('combined_output_tmp.csv', 'r') as tmp:
    with open('combined_output.csv', 'w') as outfile:
        outfile.write('sep=,\n')
        for line in tmp:
            outfile.write(line)

# Remove temporary file
os.remove("combined_output_tmp.csv")


print("Combined CSV file saved as 'combined_output.csv'")
