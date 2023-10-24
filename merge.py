import pandas as pd
import glob
import sys

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
combined_data.to_csv('combined_output.csv', index=False, header=False)

# Prepend delimitor
with open('combined_output.csv', 'r+') as f:
    content = f.read()
    f.seek(0, 0)
    f.write('sep=,\n' + content)

print("Combined CSV file saved as 'combined_output.csv'")
