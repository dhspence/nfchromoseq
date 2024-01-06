#!/usr/bin/env python3

import argparse
import pandas as pd
import sys, json

def main():
    # Set up argparse to handle the -i option and multiple CSV file arguments
    parser = argparse.ArgumentParser(description='Process CSV files')
    parser.add_argument('-i', '--id', required=True, help='ID parameter')
    parser.add_argument('-u', '--uniquify', action='store_true', default=False, help='only keep unique items')
    parser.add_argument('csv_files', nargs='+', help='CSV file paths')
    args = parser.parse_args()

    # Read and concatenate CSV files
    dataframes = [pd.read_csv(csv_file) for csv_file in args.csv_files]
    concatenated_df = pd.concat(dataframes, ignore_index=True)

    # Convert the concatenated DataFrame to a dictionary
    data_dict = concatenated_df.to_dict(orient='list')

    if (args.uniquify):
        data_dict = {key: list(set(value)) for key, value in data_dict.items()}

    output_file_name = "runinfo.json"
    with open(output_file_name, 'w') as file:
        json.dump(data_dict, file, indent=4)

if __name__ == '__main__':
    main()
