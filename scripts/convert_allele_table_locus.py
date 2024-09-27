#!/usr/bin/env python3

import pandas as pd
import argparse


def correct_pool(pool_value):
    pool_corrections = {
        '1AB': '1A,1B',
        '1B2': '1B,2'
    }
    return pool_corrections.get(pool_value, pool_value)


def convert_locus_and_pool(input_file, output_file):
    # Load the old output
    df_old = pd.read_csv(input_file, sep='\t')

    # Split the 'Locus' column by the last hyphen to separate the Pool
    df_old[['Locus', 'Pool']] = df_old['Locus'].str.rsplit(
        '-', n=1, expand=True)

    # Correct specific pool values
    df_old['Pool'] = df_old['Pool'].apply(correct_pool)

    df_old.drop(columns=['Allele'], inplace=True)

    # Save the updated DataFrame to the new format
    df_old.to_csv(output_file, sep='\t', index=False)


def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description='Convert old output to new format with Locus and Pool separated.')
    parser.add_argument('input_file', help='Path to the old output file')
    parser.add_argument('output_file', help='Path to save the new output file')

    args = parser.parse_args()

    # Convert old output to new format
    convert_locus_and_pool(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
