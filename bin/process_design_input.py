#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def process_design(file) :

    df = pd.read_csv(file, sep = ",")

    # make csv for fastq channel
    df['simple_name'] = df.apply(lambda row: os.path.basename(row['R1']).split(".")[0], axis = 1)
    df['group'] = df.apply(lambda row: f"{row['condition']}_{row['antibody']}", axis = 1)
    df['control_group'] = df.apply(lambda row: f"{row['control']}_{row['antibody']}", axis = 1)

    fastq = df[['simple_name', 'R1', 'R2']]

    # to join replicates
    replicates = df[['simple_name', 'group']]

    # once replicates are joined need to compare to control 
    print(df)
    controls = df[ df['control'].notna() ].reset_index(drop = True)
    controls = controls[['control_group', 'group']]

    # control_association
    fastq.to_csv("fastq.csv", sep = ',', header = True, index = False)
    replicates.to_csv("replicates.csv", sep = ',', header = True, index = False)
    controls.to_csv("controls.csv", sep = ",", header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-input", type = str, required = False)

    return parser.parse_args()

def main() : 

    args = get_args() 
    process_design(args.input)

if __name__ == "__main__" : 

    main()