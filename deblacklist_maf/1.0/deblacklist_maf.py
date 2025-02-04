#!/usr/bin/env python

"""
Remove blacklisted variants from MAF files

This script was created by adapting parts of the `annotate_ssm_blacklist()` GAMBLR function
to support existing blacklist scripts that exist within GAMBL. 
Generally, the minimum column required is `chrpos`
Ideally, columns should include `Reference_Allele` and `Tumor_Seq_Allele2` for the reference and alt alleles
All "filterable" columns are: `chrpos` `Reference_Allele` `Tumor_Seq_Allele2`
"""

import sys
import logging
import argparse
import pandas as pd
import csv

# Set up logging

logger = logging.getLogger("deblacklist_logger")
stdout_h = logging.StreamHandler(sys.stdout)
stderr_h = logging.StreamHandler()
stdout_h.setLevel(logging.INFO)
stderr_h.setLevel(logging.ERROR)
stdout_h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
stdout_h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

logger.addHandler(stdout_h)
logger.addHandler(stderr_h)

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Input MAF file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Deblacklisted output MAF')
    parser.add_argument('-d','--drop-threshold',
                        type=int, default=4,
                        help='Minimum number of blacklist occurrences required to drop mutation')
    parser.add_argument('-b', '--blacklists',
                        type=str, nargs='*',
                        help='Blacklist tables. Can provide multiple filepaths separated by spaces')

    return parser.parse_args()

def deblacklist(maf, blacklist, drop_threshold, minimum_cols, additional_cols):

    logger.warning(f"Applying {blacklist} blacklist...")
    logger.warning(f"MAF file has {len(maf)} rows prior to blacklist variant removal")

    blacklist_df = pd.read_table(blacklist, comment="#", sep="\t")

    cols = list(blacklist_df.columns)

    # Skip blacklisting if minimum columns not present
    if not set(minimum_cols).issubset(cols):
        logger.error(f"Minimum columns {' '.join(minimum_cols)} columns not found in blacklist file {blacklist}")
        logger.warning(f"Skipping blacklist file {blacklist}")
        return maf

    blacklist_df["Chromosome"] = blacklist_df["chrpos"].apply(lambda x: x.split(":")[0])
    blacklist_df["Start_Position"] = blacklist_df["chrpos"].apply(lambda x: int(x.split(":")[1]))

    filter_cols = ["Chromosome","Start_Position"] + additional_cols if set(additional_cols).issubset(cols) else ["Chromosome","Start_Position"]

    logger.warning(f"Filtering by columns: {' '.join(filter_cols)}")

    # Get unique combinations in blacklist
    logger.warning(f"{blacklist}: {len(blacklist_df)} rows in blacklist.")

    values = set(map(tuple, blacklist_df[filter_cols].values))

    if not values:
        variants_dropped = 0
    else:
        blacklisted_rows = maf[filter_cols].apply(lambda x, values: tuple(x) in values, args=(values,), axis=1)
        variants_dropped = sum(blacklisted_rows)
        maf = maf[~blacklisted_rows]
    
    logger.warning(f"{variants_dropped} variants were dropped")
    logger.warning(f"MAF file contains {len(maf)} rows after blacklisting")

    return maf

def main(args):
    maf = pd.read_table(args.input, comment="#", sep="\t")

    minimum_cols = ["chrpos"]
    additional_cols = ["Reference_Allele","Tumor_Seq_Allele2"]

    for blacklist in args.blacklists:
        maf = deblacklist(maf, blacklist, args.drop_threshold, minimum_cols, additional_cols)

    maf.to_csv(args.output, sep="\t", na_rep="NA", index=False)

if __name__ == '__main__':

    args = parse_arguments()
    main(args)
