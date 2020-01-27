import pandas as pd
from dateutil.parser import parse
import re, os
import shutil


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='find new sequences',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input-new', help="new meta file from which to exclude duplicates")
    parser.add_argument('--exclude', nargs='+', help="meta files whose accession #s to exclude")
    parser.add_argument('--output', help="meta with duplicates excluded")
    args = parser.parse_args()


    new_meta = pd.read_csv(args.input_new, sep="\t", index_col=False)
    newM = set(new_meta.accession)

    for excl in args.exclude:
        if os.path.isfile(excl):
            excl_meta = pd.read_csv(excl, sep='\t' if excl[-3:]=='tsv' else ',', index_col=False)
            newM = newM - set(excl_meta.accession)
        else:
            print("This is a first run - not comparing to current database.")

    new_actual_meta = new_meta.loc[new_meta['accession'].isin(newM)]

    print("{} new accession numbers were found and will be downloaded.".format(len(newM)))

    if not newM:
        print("No new accession numbers were added. Run will not proceed!")
    else:
        new_actual_meta.to_csv(args.output, sep='\t', index=False)
        


