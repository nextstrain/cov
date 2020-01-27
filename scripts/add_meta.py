import pandas as pd
import numpy as np
import re

# This script adds a manually-curated metadata file (from info that can't
# be scraped from Genbank) to the metadata file.


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='add extra metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--extra-meta-in', help="input extra meta tsv")           
    parser.add_argument('--meta-in', help="input metadata file")
    parser.add_argument('--meta-out', help="final output metadata with age categories too")
    args = parser.parse_args()

    #################################################################
    ##      Part 1
    #################################################################

    new_meta = pd.read_csv(args.extra_meta_in, sep='\t', index_col=False)
    meta = pd.read_csv(args.meta_in, sep='\t', index_col=False)

    # Uncomment and modify if need new column
    #create new column for AFM (& reorder)
    #meta = meta.reindex(columns = ['strain', 'accession', 'date', 'sex', 'age', 
    #    'symptom', 'country', 'collab_country', 'region', 'host',
    #    'subgenogroup', 'Lab-ID', 'orig_strain', 'seq-len', 'date_added', 'genbank-host', 'moltype',
    #    'virus', 'authors', 'title', 'url', 'paper_url'])


    for i, row in new_meta.iterrows():
        #only do this if the accession is in the metadata
        if row.accession in meta.accession.values:
            
            if not pd.isnull(row.host):
                meta.loc[meta.accession == row.accession, 'host'] = row.host.strip()
            if not pd.isnull(row.date):
                meta.loc[meta.accession == row.accession, 'date'] = row.date.strip()
            #uncomment & change to add more fields
            #if not pd.isnull(row.sex):
            #    meta.loc[meta.accession == row.accession, 'sex'] = row.sex.strip()
            #if not pd.isnull(row.symptom):
            #    meta.loc[meta.accession == row.accession, 'symptom'] = row.symptom.strip()


    meta.to_csv(args.meta_out, sep='\t', index=False)
