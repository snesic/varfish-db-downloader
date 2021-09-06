import pandas as pd
import vcfpy
from os import walk

# -*- coding: utf-8 -*-

path = '/Users/snesic/sciebo/Development/WES/varfish/varfish-db-downloader/files'

def get_ref_from_filename(file):
    return {'filename' : file.split('/')[-1], 'reference' : file.split('/')[1]}

def check_vcf(file):
    # ------ read line by line
    # reader = vcfpy.Reader.from_path(file)
    # counts = dict()
    # for record in reader:
    #     i = record.CHROM
    #     counts[i] = counts.get(i, 0) + 1
    # counts
    # df = pd.DataFrame.from_dict(counts, orient='index').reset_index(level=0)
    #------use pandas----------
    compression = 'gzip' if file.endswith('gz') else None

    df = pd.read_csv(file, comment="#", sep='\t', compression=compression, header=None, usecols=[0])

    df.columns =['chr']
    df = df.chr.value_counts()
    df = df.reset_index(level=0)

    ##----add file name and output

    df.columns =['chr', 'count']
    df.loc[len(df.index)] = ['.', df['count'].sum()]
    df['file'] = file

    return df[['file', 'chr', 'count']]



def check_tsv(file):

    ## read TSV
    df = pd.read_csv(file, comment="#", sep='\t', low_memory=False)

    ## create a list of chromosomes to intersect with the columns
    chrs = list([str(x) for x in range(1,22)])
    chrs.extend(['X', 'Y'])
    chrs = pd.Series(chrs)
    ## take only unique values from each column and see which ones intersect with chrs
    a = df.apply(lambda x: x.unique(), axis=0)
    a = a.apply(lambda x: [str(i).lower().replace('chr', '').upper() for i in x])
    a = a.apply(lambda x: len(list(set(chrs) - set(x))))
    chr_columns = a[a==0].index.tolist()

    ## take only chr columns and counts frequency
    out = df.loc[:,chr_columns].apply(pd.Series.value_counts)
    ## Pivot from wide to long format
    out.reset_index(level=0, inplace=True)
    out = pd.melt(out, id_vars=['index'], value_vars=chr_columns)
    out.columns = ['chr', 'column_name', 'count']

    # sum all chromosomes
    out_all = out.groupby('column_name').agg({'count' : 'sum'})
    out_all['chr'] = '.'
    out_all.reset_index(level=0, inplace=True)
    out = pd.concat([out, out_all])

    ## add filename
    out['filename'] = file

    return out.loc[:, ['filename', 'chr', 'column_name', 'count']]

# ---- Collect all the files

files = []
for (dirpath, dirnames, filenames) in walk(path):
    files.extend([dirpath + '/' + i for i in filenames])


# ---- Check VCF files ------
vcfs = [i for i in f if i.endswith('vcf.gz') or i.endswith('vcf.bgz')]

df = pd.concat([check_vcf(i) for i in vcfs])


# ---- Check TSV files ------
tsvs = [i for i in files if i.endswith('.tsv')]

df = pd.concat([check_tsv(i) for i in tsvs])

# ---- save output ----
df.to_csv('sanity_check_summary.tsv', sep='\t', index=False, mode='w')


#### ----- main - once everything finish copy the code here!

#def main(argv=None):
#    check_vcf(file)

#if __name__ == "__main__":
#    sys.exit(main())
