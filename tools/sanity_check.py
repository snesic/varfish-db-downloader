import pandas as pd
import vcfpy
from os import walk

# -*- coding: utf-8 -*-

path = '/Users/snesic/sciebo/Development/WES/varfish/varfish-db-downloader/files'

def extract_reference_from_filename(file):

    ref = [i for i in ['GRCh37', 'GRCh38'] if file.lower().find(i.lower()) > 0]
    return 'No Ref' if len(ref) == 0 else ', '.join(ref)

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
    df['Reference'] = extract_reference_from_filename(file)

    return df[['file', 'Reference', 'chr', 'count']]



# read a file and take only unique values
# for each column
# create a long format for easier manipulation

def file_to_dataframe(file):

    print('Process: ' + file)
    df = pd.read_csv(file, comment="#", sep='\t', low_memory=False)
    ##TODO: check if pandas handles na values properly, maybe redefine initial na_values list..
    # create a dataframe where each row is one column name and its type
    df_types = pd.DataFrame(df.dtypes).reset_index()
    df_types.columns = ['col_name', 'type']

    # in case of numerical columns, calculate MIN and MAX
    df_stat = df.select_dtypes(exclude=['object'])
    if df_stat.shape[1] == 0:
        df_stat = pd.DataFrame(columns=['col_name', 'stat'])
    else:
        df_stat = df_stat.apply(lambda x: str(x.dropna().min()) + ', ' + str(x.dropna().max()))
        df_stat = pd.DataFrame(df_stat).reset_index()
        df_stat.columns = ['col_name', 'stat']

    # Pivot columns to long format to count the elements per column
    df.reset_index(level=0, inplace=True)
    df = pd.melt(df, id_vars=['index'], var_name='col_name', value_name='col_item').drop(columns='index')
    df.columns = ['col_name', 'col_item']

    # replace continuous values with "continuous"
    df_unique = df.drop_duplicates(['col_name', 'col_item'], keep='first')
    df_unique = df_unique.groupby('col_name').size().to_frame('unique').reset_index()
    df_unique['continuous'] = df_unique.unique > 50

    continuous_columns = df_unique[df_unique.continuous].col_name
    df.loc[df.col_name.isin(continuous_columns), 'col_item'] = 'continuous' # rename continuous variables


    # count occurances of each element of each column
    out = df.groupby(['col_name', 'col_item']).size().to_frame('count').reset_index()

    # add column type and column stat(in case of numeric)
    out = out.merge(df_types, on='col_name', how='left')
    out = out.merge(df_stat, on='col_name', how='left')

    out.col_item = out.col_item.astype(str)
    out['filename'] = file
    out['Reference'] = extract_reference_from_filename(file)

    return out



# in case files are read in chunks
# append all of them
def merge_tsv_dataframes(list_df):

    df = pd.concat(list_df)

    df = df.groupby(['filename', 'Reference', 'col_name', 'col_item'])
    df = df.aggregate({'count' : sum,
                       'type' : lambda x: ', '.join(x.astype('str').unique()),
                       'stat' : lambda x: ', '.join(x.fillna(''))}).reset_index()

    return df



def chr_summary_dataframe(df):

    ## create a list of chromosomes to intersect with the columns
    chrs = list([str(x) for x in range(1,22)])
    chrs.extend(['x', 'y'])
    chrs = pd.Series(chrs)

    # intersect values of each column with chrs
    df_chr = df[df.col_item != 'continuous']
    df_chr = df_chr.groupby(['filename', 'Reference', 'col_name'])
    df_chr = df_chr.aggregate({'col_item' : lambda x: x.str.lower().replace('chr', '').unique()}).reset_index()
    df_chr = df_chr[~df_chr.col_item.isna()] # boolinas are set to na after converstion to string

    chr_cols = df_chr.col_item.apply(lambda x: len(list(set(chrs) - set(x)))) == 0

    chr_cols = df_chr[chr_cols].col_name

    known_chr_columns = ['chromosome']    # add columns that are known to represent chromosomes
    chr_cols = chr_cols.append(pd.Series(known_chr_columns))

    out = df[df.col_name.isin(chr_cols)]
    out = out.drop(columns=['type', 'stat'])
    out.columns = ['filename', 'Reference', 'col_name', 'col_item', 'count']
    out
    return out



def group_values_per_column(df):

    dd = df.groupby(['filename', 'Reference', 'col_name', 'type', 'stat'])
    dd = dd.aggregate({'col_item' : lambda x: ', '.join(x.str.replace('"""', '').str.strip('{}').str.split(r',|\|').explode().unique())}).reset_index()
    return dd[['filename', 'Reference', 'col_name', 'col_item', 'type', 'stat']]




# ---- Collect all the files

files = []
for (dirpath, dirnames, filenames) in walk(path):
    files.extend([dirpath + '/' + i for i in filenames])

# ---- Check VCF files ------
#vcfs = [i for i in files if i.endswith('vcf.gz') or i.endswith('vcf.bgz')]
#df = pd.concat([check_vcf(i) for i in vcfs])
#df


# ---- Check TSV files - all columns info ------
tsvs = [i for i in files if i.endswith('.tsv')]



list_df = [file_to_dataframe(i) for i in tsvs]

df = merge_tsv_dataframes(list_df)
# only chromosome summary
df_chr = chr_summary_dataframe(df)

# General column summary
df = group_values_per_column(df)


# ---- save output ----
df_chr.to_csv('sanity_check_chr_summary.tsv', sep='\t', index=False, mode='w')
df.to_csv('sanity_check_general_summary.tsv', sep='\t', index=False, mode='w')


#### ----- main - once everything finish copy the code here!

#def main(argv=None):
#    check_vcf(file)

#if __name__ == "__main__":
#    sys.exit(main())
