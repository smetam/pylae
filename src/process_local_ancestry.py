# Usage:
# python3 process_local_ancestry.py <filename1> <filename2> ...

# Supposed file layout:
# CHR          SNP     CLST   A1   A2      MAF    MAC  NCHROBS
#   1   rs12562034  Chane_0    1    3        0      0        2
#   1   rs12562034 Mediterranean    1    3  0.06667      2       30
#   1   rs12562034 NativeAmerican    1    3      0.1      3       30
#   1   rs12562034 NorthEastAsian    1    3   0.5333     16       30
#   1   rs12562034 NorthernEuropean    1    3        0      0       30
#   1   rs12562034 Oceanian    1    3   0.2667      8       30
#   1   rs12562034 SouthAfrican    1    3   0.4333     13       30
#   1   rs12562034 SouthEastAsian    1    3   0.3667     11       30
#   1   rs12562034 SouthWestAsian    1    3   0.2667      8       30
#   1   rs12562034 SubsaharanAfrican    1    3  0.03333      1       30

import time
import argparse
import pandas as pd
import numpy as np

from operator import itemgetter
from concurrent.futures import ThreadPoolExecutor


def get_sample_name(file):
    sample = file[file.rfind('/'):]
    return sample[sample.find('.') + 1:sample.find('.txt')]


def fast_best_pop(db, maf):
    diff = zip(db['CLST'], abs(db['MAF'] - maf))
    best_pop, d = min(diff, key=itemgetter(1))
    return best_pop


def df_best_pop(db, maf):
    maf_df = pd.Series(np.abs(db['MAF'].values - maf), index=db['CLST'])
    best_pop = maf_df.idxmin(axis=0, skipna=True)
    return best_pop


def main(files, output):
    start_time = time.monotonic()

    for file in files:
        sample = get_sample_name(file)
        print(f'Running on group: {sample}')

        df_it = pd.read_csv(file, chunksize=10, sep='\s+')
        snp_i = 0
        with open(output + sample + '.local_an.txt', 'w') as f:
            for df in df_it:
                query = df[df['CLST'] == sample].iloc[0]
                snp = query['SNP']
                maf = query['MAF']
                db = df[df['CLST'] != sample]

                best_pop = fast_best_pop(db, maf)
                # best_pop = best_pop_like_r(db, maf)
                # best_pop = df_best_pop(db, maf)
                f.write(f'{snp}|{best_pop}|{best_pop}\n')
                # snp_i += 1
                # if snp_i > 1000:
                #     break

    print(f'Finished in {time.monotonic() - start_time} sec.')


def process_file(file):
    output = '../data/'
    sample = get_sample_name(file)
    print(f'Running on group: {sample}')

    df_it = pd.read_csv(file, chunksize=10, sep='\s+')
    with open(output + sample + '.local_an.txt', 'w') as f:
        for df in df_it:
            query = df[df['CLST'] == sample].iloc[0]
            snp = query['SNP']
            maf = query['MAF']
            db = df[df['CLST'] != sample]

            best_pop = fast_best_pop(db, maf)
            f.write(f'{snp}|{best_pop}|{best_pop}\n')


def concurrent_files_processing(files, output):
    start_time = time.monotonic()

    with ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(process_file, files)

    print(f'Finished in {time.monotonic() - start_time} sec.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Input filenames", nargs='+')
    parser.add_argument("-o", "--output", help="Output directory", default='../data/')
    args = parser.parse_args()

    if len(args.files) > 0:
        main(args.files, args.output)
    else:
        print('No filenames provided, please specify files to run on.')
