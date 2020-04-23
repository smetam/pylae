# Usage:
# python3 process_local_ancestry.py <filename1> <filename2> ...

# Supposed file layout:
# CHR          SNP     CLST   A1   A2      MAF    MAC  NCHROBS
#   1   rs12562034 PopGroup    1    3        0      0        2
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

from scipy.special import softmax
from collections import Counter
from operator import itemgetter
from concurrent.futures import ThreadPoolExecutor, as_completed


OUTPUT_DIR = '../data/'
POPS = [
    'Mediterranean', 'NativeAmerican', 'NorthEastAsian', 'NorthernEuropean',
    'Oceanian', 'SouthAfrican', 'SouthEastAsian', 'SouthWestAsian', 'SubsaharanAfrican'
]


def get_group_name(file):
    sample = file[file.rfind('/'):]
    return sample[sample.find('.') + 1:sample.find('.txt')]


def fast_best_pop(db, maf):
    diff = zip(db['CLST'], abs(db['MAF'] - maf))
    best_pop, d = min(diff, key=itemgetter(1))
    return best_pop


def process_file_old(file):
    output = '../data/'
    group = get_group_name(file)
    print(f'Running on group: {group}')

    df_it = pd.read_csv(file, chunksize=10, sep='\s+')
    with open(output + group + '.local_an.txt', 'w') as f:
        for df in df_it:
            query = df[df['CLST'] == group].iloc[0]
            snp = query['SNP']
            maf = query['MAF']
            db = df[df['CLST'] != group]

            best_pop = fast_best_pop(db, maf)
            f.write(f'{snp}|{best_pop}|{best_pop}\n')


def softmax_prob(db, maf, scale=5):
    """Softmax method for probabilities"""
    db = pd.Series(np.abs(db['MAF'].values - maf), index=db['CLST'].values)
    return np.around(softmax(scale * (1 - db.values)), 2)


def process_file(file):
    group = get_group_name(file)
    output = OUTPUT_DIR + group + '.local_an.txt'
    print(f'Running on group: {group}\n')

    with open(output, 'w') as f:
        f.write('SNP_ID\t' + '\t'.join(POPS) + '\n')
        for df in pd.read_csv(file, chunksize=10, sep='\s+'):
            query = df[df['CLST'] == group].iloc[0]
            snp = query['SNP']
            maf = query['MAF']
            db = df[df['CLST'] != group]

            pop_probabilities = softmax_prob(db, maf)
            f.write(f'{snp}\t' + '\t'.join(map(str, pop_probabilities)) + '\n')

    return output


def process_probabilities(file, window_len=300):
    output = f'{file}_{window_len}.csv'
    result = []

    with open(output, 'w') as f_out:
        for df in pd.read_csv(file, sep='\t', index_col=0, chunksize=window_len):
            s = pd.Series([df[pop].mean() for pop in POPS], index=POPS)
            b = s.idxmax()
            result.append(str(b))

        f_out.write('\n'.join(result))


def concurrent_pipeline(files, window_len, n_threads):
    mid_files = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(process_file, file) for file in files]
        for future in as_completed(futures):
            mid_files.append(future.result())

        executor.map(process_probabilities, mid_files)


def main(files, n_threads, window_len):
    start_time = time.monotonic()
    if n_threads == 1:
        # file_in = files[0]
        # file_mid = process_file(file_in)
        file_mid = '../data/QuechuaCandelaria_3.local_an.txt'
        process_probabilities(file_mid, window_len=window_len)
        with open(f'../data/QuechuaCandelaria_3.local_an.txt_{window_len}.csv') as f:
            l = f.readlines()
            c = Counter([a.strip() for a in l])
            print(*[(a, c[a] / len(l)) for a in c], sep='\n')
            print(sum(c.values()) / len(l))
    else:
        concurrent_pipeline(files, window_len, n_threads)

    print(f'Finished in {time.monotonic() - start_time} sec.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Input filenames", nargs='+')
    parser.add_argument("-o", "--output", help="Output directory.", default='../data/')
    parser.add_argument("--n-threads", help="Number of threads to use.", type=int, default=1)
    parser.add_argument("--window-len", help="Window length to use.", type=int, default=250)

    args = parser.parse_args()
    # group = 'Zapotec_0',
    groups = ['QuechuaCandelaria_3', 'Atacama_3']
    filenames = [f'../data/STAT_America.{group}.txt_GENO.frq.strat.gz' for group in groups]
    main(filenames, args.n_threads, args.window_len)
    # if args.output:
    #     global OUTPUT_DIR
    #     OUTPUT_DIR = args.output
    #
    # if len(args.files) > 0:
    #     main(args.files, args.n_threads)
    # else:
    #     print('No filenames provided, please specify files to run on.')
