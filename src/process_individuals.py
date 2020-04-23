# Usage:
# python3 process_individuals.py --mode fb --window-len 200  <filename> ...

import time
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from copy import deepcopy
from collections import Counter
from scipy.special import softmax
from concurrent.futures import ThreadPoolExecutor, as_completed


POPS = sorted([
    'Mediterranean', 'NativeAmerican', 'NorthEastAsian', 'NorthernEuropean',
    'Oceanian', 'SouthAfrican', 'SouthEastAsian', 'SouthWestAsian', 'SubsaharanAfrican'
])
N_POPS = len(POPS)
state_matrix = np.array(
    [[0.2  , 0.09 , 0.094, 0.164, 0.084, 0.083, 0.089, 0.114, 0.083],
     [0.11 , 0.171, 0.126, 0.11 , 0.095, 0.089, 0.101, 0.108, 0.089],
     [0.112, 0.123, 0.166, 0.112, 0.096, 0.09 , 0.102, 0.11 , 0.09 ],
     [0.164, 0.09 , 0.094, 0.2  , 0.084, 0.083, 0.089, 0.114, 0.083],
     [0.112, 0.104, 0.107, 0.112, 0.161, 0.092, 0.108, 0.11 , 0.092],
     [0.109, 0.096, 0.099, 0.109, 0.091, 0.172, 0.095, 0.107, 0.123],
     [0.114, 0.105, 0.109, 0.114, 0.102, 0.092, 0.16 , 0.112, 0.092],
     [0.128, 0.099, 0.103, 0.128, 0.092, 0.091, 0.099, 0.169, 0.091],
     [0.109, 0.096, 0.099, 0.109, 0.091, 0.123, 0.095, 0.107, 0.172]]
)


def emission_prob(snp: str, pop: int, row):
    """Probability of open state emission for population"""
    if snp == '1':
        p = row[pop] ** 2
    elif snp == '0':
         p = (1 - row[pop]) ** 2
    else:
        p = 2 * row[pop] * (1 - row[pop])
    return p + 0.0001


def switch_prob(i: int, j: int):
    """Probability of hidden state change from i to j"""
    return 2 / (N_POPS + 1) if i == j else 1 / (N_POPS + 1)


def fb_chunk(df, ind, outfile):
    """Apply Forward-Backward algorithm to a chunk of data with win_len <= 500 nt"""
    n = len(df)
    a, b = np.zeros((N_POPS, n)), np.zeros((N_POPS, n))

    start_p = 1 / N_POPS
    a_row = df.iloc[0]
    snp = a_row[ind]
    for pop_idx in range(N_POPS):
        a[pop_idx, 0] = start_p * emission_prob(snp, pop_idx, a_row[4:].values)
        b[pop_idx, n - 1] = 1

    for i in range(1, n):
        a_row = df.iloc[i]
        a_snp = a_row[ind]
        a_prev_layer = a[:, i - 1]
        for pop_idx in range(N_POPS):
            a[pop_idx, i] = sum(
                [a_prev_layer[j] * switch_prob(j, pop_idx) * emission_prob(a_snp, pop_idx, a_row[4:].values)
                 for j in range(N_POPS)]
            )

        b_row = df.iloc[n - i]
        b_snp = b_row[ind]
        b_prev_layer = b[:, n - i]
        for pop_idx in range(N_POPS):
            b[pop_idx, n - i - 1] = sum(
                [b_prev_layer[j] * switch_prob(pop_idx, j) * emission_prob(b_snp, j, b_row[4:].values)
                 for j in range(N_POPS)]
            )

    p = np.zeros((N_POPS, n))
    for i in range(n):
        for pop_idx in range(N_POPS):
            p[pop_idx, i] = a[pop_idx, i] * b[pop_idx, i] / np.sum(a[:, n - 1])

    with open(outfile, 'a') as f_out:
        for i in range(n):
            snp_id = df.iloc[i, 2]
            line = p[:, i]
            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', line)) + '\n')


def fb_prob(file, outfile, group, win_len=200):
    """Apply FB algorithm to input file by chunks"""
    ind = 'AF_' + group
    part = 0
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(POPS) + '\n')

    for df in pd.read_csv(file, sep=' ', chunksize=win_len, dtype={ind: object}):
        fb_chunk(df, ind, outfile)
        part += 1

        if part % 20 == 0:
            print(f'Processed {part * win_len}.')


def softmax_prob(file, outfile, group, scale=5):
    df = pd.read_csv(file, sep=' ')
    df[f'AF_{group}'] = pd.to_numeric(df[f'AF_{group}'], errors='coerce').fillna(0.5)
    n_snp = len(df)
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(POPS) + '\n')
        for i, row in df.iterrows():
            snp_id = row['ID']
            maf = row[f'AF_{group}']
            prob = row[4:].values.astype(float)
            p = softmax(scale * (1 - np.abs(prob - maf)))
            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', p)) + '\n')

            if i % 10000 == 0:
                print(f'Processed {i} / {n_snp}.')


def bayes_prob(file, outfile, group):
    ind = 'AF_' + group
    df = pd.read_csv(file, sep=' ')
    df[f'AF_{group}'] = pd.to_numeric(df[f'AF_{group}'], errors='coerce').fillna(0.5)
    n_snp = len(df)
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(POPS) + '\n')
        for i, row in df.iterrows():
            p_snp = np.sum(row[4:].values) / N_POPS + 0.0001
            snp_id = row['ID']
            snp = row[ind]
            if snp == 1.:
                p = (row[4:].values / N_POPS / p_snp) ** 2
            elif snp == 0.:
                p = ((1 - row[4:].values) / N_POPS / (1 - p_snp)) ** 2
            elif snp == 0.5:
                p = row[4:].values * (1 - row[4:].values) / (1 - p_snp) / p_snp / (N_POPS ** 2)
            else:
                continue

            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', p)) + '\n')
            if i % 10000 == 0:
                print(f'Processed {i} / {n_snp}.')


def process_probabilities(file, output, window_len=300):
    result = []

    with open(output, 'w') as f_out:
        for df in pd.read_csv(file, sep='\t', index_col=0, chunksize=window_len):
            s = pd.Series([df[pop].sum() for pop in POPS], index=POPS)
            b = s.idxmax()
            result.append(str(b))

        f_out.write('\n'.join(result))


def main(file, win_len, mode, n_threads=4):
    start_time = time.monotonic()

    file_path = Path(file)
    filename = file_path.stem
    group, ind = filename.split('.')
    base_path = Path(file).parent

    input_file = str(file_path)
    snp_prob_file = f'{base_path}/{group}_{mode}_{win_len}_snp_prob.txt'
    prediction_file = f'{base_path}/{group}_{mode}_{win_len}_predictions.csv'
    stats_file = f'{base_path}/{group}_{mode}_{win_len}_stats.csv'

    if mode in ['fb', 'hmm']:
        fb_prob(input_file, snp_prob_file, group)
    elif mode == 'bayes':
        bayes_prob(input_file, snp_prob_file, group)
    elif mode == 'softmax':
        softmax_prob(input_file, snp_prob_file, group)
    else:
        print("Incorrect mode, please choose from 'bayes', 'fb', 'softmax'.")
        return

    print(f"Probabilities for each SNP are available at {snp_prob_file}")

    process_probabilities(snp_prob_file, prediction_file, window_len=win_len)
    print(f"Predictions for each window are available at {prediction_file}")

    with open(prediction_file) as f:
        lines = f.readlines()
    counts = Counter([line.strip() for line in lines])
    res = [(pop, counts[pop] / len(lines)) for pop in counts]
    series = pd.Series([x[1] for x in res], index=[x[0] for x in res])
    series.to_csv(stats_file)
    print(f"Overall stats  are available at {stats_file}")

    print(f'Finished in {time.monotonic() - start_time} sec.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input filename")
    # parser.add_argument("-o", "--output", help="Output directory.", default='../data/')
    # parser.add_argument("--n-threads", help="Number of threads to use.", type=int, default=1)
    parser.add_argument("--window-len", help="Window length to use.", type=int, default=250)
    parser.add_argument("--mode", choices=['bayes', 'fb', 'softmax'], default='fb',
                        help="Calculation mode, should be one of 'bayes', 'fb', 'softmax'.")

    args = parser.parse_args()

    main(args.file, args.window_len, args.mode.lower())
