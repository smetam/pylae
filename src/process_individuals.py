# Usage:
# python3 src/process_individuals.py --distance-matrix configs/matrix.csv --mode fb --window-len 200 data/quechua_candelaria/QuechuaCandelaria_3.GA002786.txt

import time
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from collections import Counter
from scipy.special import softmax, xlogy


class ModelConfig:
    def __init__(self, args, populations):
        self.populations = sorted(populations)
        self.n_pops = len(self.populations)
        self.mode = args.mode
        self.win_len = args.window_len
        self.matrix = self._build_matrix(args.transition_matrix, args.distance_matrix)
        self._init_paths(args.file, args.output)

    def _init_paths(self, file, output):
        self.file_path = Path(file)
        self.filename = self.file_path.stem
        self.group, self.ind = self.filename.split('.')
        self.base_path = Path(output) if output else Path(file).parent
        self.base_path.mkdir(exist_ok=True, parents=True)

        self.input_file = str(self.file_path)
        self.snp_file = f'{self.base_path}/{self.group}_{self.mode}_{self.win_len}_snp_prob.tsv'
        self.prediction_file = f'{self.base_path}/{self.group}_{self.mode}_{self.win_len}_predictions.csv'
        self.stats_file = f'{self.base_path}/{self.group}_{self.mode}_{self.win_len}_stats.csv'

    @property
    def header(self):
        return (f"CHROM POS ID AF_{self.group} AF_" + " AF_".join(self.populations)).split()

    def _build_matrix(self, transition_matrix, distance_matrix):
        if transition_matrix:
            df = pd.read_csv(transition_matrix, index_col=0)
            sorted_header = list(sorted(list(df)))
            return df.sort_index()[sorted_header].values
        elif distance_matrix:
            df = pd.read_csv(distance_matrix, index_col=0)
            sorted_header = list(sorted(list(df)))
            df = df.sort_index()[sorted_header]
            return self._build_transition_matrix(df.values)
        else:
            return (np.ones((self.n_pops, self.n_pops)) + np.eye(self.n_pops)) / (self.n_pops + 1)

    def _build_transition_matrix(self, a):
        a_min = np.min(a)
        a_max = np.max(a)
        a = (a - a_min) / (a_max - a_min)
        a = np.exp(-3 * a)
        a = a / a.sum(axis=0, keepdims=1)
        return np.around(a, 3)


def convert_vcf(file):
    return file


def xlogx(x):
    return xlogy(x, x)


def emission_prob(snp: str, pop: int, row):
    """Probability of open state emission for population"""
    if snp == '1':
        p = row[pop] ** 2
    elif snp == '0':
         p = (1 - row[pop]) ** 2
    else:
        p = 2 * row[pop] * (1 - row[pop])
    return p + 0.0001


# def switch_prob(i: int, j: int, matrix):
#     """Probability of hidden state change from i to j"""
#     return matrix[i, j]


def fb_chunk(df, ind, outfile, n_pops, matrix):
    """Apply Forward-Backward algorithm to a chunk of data with win_len <= 500 nt"""
    n = len(df)
    a, b = np.zeros((n_pops, n)), np.zeros((n_pops, n))

    start_p = 1 / n_pops
    a_row = df.iloc[0]
    snp = a_row[ind]
    for pop_idx in range(n_pops):
        a[pop_idx, 0] = start_p * emission_prob(snp, pop_idx, a_row[4:].values)
        b[pop_idx, n - 1] = 1

    for i in range(1, n):
        a_row = df.iloc[i]
        a_snp = a_row[ind]
        a_prev_layer = a[:, i - 1]
        for pop_idx in range(n_pops):
            a[pop_idx, i] = sum(
                [a_prev_layer[j] * matrix[j, pop_idx] * emission_prob(a_snp, pop_idx, a_row[4:].values)
                 for j in range(n_pops)]
            )

        b_row = df.iloc[n - i]
        b_snp = b_row[ind]
        b_prev_layer = b[:, n - i]
        for pop_idx in range(n_pops):
            b[pop_idx, n - i - 1] = sum(
                [b_prev_layer[j] * matrix[j, pop_idx] * emission_prob(b_snp, j, b_row[4:].values)
                 for j in range(n_pops)]
            )

    p = np.zeros((n_pops, n))
    for i in range(n):
        for pop_idx in range(n_pops):
            p[pop_idx, i] = a[pop_idx, i] * b[pop_idx, i] / np.sum(a[:, n - 1])

    with open(outfile, 'a') as f_out:
        for i in range(n):
            snp_id = df.iloc[i, 2]
            line = p[:, i]
            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', line)) + '\n')


def fb_prob(file, outfile, group, populations, names, matrix, win_len=200):
    """Apply FB algorithm to input file by chunks"""
    ind = 'AF_' + group
    part = 0
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(populations) + '\n')

    for df in pd.read_csv(file, sep=' ', skiprows=1, names=names, chunksize=win_len, dtype={ind: object}):
        fb_chunk(df, ind, outfile, len(populations), matrix)
        part += 1

        if part % 20 == 0:
            print(f'Processed {part * win_len}.')


def softmax_prob(file, outfile, group, populations, names, scale=5):
    df = pd.read_csv(file, sep=' ', skiprows=1, names=names)
    df[f'AF_{group}'] = pd.to_numeric(df[f'AF_{group}'], errors='coerce').fillna(0.5)
    n_snp = len(df)
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(populations) + '\n')
        for i, row in df.iterrows():
            snp_id = row['ID']
            maf = row[f'AF_{group}']
            prob = row[4:].values.astype(float)
            p = softmax(scale * (1 - np.abs(prob - maf)))
            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', p)) + '\n')

            if i % 10000 == 0:
                print(f'Processed {i} / {n_snp}.')


def bayes_prob(file, outfile, group, populations, names):
    n_pops = len(populations)
    ind = 'AF_' + group
    df = pd.read_csv(file, sep=' ', skiprows=1, names=names)
    df[f'AF_{group}'] = pd.to_numeric(df[f'AF_{group}'], errors='coerce').fillna(0.5)
    n_snp = len(df)
    with open(outfile, 'w') as f_out:
        f_out.write('SNP_ID\t' + '\t'.join(populations) + '\n')
        for i, row in df.iterrows():
            p_snp = np.sum(row[4:].values) / n_pops + 0.0001
            snp_id = row['ID']
            snp = row[ind]
            if snp == 1.:
                p = (row[4:].values / n_pops / p_snp) ** 2
            elif snp == 0.:
                p = ((1 - row[4:].values) / n_pops / (1 - p_snp)) ** 2
            elif snp == 0.5:
                p = row[4:].values * (1 - row[4:].values) / (1 - p_snp) / p_snp / (n_pops ** 2)
            else:
                continue

            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.2f}', p)) + '\n')
            if i % 10000 == 0:
                print(f'Processed {i} / {n_snp}.')


def process_probabilities(file, output, populations, window_len=200):
    result = []

    with open(output, 'w') as f_out:
        for df in pd.read_csv(file, sep='\t', index_col=0, chunksize=window_len):
            s = pd.Series([-xlogx(df[pop].values).sum() for pop in populations], index=populations)
            b = s.idxmax()
            result.append(str(b))

        f_out.write('\n'.join(result))


def main(config):
    start_time = time.monotonic()

    if config.mode in ['fb', 'hmm']:

        fb_prob(config.input_file, config.snp_file, config.group,
                config.populations, config.header, config.matrix, config.win_len)
    elif config.mode == 'bayes':
        bayes_prob(config.input_file, config.snp_file, config.group,
                   config.populations, config.header)
    elif config.mode == 'softmax':
        softmax_prob(config.input_file, config.snp_file, config.group,
                     config.populations, config.header)
    else:
        print("Incorrect mode, please choose from 'bayes', 'fb', 'softmax'.")
        return

    print(f"Probabilities for each SNP are available at {config.snp_file}")

    process_probabilities(config.snp_file, config.prediction_file,
                          config.populations, window_len=config.win_len)
    print(f"Predictions for each window are available at {config.prediction_file}")

    with open(config.prediction_file) as f:
        lines = f.readlines()
    counts = Counter([line.strip() for line in lines])
    res = [(pop, counts[pop] / len(lines)) for pop in config.populations]
    series = pd.Series([x[1] for x in res], index=[x[0] for x in res])
    series.to_csv(config.stats_file, header=False)
    print(f"Overall stats  are available at {config.stats_file}")

    print(f'Finished in {time.monotonic() - start_time} sec.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input filename")
    parser.add_argument("-o", "--output", default=None,
                        help="Output directory. Will be created automatically. "
                             "If already exists, some files may be modified")
    # parser.add_argument("--n-threads", help="Number of threads to use.", type=int, default=1)
    # parser.add_argument("--convert_vcf", help="Use vcf as an input file. It will be transformed.",
    # action='store_true', default=False)
    parser.add_argument("--window-len", help="Window length to use.", type=int, default=250)
    parser.add_argument("--transition-matrix", help="Path to transition matrix for HMM in csv format.")
    parser.add_argument("--distance-matrix", help="Path to distance matrix for HMM in csv format. "
                                                  "If transition matrix specified this option is ignored.")
    parser.add_argument("--mode", choices=['bayes', 'fb', 'softmax'], default='fb',
                        help="Calculation mode, should be one of 'bayes', 'fb', 'softmax'.")

    args = parser.parse_args()

    populations = [
        'Mediterranean', 'NativeAmerican', 'NorthEastAsian', 'NorthernEuropean',
        'Oceanian', 'SouthAfrican', 'SouthEastAsian', 'SouthWestAsian', 'SubsaharanAfrican'
    ]
    model_config = ModelConfig(args, populations)
    main(model_config)
