# Usage:
# python3 src/process_individuals.py --distance-matrix configs/matrix.csv --mode fb --window-len 200 data/quechua_candelaria/QuechuaCandelaria_3.GA002786.txt

import time
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from collections import Counter, namedtuple


Record = namedtuple('Record', ('chrom', 'start', 'end', 'confidence', 'prediction'))


class ModelConfig:
    def __init__(self, args, populations):
        self.populations = sorted(populations)
        self.n_pops = len(self.populations)
        self.mode = 'bayes'
        self.window_len = args.window_len
        self.admixtures_file = args.admixtures
        self._init_paths(args.file, args.output)
        self.start_time = time.monotonic()

    def _init_paths(self, file, output):
        self.file_path = Path(file)
        self.filename = self.file_path.stem
        self.group, sample = self.filename.split('.')
        self.sample = args.sample or sample
        self.base_path = Path(output) if output else Path(file).parent
        self.base_path.mkdir(exist_ok=True, parents=True)

        self.input_file = str(self.file_path)
        self.snp_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_snp_prob.tsv'
        self.prediction_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_predictions.csv'
        self.results_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_result.csv'
        self.stats_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_stats.csv'

    @property
    def header(self):
        return (f"CHROM POS AF_{self.sample} AF_" + " AF_".join(self.populations)).split()


def run_bayes(config, alpha=0.0001):
    n_pops = config.n_pops
    sample_frequency = 'AF_' + config.sample
    df = pd.read_csv(config.input_file, sep=' ', skiprows=1, names=config.header)
    df[sample_frequency] = pd.to_numeric(df[sample_frequency], errors='coerce').fillna(0)
    print(df.info())
    n_snp = len(df)
    with open(config.snp_file, 'w') as f_out:
        f_out.write('CHROM\tPOS\t' + '\t'.join(config.populations) + '\n')
        for i, row in df.iterrows():
            snp_id = f"{int(row['CHROM'])}\t{int(row['POS'])}"
            snp = row[sample_frequency]
            if snp > 0.99:
                p = (row[3:].values ** 2 + alpha) / (1 + n_pops * alpha)
            elif snp < 0.01:
                p = ((1 - row[3:].values) ** 2 + alpha) / (1 + n_pops * alpha)
            else:
                p = (2 * row[3:].values * (1 - row[3:].values) + alpha) / (1 + n_pops * alpha)

            f_out.write(f'{snp_id}\t' + '\t'.join(map(lambda x: f'{x:.6f}', p / np.sum(p))) + '\n')

            if i % 100000 == 0:
                print(f'Processed {i} / {n_snp}.')

    print(f"Probabilities for each SNP are available at {config.snp_file}")


def solve_region_smoothed(df, populations, pop_prob, prev, chrom, alpha=100):
    start = df.iloc[0, 1]
    end = df.iloc[-1, 1]
    s = pd.Series([np.log(df[pop].values).sum()
                   + np.log(pop_prob[pop] / (1 + alpha) if pop != prev else (pop_prob[pop] + alpha) / (1 + alpha))
                   for pop in populations], index=populations)
    first, second = sorted(s.nlargest(2))
    best = s.idxmax()
    return f'{chrom},{start},{end},{first/second:.5f},{best}', best


def solve_region(df, populations, pop_prob, prev, chrom):
    start = df.iloc[0, 1]
    end = df.iloc[-1, 1]
    s = pd.Series([np.log(df[pop].values).sum() + np.log(pop_prob[pop])
                   for pop in populations], index=populations)
    first, second = sorted(s.nlargest(2))
    best = s.idxmax()
    return f'{chrom},{start},{end},{first/second:.5f},{best}', best


def process_probabilities(config):
    result = []
    q = 0
    prev = None
    df = pd.read_csv(config.admixtures_file, sep=',', index_col=0)
    pop_prob = df.loc[config.sample, :]

    for df in pd.read_csv(config.snp_file, sep='\t', chunksize=config.window_len, dtype={"CHROM": int, "POS": int}):
        chrom_start = df.iloc[0, 0]
        chrom_end = df.iloc[-1, 0]
        if chrom_start == chrom_end:
            rec, prev = solve_region(df, config.populations, pop_prob, prev, chrom_start)
            result.append(rec)
        else:
            for chrom in (chrom_start, chrom_end):
                tdf = df[df['CHROM'] == chrom]
                rec, prev = solve_region(tdf, config.populations, pop_prob, prev, chrom)
                result.append(rec)

        q += 1
        if q % 1000 == 0:
            print(q, 'windows processed.')

    with open(config.prediction_file, 'w') as f_out:
        f_out.write('\n'.join(result))

    print(pop_prob)


def calculate_stats(config):
    with open(config.prediction_file) as f:
        lines = f.readlines()
    counts = Counter([line.strip().split(',')[4] for line in lines])
    res = [(pop, counts[pop] / len(lines)) for pop in config.populations]

    series = pd.Series([x[1] for x in res], index=[x[0] for x in res])
    series.to_csv(config.stats_file, header=False)
    print(f"Overall stats are available at {config.stats_file}")


def merge_windows(config):
    res = []
    with open(config.prediction_file, 'r') as f_in:
        line = f_in.readline()
        prev_record = Record(*line.strip().split(','))
        for line in f_in.readlines():
            record = Record(*line.strip().split(','))
            if prev_record.chrom == record.chrom \
                    and prev_record.prediction == record.prediction:
                conf = float(prev_record.confidence) * float(record.confidence)
                prev_record = Record(record.chrom, prev_record.start, record.end, conf, record.prediction)
            else:
                res.append(f'{prev_record.chrom},{prev_record.start},{prev_record.end},'
                           f'{prev_record.confidence},{prev_record.prediction}')
                prev_record = record
    with open(config.results_file, 'w') as f_out:
        f_out.write('\n'.join(res))


def main(config):

    run_bayes(config)

    process_probabilities(config)

    calculate_stats(config)

    merge_windows(config)

    print(f'Finished in {time.monotonic() - config.start_time} sec.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input filename")
    parser.add_argument("-o", "--output", default=None,
                        help="Output directory. Will be created automatically. "
                             "If already exists, some files may be modified")
    parser.add_argument("--sample", type=str, default=None,
                        help="Sample name. If not specified will be inferred form input filename")
    parser.add_argument("--window-len", help="Window length to use.", type=int, default=250)
    parser.add_argument("--admixtures", help="Csv file with admixture vectors for sample.",
                        type=str, default='configs/admixtures.csv')

    args = parser.parse_args()

    # populations = [
    #     'Mediterranean', 'NativeAmerican', 'NorthEastAsian', 'NorthernEuropean',
    #     'Oceanian', 'SouthAfrican', 'SouthEastAsian', 'SouthWestAsian', 'SubsaharanAfrican'
    # ]

    populations = [
        'Amazonian', 'Andamanese', 'Austronesian', 'BrazilianYanomami', 'Dravidian', 'EastAsian', 'EastIndian',
        'Eskimo', 'Malaysian', 'NearEastern', 'NorthernEuropean', 'Papuan', 'PapuanBaining', 'PlillippinoNegrito',
        'SouthAmerican_Chaco', 'SouthAmerican_Quechua', 'SubSaharanAfrican', 'Yeniseyan'
    ]

    model_config = ModelConfig(args, populations)
    main(model_config)
