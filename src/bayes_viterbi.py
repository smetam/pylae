# Usage:
# python3 src/bayes_hmm.py --sample GA001371 --window-len 500 data/cdx0.GA001371.txt -o data/res --admixtures configs/admixtures.csv


import time
import argparse
import pandas as pd
import numpy as np

from copy import copy
from pathlib import Path
from collections import Counter, namedtuple


Record = namedtuple('Record', ('chrom', 'start', 'end', 'confidence', 'prediction'))


class ModelConfig:
    def __init__(self, args, populations):
        self.populations = sorted(populations)
        self.n_pops = len(self.populations)
        self.mode = 'bayes_hmm'
        self.window_len = args.window_len
        self._init_paths(args.file, args.output, args.sample)
        self._set_admixtures(args.admixtures)
        self.start_time = time.monotonic()

    def _init_paths(self, file, output, sample):
        self.file_path = Path(file)
        self.filename = self.file_path.stem
        if sample:
            self.sample = sample
        else:
            self.group, self.sample = self.filename.split('.')

        self.base_path = Path(output) if output else Path(file).parent
        self.base_path.mkdir(exist_ok=True, parents=True)

        self.input_file = str(self.file_path)
        self.snp_file = f'{self.base_path}/{self.sample}_bayes_snp_prob.tsv'
        self.hmm_input_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_viterbi_input.csv'
        self.prediction_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_predictions.csv'
        self.results_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_result.csv'
        self.stats_file = f'{self.base_path}/{self.sample}_{self.mode}_{self.window_len}_stats.csv'

    def _set_admixtures(self, admixtures_file):
        if admixtures_file:
            self.admixtures_file = admixtures_file
            df = pd.read_csv(self.admixtures_file, sep=',', index_col=0)
            self.admixtures = df.loc[self.sample, :]
        else:
            self.admixtures_file = None
            self.admixtures = pd.Series([1/self.n_pops] * self.n_pops, index=self.populations)

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


def solve_region_smoothed(df, populations, pop_prob, chrom, alpha=100):
    start = df.iloc[0, 1]
    end = df.iloc[-1, 1]
    return f'{chrom},{start},{end},' + ','.join([f'{np.log(df[pop].values).sum() + np.log(pop_prob[pop]):.2f}'
                                                for pop in populations])


def process_probabilities(config):
    result = []
    q = 0
    pop_prob = config.admixtures

    for df in pd.read_csv(config.snp_file, sep='\t', chunksize=config.window_len, dtype={"CHROM": int, "POS": int}):
        chrom_start = df.iloc[0, 0]
        chrom_end = df.iloc[-1, 0]
        if chrom_start == chrom_end:
            rec = solve_region_smoothed(df, config.populations, pop_prob, chrom_start)
            result.append(rec)
        else:
            for chrom in (chrom_start, chrom_end):
                tdf = df[df['CHROM'] == chrom]
                rec = solve_region_smoothed(tdf, config.populations, pop_prob, chrom)
                result.append(rec)

        q += 1
        if q % 1000 == 0:
            print(q, 'windows processed.')

    with open(config.hmm_input_file, 'w') as f_out:
        f_out.write('\n'.join(result))

    print(pop_prob)


def run_viterbi(config, alpha=1):
    df = pd.read_csv(config.hmm_input_file, names=['chrom', 'start', 'end'] + config.populations)
    n = len(df)
    d = []
    line = np.zeros(config.n_pops)
    new_line = copy(line)
    last = np.arange(config.n_pops)

    for i, row in df.iterrows():

        for j in range(config.n_pops):
            transition_penalty = np.log(np.ones(config.n_pops) / (config.n_pops + alpha))
            transition_penalty[j] = np.log((1 + alpha) / (config.n_pops + alpha))
            emission = row.values[j + 3]

            p = line + transition_penalty + emission
            # p = line + emission
            index = np.argmax(p)
            new_line[j] = p[index]
            last[j] = index

        d.append(copy(last))
        line = copy(new_line)

    index = np.argmax(line)
    trail = [config.populations[index]]
    for i in range(n - 1):
        index = d.pop()[index]
        trail.append(config.populations[index])
    df['trail'] = list(reversed(trail))
    df['conf'] = 1

    df[['chrom', 'start', 'end', 'conf', 'trail']].to_csv(config.prediction_file, header=False, index=False)


def calculate_stats(config):
    with open(config.prediction_file) as f:
        lines = f.readlines()
    counts = Counter([line.strip().split(',')[4] for line in lines])
    res = [(pop, counts[pop] / len(lines)) for pop in config.populations]

    series = pd.Series([x[1] for x in res], index=[x[0] for x in res])
    # series.to_csv(config.stats_file, header=False)
    df = pd.DataFrame({'Predicted': series, 'Prior': config.admixtures})
    df.to_csv(config.stats_file, sep='\t')
    print('Total error: ', np.sum((df['Predicted'] - df['Prior']) ** 2))
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
                conf = round(float(prev_record.confidence) * float(record.confidence), 5)
                prev_record = Record(record.chrom, prev_record.start, record.end, conf, record.prediction)
            else:
                res.append(f'{prev_record.chrom},{prev_record.start},{prev_record.end},'
                           f'{prev_record.confidence},{prev_record.prediction}')
                prev_record = record
    res.append(f'{prev_record.chrom},{prev_record.start},{prev_record.end},'
               f'{prev_record.confidence},{prev_record.prediction}')

    with open(config.results_file, 'w') as f_out:
        f_out.write('\n'.join(res))


def main(config):

    if not Path(config.snp_file).exists():
        run_bayes(config)

    if not Path(config.hmm_input_file).exists():
        process_probabilities(config)

    run_viterbi(config)

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
                        type=str, default=None)

    args = parser.parse_args()

    # populations = [
    #     'Mediterranean', 'NativeAmerican', 'NorthEastAsian', 'NorthernEuropean',
    #     'Oceanian', 'SouthAfrican', 'SouthEastAsian', 'SouthWestAsian', 'SubsaharanAfrican'
    # ]

    population_list = [
        'Amazonian', 'Andamanese', 'Austronesian', 'BrazilianYanomami', 'Dravidian', 'EastAsian', 'EastIndian',
        'Eskimo', 'Malaysian', 'NearEastern', 'NorthernEuropean', 'Papuan', 'PapuanBaining', 'PlillippinoNegrito',
        'SouthAmerican_Chaco', 'SouthAmerican_Quechua', 'SubSaharanAfrican', 'Yeniseyan'
    ]

    model_config = ModelConfig(args, population_list)
    main(model_config)
