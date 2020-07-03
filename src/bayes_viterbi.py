# Usage:
# python3 src/bayes_viterbi.py --sample GA001371 --window-len 500 data/cdx0.GA001371.txt -o data/res --admixtures configs/admixtures.csv


import time
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcol

from matplotlib import path
from copy import copy
from pathlib import Path
from collections import Counter, namedtuple


Record = namedtuple('Record', ('chrom', 'start', 'end', 'confidence', 'prediction'))


class ModelConfig:
    def __init__(self, args, populations):
        self.populations = populations
        self.n_pops = len(self.populations)
        self.mode = 'bayes_viterbi'
        self.plot = args.plot
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
        self.snp_file = '{}/{}_bayes_snp_prob.tsv'.format(self.base_path, self.sample)
        self.hmm_input_file = '{}/{}_{}_{}_viterbi_input.csv'.format(self.base_path, self.sample, self.mode, self.window_len)
        self.prediction_file = '{}/{}_{}_{}_predictions.csv'.format(self.base_path, self.sample, self.mode, self.window_len)
        self.results_file = '{}/{}_{}_{}_result.csv'.format(self.base_path, self.sample, self.mode, self.window_len)
        self.stats_file = '{}/{}_{}_{}_stats.csv'.format(self.base_path, self.sample, self.mode, self.window_len)
        self.plot_file = '{}/{}_{}_{}_plot.png'.format(self.base_path, self.sample, self.mode, self.window_len)

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
        return ("CHROM POS AF_{} AF_".format(self.sample) + " AF_".join(self.populations)).split()


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
            snp_id = "{}\t{}".format(int(row['CHROM']), int(row['POS']))
            snp = row[sample_frequency]
            if snp > 0.99:
                p = (row[3:].values ** 2 + alpha) / (1 + n_pops * alpha)
            elif snp < 0.01:
                p = ((1 - row[3:].values) ** 2 + alpha) / (1 + n_pops * alpha)
            else:
                p = (2 * row[3:].values * (1 - row[3:].values) + alpha) / (1 + n_pops * alpha)

            f_out.write(str(snp_id) + '\t' + '\t'.join(map(lambda x: str(round(x, 6)), p / np.sum(p))) + '\n')

            if i % 100000 == 0:
                print('Processed {} / {}.'.format(i, n_snp))

    print("Probabilities for each SNP are available at {}".format(config.snp_file))


def solve_region_smoothed(df, populations, pop_prob, chrom, alpha=100):
    start = df.iloc[0, 1]
    end = df.iloc[-1, 1]
    prob = [str(np.round(np.log(df[pop].values).sum() + np.log(pop_prob[pop]), 2)) for pop in populations]
    return '{},{},{},'.format(chrom, start, end) + ','.join(prob)


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
    print('Total error: ', np.sum((df['Predicted'] - df['Prior']) ** 2) / len(df['Predicted']))
    print("Overall stats are available at ", config.stats_file)


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
                res.append('{},{},{},{},{}'.format(prev_record.chrom, prev_record.start, prev_record.end,
                                                   prev_record.confidence, prev_record.prediction))
                prev_record = record
    res.append('{},{},{},{},{}'.format(prev_record.chrom, prev_record.start, prev_record.end,
                                       prev_record.confidence, prev_record.prediction))

    with open(config.results_file, 'w') as f_out:
        f_out.write('\n'.join(res))


def plot_rects(anc, chrom, start, stop, pop_order, colors, ax, chrX=False, conf=1):
    conf *= 0.7
    verts = [
            (float(start), chrom),  # left, bottom
            (float(start), chrom + conf),  # left, top
            (float(stop), chrom + conf),  # right, top
            (float(stop), chrom),  # right, bottom
            (0, 0),  # ignored
        ]

    codes = [
        path.Path.MOVETO,
        path.Path.LINETO,
        path.Path.LINETO,
        path.Path.LINETO,
        path.Path.CLOSEPOLY,
    ]

    clip_path = path.Path(verts, codes)
    if anc in pop_order:
        col = mcol.PathCollection([clip_path], facecolor=colors[pop_order.index(anc)], linewidths=0)
    else:
        col = mcol.PathCollection([clip_path], facecolor=colors[-1], linewidths=0)
    return col


def plot_admixture(admix_file, stats_file=None, save_name=None):
    df = pd.read_csv(admix_file, header=None,
                     names=['chr', 'start', 'end', 'conf', 'pop'])
    if stats_file:
        stat_df = pd.read_csv(stats_file, sep='\t', index_col=0)
        largest_pops = list(stat_df['Predicted'].nlargest(5).index)
        largest_pops.append('Other')
    else:
        largest_pops = list(df['pop'].unique())

    colors = ['#4ECBF5', '#368CA8', '#F5AF4E', '#A86636', '#383838', '#C7C7C7',
              '#F2EB61']
    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111)
    ax.set_xlim(-5, 260)
    ax.set_ylim(23, 0)
    plt.yticks(range(1, 24))
    plt.xlabel('Genomic position (Mbp)')
    plt.ylabel('Chromosome')
    plt.title('Local ancestry')

    p = []
    for i in range(len(largest_pops)):
        p.append(plt.Rectangle((0, 0), 1, 1, color=colors[i]))
    p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
    labs = list(largest_pops)

    leg = ax.legend(p, labs, loc=4, fancybox=True)
    leg.get_frame().set_alpha(0)

    for i, row in df.iterrows():
        popul = row['pop']
        if popul not in largest_pops:
            popul = 'Other'
        col = plot_rects(popul, row['chr'], row['start'] / 10 ** 6,
                         row['end'] / 10 ** 6, largest_pops, colors, ax)
        ax.add_collection(col)
    if save_name:
        fig.savefig(save_name)
    print("Admixture plot available at ", save_name)
    # plt.show()


def main(config):

    if not Path(config.snp_file).exists():
        run_bayes(config)

    if not Path(config.hmm_input_file).exists():
        process_probabilities(config)

    run_viterbi(config)

    calculate_stats(config)

    merge_windows(config)

    if config.plot:
        plot_admixture(config.results_file, stats_file=config.stats_file,
                       save_name=config.plot_file)

    print('Finished in {} sec.'.format(time.monotonic() - config.start_time))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input filename")
    parser.add_argument("-o", "--output", default=None,
                        help="Output directory. Will be created automatically. "
                             "If already exists, some files may be modified")
    parser.add_argument("--sample", type=str, default=None,
                        help="Sample name. If not specified will be inferred form input filename")
    parser.add_argument("--window-len", help="Window length to use.", type=int, default=250)
    parser.add_argument("--plot", help="Visualize results", default=False, action='store_true')
    parser.add_argument("--admixtures", help="Csv file with admixture vectors for sample.",
                        type=str, default=None)

    args = parser.parse_args()

    with open(args.file, 'r') as f_in:
        header = f_in.readline()
    s = header.split()[4:]
    pop_list = [x.split('AF_')[1] for x in s]
    print('Population list is inferred from header, using following populations: ')
    print(", ".join(pop_list))

    model_config = ModelConfig(args, pop_list)
    main(model_config)
