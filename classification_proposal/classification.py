import argparse
import glob
import json
import logging
import os
import subprocess
from subprocess import PIPE
import sys
from ACF import analysis_complete_flag as acf


class ClassificationProposal:
    def __init__(self, intsv, outdir, indir):
        self.intsv = intsv
        self.outdir = outdir
        self.indir = indir
        os.makedirs(self.outdir, exist_ok=True)
        self.log = os.path.join(outdir, 'log.txt')
        self.setup_logger()

    def setup_logger(self, name=None):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(self.log)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s '
            '%(filename)s %(funcName)s :\n%(message)s')
        fh.setFormatter(fh_formatter)
        # creates a file handler that logs messages above INFO level
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        sh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s : %(message)s',
            '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        self.logger = logger

    def execute_cmd(self, cmd, shell=True):
        self.logger.info('[cmd] {}'.format(cmd))
        proc = subprocess.run(
            cmd, shell=shell, stdout=PIPE,
            stderr=PIPE, text=True)
        stdout = proc.stdout
        stderr = proc.stderr
        if len(stdout) > 0:
            self.logger.info(stdout)
        if len(stderr) > 0:
            self.logger.info(stderr)
        return stdout, stderr

    def logging_title(self, title):
        title = f'#### {title.strip()} ####'
        sharps = '#' * len(title)
        msg = f'\n{sharps}\n{title}\n{sharps}'
        self.logger.info(msg)

    def classify_method1(self):
        self.count_dict = {}
        query_contig = ''
        self.ref_set = set()
        test_len = set()
        with open(self.intsv, 'r') as intsv:
            lines = intsv.readlines()
        # load the most similar contig(s)
        for line in lines:
            data = [i.strip() for i in line.split(sep='\t')]
            ref = '__'.join(data[1].split(sep='__')[:-1])
            contig = data[0]
            self.ref_set.add(ref)
            test_len.add(data[0])
            if len(data) < 12:
                continue
            if query_contig != contig:
                max_score = float(data[-1])
                self.count_dict[contig] = {'count': [ref]}
                query_contig = contig
                continue
            if max_score == float(data[-1]):
                self.count_dict[contig]['count'].append(ref)
        # scoring per ref
        print(len(test_len))
        self.ref_score_dict = {
            i: {'count': 0, 'ratio': 0} for i in self.ref_set}
        query_contigs = 0
        for contig in self.count_dict.keys():
            all_contigs = self.count_dict[contig]['count']
            query_contigs += 1
            for ref in self.ref_set:
                frac = all_contigs.count(ref) / len(all_contigs)
                self.ref_score_dict[ref]['count'] += frac
        for ref in self.ref_set:
            count = self.ref_score_dict[ref]['count']
            self.ref_score_dict[ref]['ratio'] = count / query_contigs
        # sample name
        sample_name = os.path.basename(self.intsv)[:-4]
        # output json
        outj = os.path.join(self.outdir, f'{sample_name}.json')
        with open(outj, 'w') as outjp:
            json.dump(
                {
                    'query_contig_num': query_contigs,
                    'data': self.ref_score_dict
                    },
                outjp,
                indent=4)
        # output summary tsv
        outobj = [
            [
                i, self.ref_score_dict[i]['count'],
                self.ref_score_dict[i]['ratio']
                ]
            for i in self.ref_score_dict.keys()]
        outobj.sort(key=lambda x: x[2], reverse=True)
        outobj = ['#sample\tcount\tratio'] + \
            [f'{i[0]}\t{i[1]}\t{i[2]}' for i in outobj]
        outt = os.path.join(self.outdir, f'{sample_name}_summary.tsv')
        with open(outt, 'w') as outtp:
            outtp.write('\n'.join(outobj))

    def run(self):
        if self.indir is None or self.indir == 'None':
            self.classify_method1()
        else:
            indir_star = os.path.join(self.indir, '*')
            tsvs = set(glob.glob(indir_star))
            for tsv in tsvs:
                self.intsv = tsv
                self.classify_method1()


def argument_parser():
    parser = argparse.ArgumentParser(description="blast res classifier")
    parser.add_argument(
        '-i', '--intsv', dest='intsv', type=str,
        help='blast tsv', required=False)
    parser.add_argument(
        '-o', '--output_dir', dest='outdir', type=str,
        help='output directory', required=True)
    parser.add_argument(
        '-indir', '--input_dir', dest='indir', type=str,
        help='input tsv directory', required=False,
        default=None)
    return parser.parse_args()


def main():
    args = argument_parser()
    CP = ClassificationProposal(
        intsv=args.intsv,
        outdir=args.outdir,
        indir=args.indir)
    CP.run()
    # analysis complete flag
    acf(outdir=args.outdir)


if __name__ == '__main__':
    main()
