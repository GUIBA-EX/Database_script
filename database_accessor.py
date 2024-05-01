import argparse
import configparser
import glob
import gzip
import logging
import os
import subprocess
from subprocess import PIPE
import sys
from ACF import analysis_complete_flag as acf


class DatabaseAccessor:
    def __init__(self, conf):
        self.conf = conf
        self.parse_config()
        self.indir = self.settings['input_dir']
        self.outdir = self.settings['output_dir']
        os.makedirs(self.outdir, exist_ok=True)
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

    def parse_config(self):
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(self.conf, 'UTF-8')
        self.tools = dict(conf.items('tools'))
        self.settings = dict(conf.items('settings'))
        self.log = os.path.join(self.settings['output_dir'], 'log.txt')

    def logging_title(self, title):
        title = f'#### {title.strip()} ####'
        sharps = '#' * len(title)
        msg = f'\n{sharps}\n{title}\n{sharps}'
        self.logger.info(msg)

    def copy_conf(self):
        copy_cmd = 'cp {conf} {outdir}'
        cmd = copy_cmd.format(
            outdir=self.settings['output_dir'],
            conf=self.conf)
        self.execute_cmd(cmd)

    def inputdir_to_ini(self):
        # indir
        if self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip loading fastq')
            return 0
        self.logging_title('load fastq dir')
        indir_star = os.path.join(self.indir, '*')
        fastqs = set(glob.glob(indir_star))
        self.files_to_process = {}
        outlines = ['#samplename,forward,reverse']
        for fq in fastqs:
            if fq.endswith('_R1_001.fastq.gz'):
                basename = os.path.basename(fq)
                paired = fq[:-16] + '_R2_001.fastq.gz'
                if not os.path.isfile(paired):
                    msg = f'no such a paired-end file : {paired}'
                    self.logger.info(msg)
                    continue
                sample = basename[:-21]
                basename = os.path.basename(fq)
                paired_basename = os.path.basename(paired)
                outline = f'{sample},{basename},{paired_basename}'
                outlines.append(outline)
                self.files_to_process[sample] = [fq, paired]
            else:
                continue
        with open(os.path.join(self.outdir, 'sample.csv'), 'w') as ini:
            ini.write('\n'.join(outlines))
        self.samples = self.files_to_process.keys()

    def quality_control(self):
        # quality cutoff
        if self.settings['quality_control'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip quality control')
            return 0
        self.logging_title('quality control : fastq')
        fqfd = os.path.join(self.outdir, 'quality_control')
        os.makedirs(fqfd, exist_ok=True)
        fqf = self.tools['fastq_quality_filter']
        ft = self.tools['fastx_trimmer']
        qual = self.settings['quality']
        fcl = int(self.settings['forward_cut_length'].strip())
        rcl = int(self.settings['reverse_cut_length'].strip())
        for sample in self.samples:
            # forward read
            infq1 = self.files_to_process[sample][0]
            outfq1 = os.path.join(fqfd, f'{sample}_qc_R1.fastq.gz')
            trim_fr = f'{ft} -Q 33 -f {fcl} -i - |'
            if fcl <= 0:
                trim_fr = ''
            cmd_fr = [
                f'gunzip -c {infq1} | {trim_fr} {fqf}',
                f'-v -Q 33 -q {qual} -p 40 -i - | gzip -c > {outfq1}']
            self.execute_cmd(' '.join(cmd_fr))
            # reverse read
            infq2 = self.files_to_process[sample][1]
            outfq2 = os.path.join(fqfd, f'{sample}_qc_R2.fastq.gz')
            trim_rr = f'{ft} -Q 33 -f {rcl} -i - |'
            if rcl <= 0:
                trim_rr = ''
            cmd_rr = [
                f'gunzip -c {infq2} | {trim_rr} {fqf}',
                f'-v -Q 33 -q {qual} -p 40 -i - | gzip -c > {outfq2}']
            self.execute_cmd(' '.join(cmd_rr))
            self.files_to_process[sample] = [outfq1, outfq2]

    def adapter_trim(self):
        # adapter trim
        if self.settings['adapter_trim'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip adapter trim')
            return 0
        self.logging_title('adapter trim : fastq')
        ca = os.path.join(self.outdir, 'adapter_trim')
        os.makedirs(ca, exist_ok=True)
        fra = self.settings['fra'].strip()
        rra = self.settings['rra'].strip()
        err = self.settings['error_rate']
        ml = self.settings['min_len']
        for sample in self.samples:
            # forward read
            infq1 = self.files_to_process[sample][0]
            outfq1 = os.path.join(ca, f'{sample}_at_R1.fastq.gz')
            a_fr = f'-b {fra}'
            noseq_f = len(fra.strip('ATCG'))
            if len(fra) == 0 or noseq_f > 0:
                a_fr = ''
            cmd_fr = f'cutadapt -e {err} -m {ml} -M {ml} {a_fr} -o {outfq1} {infq1}'
            self.execute_cmd(cmd_fr)
            # reverse read
            infq2 = self.files_to_process[sample][1]
            outfq2 = os.path.join(ca, f'{sample}_at_R2.fastq.gz')
            a_rr = f'-b {rra}'
            noseq_r = len(rra.strip('ATCG'))
            if len(rra) == 0 or noseq_r > 0:
                a_rr = ''
            cmd_rr = f'cutadapt -e {err} -m {ml} -M {ml} {a_rr} -o {outfq2} {infq2}'
            self.execute_cmd(cmd_rr)
            self.files_to_process[sample] = [outfq1, outfq2]

    def fastq_to_fasta(self):
        # fastq to fasta
        if self.settings['fastq_to_fasta'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip fastq_to_fasta')
            return 0
        self.logging_title('fastq to fasta')
        fq2fa = self.tools['fastq_to_fasta']
        bgzip = self.tools['bgzip']
        fasta = os.path.join(self.outdir, 'fasta')
        os.makedirs(fasta, exist_ok=True)
        for sample in self.samples:
            infq1 = self.files_to_process[sample][0]
            infq2 = self.files_to_process[sample][1]
            outfa = os.path.join(fasta, f'{sample}.fasta')
            cmds = [
                f'zcat {infq1} | {fq2fa} -Q 33 |',
                f"sed -e '/^>/s/ /--/g'> {outfa}",
                '&&',
                f'zcat {infq2} | {fq2fa} -Q 33 |',
                f"sed -e '/^>/s/ /--/g' >> {outfa}",
                '&&',
                f'{bgzip} -f {outfa}']
            self.execute_cmd(' '.join(cmds))
            self.files_to_process[sample] = [f'{outfa}.gz']

    def from_fasta(self):
        indir_star = os.path.join(self.indir, '*')
        fastas = set(glob.glob(indir_star))
        self.files_to_process = {}
        outlines = ['#samplename,fasta']
        for fa in fastas:
            if fa.endswith('.fasta.gz'):
                basename = os.path.basename(fa)
                sample = basename[:-9]
                outline = f'{sample},{basename}'
                outlines.append(outline)
                self.files_to_process[sample] = [fa]
            else:
                continue
        with open(os.path.join(self.outdir, 'sample.csv'), 'w') as ini:
            ini.write('\n'.join(outlines))
        self.samples = self.files_to_process.keys()

    def stacks(self):
        if self.settings['stacks'] != 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip stacks')
            return 0
        self.logging_title('ustacks')
        if self.settings['from_fa'] == 'T':
            self.logger.info('start from fasta')
            self.from_fasta()
        # prepare for stacks
        outdir = os.path.join(self.outdir, 'ustacks')
        os.makedirs(outdir, exist_ok=True)
        thre = self.settings['threads']
        ustacks = self.tools['ustacks']
        ustacks_opts = self.settings['ustacks_opts']
        id_sample_list = []
        if len(ustacks_opts.strip()) == 0:
            ustacks_opts = ''
        identifier = 1
        for sample in self.samples:
            infa = self.files_to_process[sample][0]
            cmd = (
                f'{ustacks} -f {infa} -o {outdir} --name {sample} '
                f'-i {identifier} {ustacks_opts} -p {thre}')
            self.execute_cmd(cmd)
            id_sample_list.append(f'{identifier},{sample}')
            outtsv = os.path.join(
                outdir, f'{sample}.tags.tsv.gz')
            self.files_to_process[sample] = [outtsv]
            identifier += 1
        # save sample_id - sample_name list
        id_sample_csv = os.path.join(outdir, 'id_sample.csv')
        with open(id_sample_csv, 'w') as isc:
            isc.write('\n'.join(id_sample_list))

    def extract_consensus_fasta(self):
        if self.settings['extract_consensus_fasta'] != 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip extracting consensus sequences')
            return 0
        self.logging_title('extract consensus fasta')
        outdir = os.path.join(self.outdir, 'consensus_fasta')
        os.makedirs(outdir, exist_ok=True)
        for sample in self.samples:
            tagsgz = self.files_to_process[sample][0]
            outlines = []
            with gzip.open(tagsgz, 'rt') as tags:
                line = tags.readline()
                while line:
                    data = [i.strip() for i in line.split(sep='\t')]
                    if len(data) < 6:
                        line = tags.readline()
                        continue
                    if data[2] != 'consensus':
                        line = tags.readline()
                        continue
                    sample_id = data[0]
                    loc_id = data[1]
                    seq = data[5]
                    fa = f'>{sample}__{sample_id}__{loc_id}\n{seq}'
                    outlines.append(fa)
                    line = tags.readline()
            outfasta = os.path.join(outdir, f'{sample}.fasta')
            with open(outfasta, 'w') as of:
                of.write('\n'.join(outlines))
            self.files_to_process[sample] = [outfasta]

    def from_consensus_fasta(self):
        indir_star = os.path.join(self.indir, '*')
        fastas = set(glob.glob(indir_star))
        self.files_to_process = {}
        outlines = ['#samplename,fasta']
        for fa in fastas:
            if fa.endswith('.fasta'):
                basename = os.path.basename(fa)
                sample = basename[:-6]
                outline = f'{sample},{basename}'
                outlines.append(outline)
                self.files_to_process[sample] = [fa]
            else:
                continue
        with open(os.path.join(self.outdir, 'sample.csv'), 'w') as ini:
            ini.write('\n'.join(outlines))
        self.samples = self.files_to_process.keys()

    def blastn_search(self):
        self.logging_title('blast search')
        if self.settings['from_cons_fa'] == 'T':
            self.from_consensus_fasta()
        db = self.settings['database']
        blastn = self.tools['blastn']
        outdir = os.path.join(self.outdir, 'blast_results')
        os.makedirs(outdir, exist_ok=True)
        for sample in self.samples:
            fa = self.files_to_process[sample][0]
            outtsv = os.path.join(outdir, f'{sample}.tsv')
            cmd = (
                f'{blastn} -db {db} -query {fa} '
                f'-outfmt 6 > {outtsv}')
            self.execute_cmd(cmd)

    def make_blastdb(self):
        self.logging_title('make blastdb')
        if self.settings['from_cons_fa'] == 'T':
            self.from_consensus_fasta()
        db = self.settings['database']
        mkblastdb = self.tools['makeblastdb']
        blastdir = os.path.dirname(db)
        os.makedirs(blastdir, exist_ok=True)
        merged_fasta = db + '.fasta'
        # make merged fasta
        fastas = ' '.join(
            [i[0] for i in self.files_to_process.values()])
        cmd = f'cat {fastas} > {merged_fasta}'
        self.execute_cmd(cmd)
        # makeblastdb
        cmd = (
            f'{mkblastdb} -in {merged_fasta} '
            f'-out {db} -dbtype nucl -parse_seqids')
        self.execute_cmd(cmd)

    def run(self):
        self.copy_conf()
        self.inputdir_to_ini()
        self.quality_control()
        self.adapter_trim()
        self.fastq_to_fasta()
        self.stacks()
        self.extract_consensus_fasta()
        mode = self.settings['mode']
        if mode == 'search':
            self.blastn_search()
        elif mode == 'makedb':
            self.make_blastdb()
        else:
            msg = (
                'no such a mode : {mode}'
                '\nchoose mode from [search, makedb]')
            self.logger.info(msg)
            sys.exit(1)


def argument_parser():
    parser = argparse.ArgumentParser(description="database accessor")
    parser.add_argument(
        '-c', '--conf', dest='conf', type=str,
        help='config file', required=True)
    return parser.parse_args()


def main():
    args = argument_parser()
    MS = DatabaseAccessor(conf=args.conf)
    MS.run()
    # analysis complete flag
    acf(config=args.conf)


if __name__ == '__main__':
    main()
