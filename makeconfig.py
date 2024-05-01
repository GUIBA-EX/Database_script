import argparse
import configparser
import os
import sys
import logging
import subprocess
from subprocess import PIPE


class MakeConfig:
    def __init__(self):
        self.argument_parser()
        self.setup_logger()

    def argument_parser(self):
        parser = argparse.ArgumentParser(description="migseq config maker")
        parser.add_argument(
            '-t', '--tmp_conf', dest='tmp', type=str,
            help='template config file to get path to tools', required=False,
            default=os.path.join(os.path.dirname(__file__), 'tmp.conf'))
        parser.add_argument(
            '-oc', '--output_conf', dest='oc', type=str,
            help='output config file', required=False,
            default='out.conf')
        parser.add_argument(
            '-g', '--genome', dest='genome', type=str,
            help='path to genome fasta',
            required=False, default='None')
        parser.add_argument(
            '-popmap', '--population_map', dest='popmap', type=str,
            help='popmap file',
            default='',
            required=False)
        parser.add_argument(
            '-init', '--initialization', dest='init',
            action='store_true',
            help='check and make index of genome for bwa execution',
            required=False)
        parser.add_argument(
            '-outdir', '--output_dir', dest='outdir', type=str,
            help='directory to output results',
            required=True, default=None)
        parser.add_argument(
            '-indir', '--input_dir', dest='indir', type=str,
            help='directory which contains input fastqs',
            required=True, default=None)
        parser.add_argument(
            '-qc', '--quality_control', dest='qc', type=str,
            choices=['T', 'F'], default='F',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-at', '--adapter_trim', dest='at', type=str,
            choices=['T', 'F'], default='F',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-rp', '--repair', dest='rp', type=str,
            choices=['T', 'F'], default='F',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-map', '--mapping', dest='map', type=str,
            choices=['T', 'F'], default='F',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-st', '--stacks', dest='st', type=str,
            choices=['T', 'F'], default='T',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-uop', '--use_original_popmap', dest='uop', type=str,
            choices=['T', 'F'], default='F',
            help='choose from T(rue) and F(alse)', required=False)
        parser.add_argument(
            '-q', '--quality', dest='q', type=str, default='30',
            help='quality cutoff : default : 30', required=False)
        parser.add_argument(
            '-frcl', '--forward_read_cut_lenght',
            dest='frcl', type=str,
            default='0',
            help='forward read cut lenght',
            required=False)
        parser.add_argument(
            '-rrcl', '--reverse_read_cut_lenght',
            dest='rrcl', type=str,
            default='0',
            help='reverse read cut lenght',
            required=False)
        parser.add_argument(
            '-fra', '--forward_read_adapter',
            dest='fra', type=str,
            default='None',
            help="forward read adapter, example : ATCGGA/CCATGA",
            required=False)
        parser.add_argument(
            '-rra', '--reverse_read_adapter',
            dest='rra', type=str,
            default='None',
            help="reverse read adapter, example : ATCGGA/CCATGA",
            required=False)
        parser.add_argument(
            '-aer', '--adapter_error_rate', dest='aer',
            type=str, default='0.05',
            help="adapter error rate : default 0.05", required=False)
        parser.add_argument(
            '-r', '--min_samples_per_pop',
            dest='r', type=str,
            default='0.7',
            help="minimum samples per population : default 0.7",
            required=False)
        parser.add_argument(
            '-R', '--min_samples_overall',
            dest='R', type=str,
            default='0.1',
            help="minimum samples overall : default 0.1",
            required=False)
        parser.add_argument(
            '-p', '--min_populations',
            dest='p', type=str,
            default='1',
            help="minimum populations : default 1",
            required=False)
        parser.add_argument(
            '-mmaf', '--min_maf',
            dest='mmaf', type=str,
            default='0.1',
            help="min-maf : default 0.1",
            required=False)
        parser.add_argument(
            '-mmac', '--min_mac',
            dest='mmac', type=str,
            default='0.1',
            help="min-mac : default 0.1",
            required=False)
        parser.add_argument(
            '-moh', '--max_obs_het',
            dest='moh', type=str,
            default='0.99',
            help="maximum observed heterozygosity : default 0.99",
            required=False)
        parser.add_argument(
            '-thre', '--threads', dest='thre', type=int,
            default=2,
            help="num of threads : default 2", required=False)
        parser.add_argument(
            '-log', '--log', dest='log', type=str,
            default='makeconfig_log.txt',
            help="makeconfig logfile : default makeconfig_log.txt",
            required=False)
        self.args = parser.parse_args()

    def setup_logger(self, name=None):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(self.args.log)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s %(funcName)s :\n%(message)s')
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
        if len(stdout) > 1:
            self.logger.info(stdout)
        if len(stderr) > 1:
            self.logger.info(stderr)
        return stdout, stderr

    def parse_config(self):
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(self.args.tmp, 'UTF-8')
        self.tools = dict(conf.items('tools'))

    def check_input_info(self):
        # check tools
        for tool in self.tools.keys():
            ToF = os.path.exists(self.tools[tool])
            if ToF is False:
                self.logger.info(
                    f'{tool} : no such a tool : {self.tools[tool]}')
                sys.exit(1)
        # check genome
        genome = self.args.genome
        ToF = os.path.isfile(genome)
        if self.args.map == 'T' and ToF is not True:
            self.logger.info(f'no such a genome file : {genome}')
            sys.exit(1)

    def initialize(self):
        if self.args.map == 'F' or self.args.init is not True:
            self.logger.info('skip checking or making index')
            return 0
        # index file check
        genome = self.args.genome
        # samtools
        samtools_index = os.path.isfile(genome + '.fai')
        # bwa
        for ind in ['.amb', '.ann', '.bwt', '.sa', '.pac']:
            bwa_index = os.path.isfile(genome + ind)
            if bwa_index is False:
                break
        # prepare for indexing
        samtools = self.tools['samtools']
        bwa = self.tools['bwa']
        # index command
        cmd_samtools = samtools + ' faidx {genome}'
        cmd_bwa = bwa + ' index {genome}'
        exec_dict = {
            'samtools': {'exists': samtools_index, 'cmd': cmd_samtools},
            'bwa': {'exists': bwa_index, 'cmd': cmd_bwa}}
        # execute command
        for tool in exec_dict.keys():
            Exists = exec_dict[tool]['exists']
            cmd = exec_dict[tool]['cmd']
            if Exists:
                self.logger.info(
                    f'Index for {tool} is already exists. Skip indexing')
            else:
                self.execute_cmd(cmd.format(genome=genome))

    def make_config(self):
        lines = []
        # make tools paths
        lines.append('[tools]')
        for tool in self.tools.keys():
            lines.append(tool + ' = ' + self.tools[tool])
        # make settings
        lines.append('[settings]')
        # append genome
        lines.append('genome_index = ' + self.args.genome)
        # popmap
        lines.append('popmap = ' + self.args.popmap)
        # input _ output dir
        lines.append('input_dir = ' + self.args.indir)
        lines.append('output_dir = ' + self.args.outdir)
        # settings for pipeline
        lines.append('quality_control = ' + self.args.qc)
        lines.append('adapter_trim = ' + self.args.at)
        lines.append('repair = ' + self.args.rp)
        lines.append('mapping = ' + self.args.map)
        lines.append('stacks = ' + self.args.st)
        # settings for each process
        if self.args.qc == 'T':
            lines.append('quality = ' + str(self.args.q))
            lines.append('forward_cut_length = ' + str(self.args.frcl))
            lines.append('reverse_cut_length = ' + str(self.args.rrcl))
        if self.args.at == 'T':
            lines.append('fra = ' + self.args.fra)
            lines.append('rra = ' + self.args.rra)
            lines.append('error_rate = ' + self.args.aer)
        if self.args.st == 'T':
            lines.append('min_samples_per_pop = ' + str(self.args.r))
            lines.append('min_samples_overall = ' + str(self.args.R))
            lines.append('min_pops = ' + str(self.args.p))
            lines.append('min_maf = ' + str(self.args.mmaf))
            lines.append('min_mac = ' + str(self.args.mmac))
            lines.append('max_obs_het = ' + str(self.args.moh))
            lines.append('use_original_popmap = ' + str(self.args.uop))
        lines.append('threads = ' + str(self.args.thre))
        with open(self.args.oc, 'w') as oc:
            oc.write('# command : ' + ' '.join(sys.argv) + '\n')
            oc.write('\n'.join(lines))

    def run(self):
        self.parse_config()
        self.check_input_info()
        self.initialize()
        self.make_config()


def main():
    MC = MakeConfig()
    MC.run()


if __name__ == '__main__':
    main()
