import argparse
import configparser
import logging
import os
import sys
import subprocess
from subprocess import PIPE


class ACF:
    def __init__(
            self, outdir='outdir',
            acftext='analysis_complete.txt',
            acfmessage=None):
        self.outdir = outdir
        self.acf = acftext
        self.acfmessage = acfmessage
        self.acftext = os.path.join(outdir, acftext)

    def setup_logger(self, name=None):
        os.makedirs(self.outdir, exist_ok=True)
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(self.acftext)
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

    def get_outdir_from_config(
            self, config='config.conf',
            section='settings', outdir='output_dir'):
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(config, 'UTF-8')
        sec = dict(conf.items(section))
        self.outdir = sec[outdir]
        self.acftext = os.path.join(self.outdir, self.acf)

    def tree_outdir(self, treefile='tree.txt'):
        outdir = self.outdir
        cmd = f'cd {outdir} && tree . > {treefile}'
        self.execute_cmd(cmd)

    def write_analysis_complete_flag(self):
        message = 'ANCAT_ANALYSIS_COMPLETED'
        if self.acfmessage is not None:
            message = self.acfmessage
        self.logger.info(message)

    def run(self):
        self.setup_logger()
        self.tree_outdir()
        self.write_analysis_complete_flag()


def analysis_complete_flag(
        config=None, outdir='output_dir',
        section='settings',
        acftext='analysis_complete.txt', tree='tree.txt'):
    Instance = ACF(outdir=outdir, acftext=acftext)
    if config is not None:
        Instance.get_outdir_from_config(
            config=config, section=section, outdir=outdir)
    Instance.run()


if __name__ == '__main__':
    analysis_complete_flag(
        config='out.conf', section='settings', outdir='output_dir')
