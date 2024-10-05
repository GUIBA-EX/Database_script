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
import shutil


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
        print(f"当前工作目录: {os.getcwd()}")
        print(f"配置文件路径: {self.conf}")
        print(f"配置文件是否存在: {os.path.exists(self.conf)}")
        
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(self.conf, 'UTF-8')
        
        print(f"读取到的配置部分: {conf.sections()}")
        
        if 'tools' not in conf.sections():
            print("警告：'tools' 部分不在配置文件中")
        
        try:
            self.tools = dict(conf.items('tools'))
            print(f"成功读取 'tools' 部分: {self.tools}")
        except configparser.NoSectionError:
            print("错误：无法读取 'tools' 部分")
            raise
        
        self.settings = dict(conf.items('settings'))
        self.log = os.path.join(self.settings['output_dir'], 'log.txt')
        # 添加参数验证
        required_settings = ['input_dir', 'output_dir', 'quality_control', 'adapter_trim', 'mode']
        for setting in required_settings:
            if setting not in self.settings:
                raise ValueError(f"配置文件中缺少必要的设置: {setting}")

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
        if self.settings['quality_control'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('skip quality control')
            return 0
        
        self.logging_title('quality control : fastq')
        fqfd = os.path.join(self.outdir, 'quality_control')
        os.makedirs(fqfd, exist_ok=True)
        
        fastp = self.tools.get('fastp')
        if not fastp:
            self.logger.error("fastp路径未在配置中找到")
            return 1

        qual = self.settings.get('quality', '20')  # 默认值为20
        fcl = int(self.settings.get('forward_cut_length', '0').strip())
        rcl = int(self.settings.get('reverse_cut_length', '0').strip())
        threads = self.settings.get('threads', '4')  # 从配置文件读取线程数，默认为4
        
        for sample in self.samples:
            self.logger.info(f"Processing sample: {sample}")
            infq1 = self.files_to_process[sample][0]
            infq2 = self.files_to_process[sample][1]
            outfq1 = os.path.join(fqfd, f'{sample}_qc_R1.fastq.gz')
            outfq2 = os.path.join(fqfd, f'{sample}_qc_R2.fastq.gz')
            
            cmd = [
                f'{fastp}',
                f'-i {infq1} -I {infq2}',
                f'-o {outfq1} -O {outfq2}',
                f'-q {qual}',
                f'-f {fcl} -F {rcl}',
                '--cut_right',
                f'--thread {threads}',
                f'--json {os.path.join(fqfd, f"{sample}_fastp.json")}',
                f'--html {os.path.join(fqfd, f"{sample}_fastp.html")}'
            ]
            
            stdout, stderr = self.execute_cmd(' '.join(cmd))
            if stderr:
                self.logger.warning(f"fastp for {sample} produced warnings: {stderr}")
            
            self.files_to_process[sample] = [outfq1, outfq2]
            self.logger.info(f"Finished processing sample: {sample}")

        self.logger.info("Quality control completed for all samples")
        return 0

    def adapter_trim(self):
        # adapter trim
        if self.settings['adapter_trim'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('跳过接头修剪')
            return 0
        
        self.logging_title('接头修剪：fastq')
        ca = os.path.join(self.outdir, 'adapter_trim')
        os.makedirs(ca, exist_ok=True)
        fra = self.settings['fra'].strip()
        rra = self.settings['rra'].strip()
        err = self.settings['error_rate']
        ml = self.settings['min_len']
        
        for sample in self.samples:
            infq1 = self.files_to_process[sample][0]
            infq2 = self.files_to_process[sample][1]
            outfq1 = os.path.join(ca, f'{sample}_at_R1.fastq.gz')
            outfq2 = os.path.join(ca, f'{sample}_at_R2.fastq.gz')
            
            # 构 fastp 命令
            cmd = f"fastp -i {infq1} -I {infq2} -o {outfq1} -O {outfq2} " \
                  f"--adapter_sequence {fra} --adapter_sequence_r2 {rra} " \
                  f"--cut_right --cut_window_size 4 --cut_mean_quality 20 " \
                  f"--length_required 36 --length_limit 150 " \
                  f"--average_qual 19 --thread 8 " \
                  f"--json {os.path.join(ca, f'{sample}_fastp.json')} " \
                  f"--html {os.path.join(ca, f'{sample}_fastp.html')}"
            
            # 执行命令
            self.execute_cmd(cmd)
            
            # 更新文件路径
            self.files_to_process[sample] = [outfq1, outfq2]
            
            # 检查输出文件
            if os.path.getsize(outfq1) == 0 or os.path.getsize(outfq2) == 0:
                self.logger.warning(f"{sample} 的接头修剪后文件为空，将使用原始文件")
                self.files_to_process[sample] = [infq1, infq2]
        
        self.logger.info('所有样本的接头修剪已完成')

    def fastq_to_fasta(self):
        if self.settings['fastq_to_fasta'] != 'T' or \
                self.settings['from_fa'] == 'T' or \
                self.settings['from_cons_fa'] == 'T':
            self.logger.info('跳过 fastq_to_fasta')
            return 0
        
        self.logging_title('fastq 转换为 fasta')
        seqkit = self.tools.get('seqkit')
        if not seqkit:
            self.logger.error("seqkit 路径未在配置中找到")
            return 1

        fasta_dir = os.path.join(self.outdir, 'fasta')
        os.makedirs(fasta_dir, exist_ok=True)
        
        for sample in self.samples:
            infq1 = self.files_to_process[sample][0]
            infq2 = self.files_to_process[sample][1]
            outfa = os.path.join(fasta_dir, f'{sample}.fasta.gz')
            
            cmd = [
                f'{seqkit} fq2fa',
                f'-o {outfa}',
                f'{infq1} {infq2}'
            ]
            
            stdout, stderr = self.execute_cmd(' '.join(cmd))
            if stderr:
                self.logger.warning(f"{sample} 的 seqkit 处理产生警告: {stderr}")
            
            self.files_to_process[sample] = [outfa]
        
        self.logger.info("所有样本的 FASTQ 到 FASTA 转换已完成")
        return 0

    def uniform_sequence_length(self, input_file, output_file, target_length):
        cmd = f"seqkit seq -m {target_length} -M {target_length} {input_file} | gzip > {output_file}"
        self.execute_cmd(cmd)

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
        self.logging_title('ustacks')
        outdir = os.path.join(self.outdir, 'ustacks')
        os.makedirs(outdir, exist_ok=True)
        ustacks = 'ustacks'  # 直接使用 ustacks 命令，因为它已在 PATH 中
        ustacks_opts = self.settings['ustacks_opts']
        thre = self.settings['threads']
        identifier = 1
        id_sample_list = []
        
        for sample in self.samples:
            infa = self.files_to_process[sample]  # 这应该是一个字符串，而不是列
            cmd = (
                f'{ustacks} -f {infa} -o {outdir} --name {sample} '
                f'-i {identifier} {ustacks_opts} -p {thre}'
            )
            
            self.logger.info(f"Running ustacks for sample: {sample}")
            self.logger.info(f"Command: {cmd}")
            
            retry_count = 0
            max_retries = 3
            while retry_count < max_retries:
                self.execute_cmd(cmd)
                if self.check_ustacks_output(sample, outdir):
                    self.logger.info(f"ustacks successfully processed sample: {sample}")
                    break
                retry_count += 1
                self.logger.warning(f"ustacks processing failed for {sample}, retry {retry_count}/{max_retries}")
            
            if retry_count == max_retries:
                self.logger.error(f"ustacks processing failed for {sample} after {max_retries} attempts, skipping this sample")
                continue
            
            id_sample_list.append(f'{identifier},{sample}')
            outtsv = os.path.join(outdir, f'{sample}.tags.tsv.gz')
            if os.path.isfile(outtsv):
                self.files_to_process[sample] = outtsv
            else:
                self.logger.error(f"ustacks 输出文件不存在: {outtsv}")
            identifier += 1
        
        # save sample_id - sample_name list
        id_sample_csv = os.path.join(outdir, 'id_sample.csv')
        with open(id_sample_csv, 'w') as isc:
            isc.write('\n'.join(id_sample_list))

    def extract_consensus_fasta(self):
        self.logging_title('extract consensus fasta')
        consensus_dir = os.path.join(self.outdir, 'consensus_fasta')
        os.makedirs(consensus_dir, exist_ok=True)
        
        for sample in self.samples:
            input_file = self.files_to_process[sample]
            output_file = os.path.join(consensus_dir, f"{sample}.fa")
            cmd = (f"zcat {input_file} | awk '{{if($3==\"consensus\") {{printf \">{sample}_%d\\n%s\\n\", NR, $4}}}}' > {output_file}")
            self.logger.info(f"[cmd] {cmd}")
            self.execute_cmd(cmd)
            
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                self.logger.info(f"成功生成共识序列文件: {output_file}")
                self.files_to_process[sample] = output_file
            else:
                self.logger.error(f"生成共识序列文件失败: {output_file}")
            
            # 检查生成的 FASTA 文件
            self.logger.info(f"检查生成的 FASTA 文件: {output_file}")
            self.execute_cmd(f"head -n 10 {output_file}")

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
        merged_fasta = os.path.join(blastdir, "merged_consensus.fasta")
        
        # 准备 FASTA 文件列表
        fastas = [self.files_to_process[sample] for sample in self.samples]
        if not fastas:
            self.logger.error("没有找到 FASTA 文件来创建 BLAST 数据库")
            return
        
        # 使用 cat 命令合并 FASTA 文件
        cat_cmd = f"cat {' '.join(fastas)} > {merged_fasta}"
        self.logger.info(f"合并 FASTA 文件到: {merged_fasta}")
        self.execute_cmd(cat_cmd)
        
        if not os.path.exists(merged_fasta) or os.path.getsize(merged_fasta) == 0:
            self.logger.error(f"合并的 FASTA 文件 {merged_fasta} 不存在或为空")
            return
        
        # 检查合并后的 FASTA 文件格式
        self.logger.info(f"检查合并后的 FASTA 文件格式")
        self.execute_cmd(f"head -n 10 {merged_fasta}")
        
        # 创建 BLAST 数据库
        db_name = os.path.join(blastdir, "acropora_db")
        cmd = (
            f'{mkblastdb} -in {merged_fasta} '
            f'-out {db_name} -dbtype nucl -parse_seqids'
        )
        self.logger.info(f"创建 BLAST 数据库: {db_name}")
        self.execute_cmd(cmd)
        
        # 修改检查逻辑
        db_files = [f"{db_name}.{ext}" for ext in ['nhr', 'nin', 'nsq']]
        if all(os.path.exists(file) for file in db_files):
            self.logger.info(f"BLAST 数据库创建成功: {db_name}")
        else:
            self.logger.error(f"BLAST 数据库创建失败: {db_name}")
            missing_files = [file for file in db_files if not os.path.exists(file)]
            self.logger.error(f"缺少以下文件: {', '.join(missing_files)}")
        
        # 列出 blastdb 目录中的文件
        self.logger.info("BLAST 数据库目录内容:")
        self.execute_cmd(f"ls -l {blastdir}")

    def check_file_format(self, file_path):
        with gzip.open(file_path, 'rt') as f:
            first_char = f.read(1)
            if first_char == '>':
                return 'fasta'
            elif first_char == '@':
                return 'fastq'
            else:
                return 'unknown'

    def process_step(self, input_file, output_file, process_func):
        try:
            process_func(input_file, output_file)
            if os.path.getsize(output_file) == 0:
                logging.warning(f"{output_file} 为空，使用输入文件继续")
                shutil.copy(input_file, output_file)
        except Exception as e:
            logging.error(f"处理 {input_file} 时出错: {str(e)}")
            shutil.copy(input_file, output_file)

    def check_ustacks_output(self, sample, ustacks_output_dir):
        expected_files = [
            f"{sample}.alleles.tsv.gz",
            f"{sample}.snps.tsv.gz",
            f"{sample}.tags.tsv.gz"
        ]
        for file in expected_files:
            full_path = os.path.join(ustacks_output_dir, file)
            if not os.path.exists(full_path):
                self.logger.error(f"Missing ustacks output file for {sample}: {file}")
                return False
        return True

    def run(self):
        try:
            self.copy_conf()
            self.inputdir_to_ini()
            self.quality_control()
            self.adapter_trim()
            self.fastq_to_fasta()
            # 在 FASTA 转换后使用
            for sample in self.samples:
                fasta_file = f"/home/data/acropora_db/fasta/{sample}.fasta.gz"
                uniform_fasta = f"/home/data/acropora_db/fasta/{sample}_uniform.fasta.gz"
                self.uniform_sequence_length(fasta_file, uniform_fasta, 80)  # 假设目标长度为80
                self.files_to_process[sample] = uniform_fasta
            self.stacks()
            self.logger.info("stacks 处理完成")
            self.logger.info(f"当前 files_to_process: {self.files_to_process}")
            self.extract_consensus_fasta()
            mode = self.settings['mode']
            if mode == 'search':
                self.blastn_search()
            elif mode == 'makedb':
                self.make_blastdb()
            else:
                raise ValueError(f'无效的模式: {mode}. 请从 [search, makedb] 中选择')
        except ValueError as ve:
            self.logger.error(f"配置错误: {str(ve)}")
            sys.exit(1)
        except IOError as ioe:
            self.logger.error(f"IO错误: {str(ioe)}")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"运行过程中发生错误: {str(e)}")
            self.logger.exception("详细错误信息:")
        finally:
            self.logger.info("处理完成")


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