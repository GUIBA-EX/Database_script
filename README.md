# CoralDatabaseBuilder
这是一个使用Migseq采集的数据进行BLAST数据库构建及搜索的程序。

## 使用docker安装
```
# build image
docker image build -t classifier .
# run container
docker run -itd -v $PWD:/home/data --name coral_classifier classifier bash

# exec container
docker exec -it coral_classifier bash
```

## Usage : database_accessor
fastq => migseq预处理 => ustacks => fasta提取 => blast搜索 or 物种序列登录
```
usage: database_accessor.py [-h] -c CONF

database accessor

options:
  -h, --help            show this help message and exit
  -c CONF, --conf CONF  config file
```
例子
```
python database_accessor.py -c config.conf
python database_accessor.py -c mkdb_acropora.conf
python database_accessor.py -c eDNA.conf
```

## config file

```
[tools]
fastp=fastp                         # 指定 fastp 可执行程序的路径
seqkit=seqkit                       # 指定 seqkit 可执行程序的路径
bgzip=bgzip                         # 指定 bgzip 可执行程序的路径
ustacks=ustacks                     # 指定 ustacks 可执行程序的路径
makeblastdb=makeblastdb             # 指定 makeblastdb 可执行程序的路径
blastn=blastn                       # 指定 blastn 可执行程序的路径

[settings]
mode=makedb                         # 指定运行模式，创建数据库为 makedb，搜索为 search
input_dir=/home/data/acropora/      # 输入文件所在目录
output_dir=/home/data/acropora_db/  # 输出文件保存目录
database=/home/data/acropora_db/blastdb/  # 指定数据库文件路径

threads=8                           # 指定使用的线程数

quality_control=T                   # 是否执行质量控制 (T: 是, F: 否)
quality=30                          # 质量控制时使用的质量阈值
reverse_cut_length=15               # 用于切割 reverse read 的长度
forward_cut_length=0                # 用于切割 forward read 的长度

adapter_trim=T                      # 是否执行适配器去除 (T: 是, F: 否)
error_rate=0.05                     # 适配器的错误允许率
fra=GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  # forward read 中需要去除的适配器序列
rra=CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC      # reverse read 中需要去除的适配器序列
min_len=80                          # 最小序列长度

fastq_to_fasta=T                    # 是否将 fastq 转换为 fasta 格式 (T: 是, F: 否)
from_fa=F                           # 是否从 fasta 格式开始分析 (T: 是, F: 否)
stacks=T                            # 是否使用 stacks 进行分析 (T: 是, F: 否)
ustacks_opts=-m 3 -M 2              # ustacks 的运行参数
extract_consensus_fasta=T           # 是否从 ustacks 结果中提取共识序列生成 fasta 文件 (T: 是, F: 否)
from_cons_fa=F                      # 是否从 ustacks 结果生成的 fasta 文件开始分析 (T: 是, F: 否)
```

| Category  | Command                | Function                                                                                                                                  |
|:----------|:-----------------------|:------------------------------------------------------------------------------------------------------------------------------------------|
| **Tools**     | **[fastp]**                | Specify the path to the `fastp` executable                                                                                               |
|           | **[seqkit]**               | Specify the path to the `seqkit` executable                                                                                              |
|           | **[bgzip]**                | Specify the path to the `bgzip` executable                                                                                                |
|           | **[ustacks]**              | Specify the path to the `ustacks` executable                                                                                              |
|           | **[makeblastdb]**          | Specify the path to the `makeblastdb` executable                                                                                          |
|           | **[blastn]**               | Specify the path to the `blastn` executable                                                                                               |
| **Settings**  | **[mode]**                 | Set the running mode: `makedb` for creating a database, `search` for searching                                                            |
|           | **[input_dir]**            | Specify the directory containing input files                                                                                              |
|           | **[output_dir]**           | Specify the directory to save output files                                                                                                |
|           | **[database]**             | Specify the path to the database file                                                                                                     |
|           | **[threads]**              | Specify the number of threads to use                                                                                                      |
|           | **[quality_control]**      | Determine whether to perform quality control (`T`: Yes, `F`: No)                                                                          |
|           | **[quality]**              | Set the quality threshold for quality control                                                                                             |
|           | **[reverse_cut_length]**   | Set the length to trim from reverse reads                                                                                                 |
|           | **[forward_cut_length]**   | Set the length to trim from forward reads                                                                                                 |
|           | **[adapter_trim]**         | Determine whether to remove adapters (`T`: Yes, `F`: No)                                                                                  |
|           | **[error_rate]**           | Set the allowed error rate for adapter trimming                                                                                           |
|           | **[fra]**                  | Specify the adapter sequence to remove from forward reads                                                                                 |
|           | **[rra]**                  | Specify the adapter sequence to remove from reverse reads                                                                                 |
|           | **[min_len]**              | Set the minimum sequence length after trimming                                                                                            |
|           | **[fastq_to_fasta]**       | Determine whether to convert FASTQ files to FASTA format (`T`: Yes, `F`: No)                                                              |
|           | **[from_fa]**              | Determine whether to start analysis from FASTA files (`T`: Yes, `F`: No)                                                                  |
|           | **[stacks]**               | Determine whether to use Stacks for analysis (`T`: Yes, `F`: No)                                                                          |
|           | **[ustacks_opts]**         | Specify parameters for running `ustacks`                                                                                                  |
|           | **[extract_consensus_fasta]** | Determine whether to extract consensus sequences from `ustacks` results to generate a FASTA file (`T`: Yes, `F`: No)                   |
|           | **[from_cons_fa]**         | Determine whether to start analysis from the FASTA file generated by `ustacks` results (`T`: Yes, `F`: No)                                |



### input_dir

输入文件需要使用pair-end测序文件，命名遵循以下规则

forward read：``[sample_name]_L001_R1_001.fastq.gz``

reverse read：``[sample_name]_L001_R2_001.fastq.gz``

#### files in input_dir (makedb from fastq)

```
ls /home/data/indir-coral/ref/
CJ_L001_R1_001.fastq.gz        Pelatius_L001_R2_001.fastq.gz
CJ_L001_R2_001.fastq.gz        Pkonojoi_L001_R1_001.fastq.gz
Pelatius_L001_R1_001.fastq.gz  Pkonojoi_L001_R2_001.fastq.gz
```

#### files in input_dir (search from fastq)

```
ls /home/data/indir-coral/Coralliidae
CJ-G1342_L001_R1_001.fastq.gz           Pelatius-IOU010_L001_R2_001.fastq.gz
CJ-G1342_L001_R2_001.fastq.gz           Pelatius-S2201_L001_R1_001.fastq.gz
CJ-G1343_L001_R1_001.fastq.gz           Pelatius-S2201_L001_R2_001.fastq.gz
CJ-G1343_L001_R2_001.fastq.gz           Pelatius-TNG002_L001_R1_001.fastq.gz
CJ-G1344_L001_R1_001.fastq.gz           Pelatius-TNG002_L001_R2_001.fastq.gz
CJ-G1344_L001_R2_001.fastq.gz           Pelatius-TNG005_L001_R1_001.fastq.gz
CJ-G1345_L001_R1_001.fastq.gz           Pelatius-TNG005_L001_R2_001.fastq.gz
CJ-G1345_L001_R2_001.fastq.gz           Pelatius-TNG006_L001_R1_001.fastq.gz
CJ-G1346_L001_R1_001.fastq.gz           Pelatius-TNG006_L001_R2_001.fastq.gz
CJ-G1346_L001_R2_001.fastq.gz           Pkonojoi-IOU002_L001_R1_001.fastq.gz
CJ-S2201_L001_R1_001.fastq.gz           Pkonojoi-IOU002_L001_R2_001.fastq.gz
CJ-S2201_L001_R2_001.fastq.gz           Pkonojoi-IOU003_L001_R1_001.fastq.gz
CJ-S2202_L001_R1_001.fastq.gz           Pkonojoi-IOU003_L001_R2_001.fastq.gz
CJ-S2202_L001_R2_001.fastq.gz           Pkonojoi-IOU004_L001_R1_001.fastq.gz
CJ-S2203_L001_R1_001.fastq.gz           Pkonojoi-IOU004_L001_R2_001.fastq.gz
CJ-S2203_L001_R2_001.fastq.gz           Pkonojoi-IOU007_L001_R1_001.fastq.gz
CJ-S2204_L001_R1_001.fastq.gz           Pkonojoi-IOU007_L001_R2_001.fastq.gz
CJ-S2204_L001_R2_001.fastq.gz           Pkonojoi-TAK001_L001_R1_001.fastq.gz
CJ-S2205_L001_R1_001.fastq.gz           Pkonojoi-TAK001_L001_R2_001.fastq.gz
CJ-S2205_L001_R2_001.fastq.gz           Pkonojoi-TAK003_L001_R1_001.fastq.gz
Pelatius-IOU005_L001_R1_001.fastq.gz    Pkonojoi-TAK003_L001_R2_001.fastq.gz
Pelatius-IOU005_L001_R2_001.fastq.gz    Pkonojoi-TAK005_L001_R1_001.fastq.gz
Pelatius-IOU009-2_L001_R1_001.fastq.gz  Pkonojoi-TAK005_L001_R2_001.fastq.gz
Pelatius-IOU009-2_L001_R2_001.fastq.gz  Pkonojoi-TAK006_L001_R1_001.fastq.gz
Pelatius-IOU010_L001_R1_001.fastq.gz    Pkonojoi-TAK006_L001_R2_001.fastq.gz
```

### output_dir

- adapter_trim : 存放通过fastx_trimmer切除以及通过cutadapt去除adapter后的fastq文件的目录
- quality_control : 存放通过cutadapt进行质量控制后的fastq文件的目录
- fasta : 存放将fastq转换为fasta后的文件的目录
- ustacks : 存放通过ustacks处理fasta文件的结果的目录
- consensus_fasta : 存放ustacks提取出的判断为共识序列的reads并保存为fasta格式的目录
- blast_results : 存放blast结果的目录（仅在搜索模式下使用，输出格式为outfmt 6的tsv文件）

#### output_dir (from fastq)

```
outdir-coral-mkdb-from-docker
├── adapter_trim
│   ├── CJ_at_R1.fastq.gz
│   ├── CJ_at_R2.fastq.gz
│   ├── Pelatius_at_R1.fastq.gz
│   ├── Pelatius_at_R2.fastq.gz
│   ├── Pkonojoi_at_R1.fastq.gz
│   └── Pkonojoi_at_R2.fastq.gz
├── analysis_complete.txt
├── blast_results
│   ├── CJ.tsv
│   ├── Pelatius.tsv
│   └── Pkonojoi.tsv
├── consensus_fasta
│   ├── CJ.fasta
│   ├── Pelatius.fasta
│   └── Pkonojoi.fasta
├── fasta
│   ├── CJ.fasta.gz
│   ├── Pelatius.fasta.gz
│   └── Pkonojoi.fasta.gz
├── log.txt
├── mkdb_tmp.conf
├── quality_control
│   ├── CJ_qc_R1.fastq.gz
│   ├── CJ_qc_R2.fastq.gz
│   ├── Pelatius_qc_R1.fastq.gz
│   ├── Pelatius_qc_R2.fastq.gz
│   ├── Pkonojoi_qc_R1.fastq.gz
│   └── Pkonojoi_qc_R2.fastq.gz
├── sample.csv
├── tree.txt
└── ustacks
    ├── CJ.alleles.tsv.gz
    ├── CJ.snps.tsv.gz
    ├── CJ.tags.tsv.gz
    ├── Pelatius.alleles.tsv.gz
    ├── Pelatius.snps.tsv.gz
    ├── Pelatius.tags.tsv.gz
    ├── Pkonojoi.alleles.tsv.gz
    ├── Pkonojoi.snps.tsv.gz
    ├── Pkonojoi.tags.tsv.gz
    └── id_sample.csv
```

## Usage : classification

通过指定 `database_accessor.py` 结果中的 `blast_results`，可以计算每个物种的 BLAST 结果，确定其归属于哪个物种。
本程序是分类的一个示例，因此建议在此基础上进一步添加其他分类方法。


```
usage: classification.py [-h] [-i INTSV] -o OUTDIR [-indir INDIR]

blast res classifier

options:
  -h, --help            show this help message and exit
  -i INTSV, --intsv INTSV
                        blast tsv
  -o OUTDIR, --output_dir OUTDIR
                        output directory
  -indir INDIR, --input_dir INDIR
                        input tsv directory
```

実行例
```
python classification.py -indir /home/data/outdir-coral-from-docker/blast_results/ -o /home/data/outdir-coral-from-docker/classified/

python ./classification_proposal/classification.py -indir /home/data/eDNA/result/blast_results/ -o /home/data/eDNA/classified/

```

### indir

```
ls /home/data/outdir-coral-from-docker/blast_results/
CJ-G1342.tsv  CJ-S2203.tsv           Pelatius-TNG002.tsv  Pkonojoi-TAK001.tsv
CJ-G1343.tsv  CJ-S2204.tsv           Pelatius-TNG005.tsv  Pkonojoi-TAK003.tsv
CJ-G1344.tsv  CJ-S2205.tsv           Pelatius-TNG006.tsv  Pkonojoi-TAK005.tsv
CJ-G1345.tsv  Pelatius-IOU005.tsv    Pkonojoi-IOU002.tsv  Pkonojoi-TAK006.tsv
CJ-G1346.tsv  Pelatius-IOU009-2.tsv  Pkonojoi-IOU003.tsv
CJ-S2201.tsv  Pelatius-IOU010.tsv    Pkonojoi-IOU004.tsv
CJ-S2202.tsv  Pelatius-S2201.tsv     Pkonojoi-IOU007.tsv
```
blast_results的文件夹。输入文件需要为[sample_name].tsv的形式。

### outdir

```
ls /home/data/outdir-coral-from-docker/classified/
CJ-G1342.json
CJ-G1342_summary.tsv
CJ-G1343.json
CJ-G1343_summary.tsv
CJ-G1344.json
CJ-G1344_summary.tsv
CJ-G1345.json
CJ-G1345_summary.tsv
CJ-G1346.json
CJ-G1346_summary.tsv
CJ-S2201.json
```
结果以 [sample_name]_summary.tsv 和 [sample_name].json 的形式输出。

### output files

- count : 匹配到该物种的reads数量（针对每个reads，提取BLAST得分最高的结果。如果该reads也匹配到其他物种的最高得分，则按等比例分配计算）
- ratio : 匹配到该物种的比例
- query_contig_num : 作为查询投入的contig数量。从输入的fastq等文件中通过ustacks提取并确定的contig数量。

#### tsv
```
#sample count   ratio
Pelatius__3     1699.6842712842697      0.7345221569940664
Pkonojoi__2     552.9871572871572       0.23897457099704286
CJ__1   61.32857142857144       0.026503272008889994
```

#### json
```
{
    "query_contig_num": 2314,
    "data": {
        "Pkonojoi__2": {
            "count": 552.9871572871572,
            "ratio": 0.23897457099704286
        },
        "CJ__1": {
            "count": 61.32857142857144,
            "ratio": 0.026503272008889994
        },
        "Pelatius__3": {
            "count": 1699.6842712842697,
            "ratio": 0.7345221569940664
        }
    }
}
```