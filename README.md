# CoralDatabaseBuilder
Migseq用に撮られたデータを用いてblastデータベースの構築および検索を行うプログラム。

## docker
```
# build image
docker image build -t classifier .
# run container
docker run -itd -v $PWD:/home/data --name coral_classifier classifier bash
# exec container
docker exec -it coral_classifier bash
```

## Usage : database_accessor
fastq => migseq前処理 => ustacks => fasta抽出 => blast検索 or 登録
```
usage: database_accessor.py [-h] -c CONF

database accessor

options:
  -h, --help            show this help message and exit
  -c CONF, --conf CONF  config file
```
利用例
```
python database_accessor.py -c config.conf
```

## config file

```
[tools]
cutadapt=cutadapt
fastq_quality_filter=fastq_quality_filter
fastx_trimmer=fastx_trimmer
fastq_to_fasta=fastq_to_fasta
bgzip=bgzip
ustacks=ustacks

makeblastdb=makeblastdb
blastn=blastn

[settings]
# choose mode from makedb and search
mode=makedb
input_dir=/home/data/indir-coral/ref/
output_dir=/home/data/outdir-coral-mkdb-from-docker/
database=/home/data/testdb/coral_from_docker

threads=2


quality_control=T
quality=30
reverse_cut_length=15
forward_cut_length=0

adapter_trim=T
error_rate=0.05
fra=GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
rra=CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC
min_len = 80

fastq_to_fasta=T
from_fa = F
stacks=T
ustacks_opts = -m 3 -M 2
extract_consensus_fasta = T
from_cons_fa = F
```

- [tools] : プログラム内で利用しているツールのパスを指定
  - cutadapt=cutadapt : cutadaptの実行可能なコマンドを指定
  - fastq_quality_filter=fastq_quality_filter : fastq_quality_filterの実行可能なコマンドを指定
  - fastx_trimmer=fastx_trimmer : fastx_trimmerの実行可能なコマンドを指定
  - fastq_to_fasta=fastq_to_fasta : fastq_to_fastaの実行可能なコマンドを指定
  - bgzip=bgzip : bgzipの実行可能なコマンドを指定
  - ustacks=ustacks : ustacksの実行可能なコマンドを指定
  - makeblastdb=makeblastdb : makeblastdbの実行可能なコマンドを指定
  - blastn=blastn : blastnの実行可能なコマンドを指定

- [settings] : プログラム実行の設定
  - mode=makedb : データベース作成の場合makedb、検索の場合はsearchを指定
  - input_dir=/home/data/indir-coral/ref/ : 入力ファイルがあるディレクトリを指定
  - output_dir=/home/data/outdir-coral-mkdb-from-docker/ : 出力ファイルを格納するディレクトリを指定（makedbではdb作成時に生成される中間ファイル群を保存する）
  - database=/home/data/testdb/coral_from_docker : databaseファイルへのパスを指定（makedbの場合にはここで指定したファイルを作成。これに拡張子がついてファイルが生成されるため、最後に/などをつけてはいけない。）
  - threads=2 : 利用するスレッド数
  - quality_control=T : quality controlを実行する(T) or しない(F)
  - quality=30 : quality control時のquality値
  - reverse_cut_length=15 : fastx_trimmerに渡すreverse readを切り出す値(15の場合14塩基切り取り)
  - forward_cut_length=0 : fastx_trimmerに渡すforward readを切り出す値(15の場合14塩基切り取り)
  - adapter_trim=T : cutadaptによるadapter trimを実行する(T) or しない(F)
  - error_rate=0.05 : adapterのエラー許容率
  - fra=GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC : forward readから除去するadapter
  - rra=CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC : reverse readから除去するadapter
  - min_len = 80 : 配列の最小長
  - fastq_to_fasta=T : fastqからfastaへの変換を行う(T) or 行わない(F)
  - from_fa = F : Fastaからスタートする(T) or しない(F)
  - stacks=T : stacksによる解析を行う(T) or 行わない(F)
  - ustacks_opts = -m 3 -M 2 : ustacksのオプション
  - extract_consensus_fasta = T : ustacksの結果からconsensusと判定された配列を抜き出したfastaを生成する(T) or しない(F)
  - from_cons_fa = F : ustacksの結果から生成したfastaからスタートする(T) or しない(F)


### input_dir

入力ファイルはpaired endでforward readは``[sample_name]_L001_R1_001.fastq.gz``、reverse readは``[sample_name]_L001_R2_001.fastq.gz``の命名規則を満たしている必要がある。

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

- adapter_trim : fastx_trimmerでの切り取り、およびcutadaptでadapterで切り取った後のfastqを格納するディレクトリ
- quality_control : cutadaptでquality controlを行った結果のfastqを格納するディレクトリ
- fasta : fastqをfastaに変換したものを格納するディレクトリ
- ustacks : fastaをustacksで処理した結果を格納するディレクトリ
- consensus_fasta : ustacksが抽出したconsensusと判断されたリードをfastaとして保存するディレクトリ
- blast_results : blastの結果を格納するディレクトリ（searchの場合のみ。outfmt 6でtsvを出力）

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

database_accessor.pyの結果中のblast_resultsを指定して実行することで、それぞれの種についてblastの結果からどの種に分類されるか計算する。
本プログラムは分類の１例であるので、ここからさらに別の分類方法も追加していくことが望まれる。


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
blast_resultsのディレクトリ。[sample_name].tsvという形になっている。

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
結果が[sample_name]_summary.tsvと[sample_name].jsonの形式で出力されている。

### output files

- count : その種に当たったリード数(各リードについて、blast scoreがtopのものを抽出。他の種にもtopで当たっている場合には当分配で割り算)
- ratio : その種に当たっている割合
- query_contig_num : クエリとして投入されたcontig数。入力したfastqなどからustacksで判定されたコンティ具を抽出して算出。

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
