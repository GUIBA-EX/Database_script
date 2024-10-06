import argparse
import glob
import json
import logging
import os
from ACF import analysis_complete_flag as acf

class ClassificationProposal:
    def __init__(self, intsv, outdir, indir):
        self.intsv = intsv
        self.outdir = outdir
        self.indir = indir
        os.makedirs(self.outdir, exist_ok=True)
        self.setup_logger()

    def setup_logger(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(asctime)s] %(levelname)s : %(message)s', '%Y-%m-%d %H:%M:%S')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def classify_method1(self):
        ref_counts = {}
        query_contigs = set()
        
        try:
            with open(self.intsv, 'r') as intsv:
                for line in intsv:
                    data = line.strip().split('\t')
                    if len(data) != 12:
                        continue
                    
                    query, ref = data[0], data[1]
                    query_contigs.add(query)
                    ref_counts[ref] = ref_counts.get(ref, 0) + 1
            
            total_matches = sum(ref_counts.values())
            ref_score_dict = {ref: {'count': count, 'ratio': count / total_matches} 
                              for ref, count in ref_counts.items()}
            
            sample_name = os.path.basename(self.intsv)[:-4]
            
            # 输出JSON
            with open(os.path.join(self.outdir, f'{sample_name}.json'), 'w') as outjp:
                json.dump({'query_contig_num': len(query_contigs), 'data': ref_score_dict}, outjp, indent=4)
            
            # 输出摘要TSV
            outobj = sorted([(ref, data['count'], data['ratio']) for ref, data in ref_score_dict.items()], 
                            key=lambda x: x[2], reverse=True)
            outobj = ['#sample\tcount\tratio'] + [f'{i[0]}\t{i[1]}\t{i[2]}' for i in outobj]
            with open(os.path.join(self.outdir, f'{sample_name}_summary.tsv'), 'w') as outtp:
                outtp.write('\n'.join(outobj))
            
            self.logger.info(f"处理了 {len(query_contigs)} 个查询序列和 {len(ref_counts)} 个参考序列")
            self.logger.info(f"分析完成,结果保存在 {self.outdir}")
        
        except IOError as e:
            self.logger.error(f"无法读取或写入文件: {str(e)}")

    def run(self):
        if self.indir is None or self.indir == 'None':
            if not self.intsv:
                self.logger.error("未提供输入文件或目录")
                return
            self.classify_method1()
        else:
            tsvs = glob.glob(os.path.join(self.indir, '*.tsv'))
            self.logger.info(f"在目录 {self.indir} 中找到 {len(tsvs)} 个TSV文件")
            for tsv in tsvs:
                self.logger.info(f"处理文件: {tsv}")
                self.intsv = tsv
                self.classify_method1()
        self.logger.info("所有文件处理完成")

def main():
    parser = argparse.ArgumentParser(description="blast res classifier")
    parser.add_argument('-i', '--intsv', type=str, help='blast tsv')
    parser.add_argument('-o', '--output_dir', type=str, help='output directory', required=True)
    parser.add_argument('-indir', '--input_dir', type=str, help='input tsv directory')
    args = parser.parse_args()

    CP = ClassificationProposal(intsv=args.intsv, outdir=args.output_dir, indir=args.input_dir)
    CP.run()
    acf(outdir=args.output_dir)

if __name__ == '__main__':
    main()