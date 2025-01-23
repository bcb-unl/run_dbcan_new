import pyhmmer
import logging
import psutil
import pandas as pd
import os


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class PyHMMERProcessor:
    def __init__(self, config):
        self.config = config
        self.hmm_file = config.get('hmm_file')
        self.input_faa = config.get('input_faa')
        self.e_value_threshold = config.get('e_value_threshold', 1e-15)
        self.coverage_threshold = config.get('coverage_threshold', 0.35)
        self.output_file = config.get('output_file')
        self.hmmer_cpu = config.get('hmmer_cpu')

    def hmmsearch(self):
        available_memory = psutil.virtual_memory().available
        target_size = os.stat(self.input_faa).st_size
        hmm_files = pyhmmer.plan7.HMMFile(self.hmm_file)
        results = []

        with pyhmmer.easel.SequenceFile(self.input_faa, digital=True) as seqs:
            targets = seqs.read_block() if target_size < available_memory * 0.1 else seqs
            for hits in pyhmmer.hmmsearch(hmm_files, targets, cpus=self.hmmer_cpu, domE=1e-15):
                for hit in hits:
                    for domain in hit.domains.included:
                        coverage = (domain.alignment.hmm_to - domain.alignment.hmm_from + 1) / domain.alignment.hmm_length
                        hmm_name = domain.alignment.hmm_name.decode('utf-8')
                        if "GT2_" in hmm_name:
                            hmm_name = "GT2.hmm"
                        hmm_length = domain.alignment.hmm_length
                        target_name = domain.alignment.target_name.decode('utf-8')
                        target_length = domain.alignment.target_length
                        i_evalue = domain.i_evalue
                        hmm_from = domain.alignment.hmm_from
                        hmm_to = domain.alignment.hmm_to
                        target_from = domain.alignment.target_from
                        target_to = domain.alignment.target_to
                        hmm_file_name=hmm_files.name.split("/")[-1].split(".")[0]
                        if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                            results.append([hmm_name, hmm_length, target_name, target_length, i_evalue, hmm_from, hmm_to, target_from, target_to, coverage, hmm_file_name])
                            
        logging.info(f"{self.hmm_file} PyHMMER search completed. Found {len(results)} hits.")
        if results:
            df = pd.DataFrame(results, columns=[
                'HMM Name', 'HMM Length', 'Target Name', 'Target Length', 'i-Evalue',
                'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name'])
            df.sort_values(by=['Target Name', 'Target From', 'Target To'], inplace=True)
            df_filtered = self.filter_overlaps(df)
            df_filtered.to_csv(self.output_file, index=False, sep='\t')
        else:
            df = pd.DataFrame(columns=[
                'HMM Name', 'HMM Length', 'Target Name', 'Target Length', 'i-Evalue',
                'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name'])
            df.to_csv(self.output_file, index=False, sep='\t')
    def filter_overlaps(self, df):
        filtered = []
        grouped = df.groupby('Target Name')

        for name, group in grouped:
            group = group.reset_index(drop=True)
            keep = []

            for i in range(len(group)):
                if not keep:
                    keep.append(group.iloc[i])
                    continue

                last = keep[-1]
                current = group.iloc[i]
                overlap = min(last['Target To'], current['Target To']) - max(last['Target From'], current['Target From'])
                if overlap > 0:
                    overlap_ratio_last = overlap / (last['Target To'] - last['Target From'])
                    overlap_ratio_current = overlap / (current['Target To'] - current['Target From'])

                    if overlap_ratio_last > 0.5 or overlap_ratio_current > 0.5:
                        if last['i-Evalue'] > current['i-Evalue']:
                            keep[-1] = current
                    else:
                        keep.append(current)
                else:
                    keep.append(current)

            filtered.extend(keep)
        
        return pd.DataFrame(filtered)

# Example usage
# if __name__ == "__main__":
#     config = {
#         'hmm_file': 'dbCAN_databases/dbCAN_sub.hmm',
#         'input_faa': 'test.faa',
#         'e_value_threshold': 1e-15
#     }
#     processor = PyHMMERProcessor(config)
#     processor.hmmsearch()
