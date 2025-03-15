import pandas as pd
import os
import re

class OverviewGenerator:
    def __init__(self, config):
        self.config = config
        self.output_dir = config.get('output_dir')

    def load_data(self):
        file_paths = {
            'diamond': os.path.join(self.output_dir, 'diamond_results.tsv'),
            'dbcan_sub': os.path.join(self.output_dir, 'dbCAN-sub.substrate.tsv'),
            'dbcan-hmm': os.path.join(self.output_dir, 'dbCAN_hmm_results.tsv')
        }
        
        
        data = {}
        for key, file_path in file_paths.items():
            if os.path.exists(file_path):
                columns = {
                    'diamond': ['Gene ID', 'CAZy ID'],
                    'dbcan_sub': ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue'],
                    'dbcan-hmm': ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']
                }
                df = pd.read_csv(file_path, sep='\t', usecols=columns[key])
                if key == 'diamond':
                    def extract_cazy_id(cazy_id):
                        parts = cazy_id.split('|')
                        for part in parts:
                            if re.match(r"^(GH|GT|CBM|AA|CE|PL)", part):
                                return '+'.join(parts[parts.index(part):])
                        return cazy_id
                    df['CAZy ID'] = df['CAZy ID'].apply(extract_cazy_id)
                else:
                    df[columns[key][1]] = df[columns[key][1]].apply(lambda x: x.split('.hmm')[0] if isinstance(x, str) and '.hmm' in x else x)
                data[key] = df
        return data

    def calculate_overlap(self, start1, end1, start2, end2):
        start_max = max(start1, start2)
        end_min = min(end1, end2)
        overlap = max(0, end_min - start_max + 1)
        length1 = end1 - start1 + 1
        length2 = end2 - start2 + 1
        return overlap / min(length1, length2) > 0.5

    def determine_best_result(self, gene_id, data):
        results = {'EC#': '-', 'dbCAN_hmm': '-', 'dbCAN_sub': '-', 'DIAMOND': '-', '#ofTools': 0, 'Recommend Results': '-'}

        if 'dbcan-hmm' in data and not data['dbcan-hmm'].empty:
            hmm_results = data['dbcan-hmm'][data['dbcan-hmm']['Target Name'] == gene_id]
            if not hmm_results.empty:
                results['dbCAN_hmm'] = '+'.join([f"{row['HMM Name']}({row['Target From']}-{row['Target To']})" for index, row in hmm_results.iterrows()])
                results['#ofTools'] += 1

        if 'dbcan_sub' in data and not data['dbcan_sub'].empty:
            sub_results = data['dbcan_sub'][data['dbcan_sub']['Target Name'] == gene_id]
            if not sub_results.empty:
                results['dbCAN_sub'] = '+'.join([f"{row['Subfam Name']}({row['Target From']}-{row['Target To']})" for index, row in sub_results.iterrows()])
                results['EC#'] = '|'.join([str(ec) if ec is not None else '-' for ec in sub_results['Subfam EC'].fillna('-').tolist()])
                results['#ofTools'] += 1

        if 'diamond' in data and not data['diamond'].empty:
            diamond_results = data['diamond'][data['diamond']['Gene ID'] == gene_id]
            if not diamond_results.empty:
                results['DIAMOND'] = '+'.join(diamond_results['CAZy ID'].tolist())
                results['#ofTools'] += 1

        if results['dbCAN_hmm'] != '-' and results['dbCAN_sub'] != '-':
            overlap_results = []
            for _, sr in sub_results.iterrows():
                sub_overlap = False
                for _, hr in hmm_results.iterrows():
                    if self.calculate_overlap(sr['Target From'], sr['Target To'], hr['Target From'], hr['Target To']):
                        if "_" in hr['HMM Name'] or sr['i-Evalue'] > hr['i-Evalue']:
                            overlap_results.append((hr['HMM Name'], hr['Target From']))
                        else:
                            overlap_results.append((sr['Subfam Name'], sr['Target From']))
                        sub_overlap = True
                if not sub_overlap:
                    overlap_results.append((sr['Subfam Name'], sr['Target From']))
            for _, hr in hmm_results.iterrows():
                if all(not self.calculate_overlap(sr['Target From'], sr['Target To'], hr['Target From'], hr['Target To']) for _, sr in sub_results.iterrows()):
                    overlap_results.append((hr['HMM Name'], hr['Target From']))

            sorted_results = sorted(overlap_results, key=lambda x: x[1])
            results['Recommend Results'] = '|'.join([str(res[0]) for res in sorted_results])

        elif results['dbCAN_hmm'] != '-':
            results['Recommend Results'] = '|'.join([name.split('(')[0] for name in results['dbCAN_hmm'].split('|')])
        elif results['dbCAN_sub'] != '-':
            results['Recommend Results'] = '|'.join([name.split('(')[0] for name in results['dbCAN_sub'].split('|')])

        return results

    def aggregate_data(self, gene_ids, data):
        aggregated_results = []
        for gene_id in sorted(gene_ids):  #sort with gene id
            result = self.determine_best_result(gene_id, data)
            aggregated_results.append([gene_id] + list(result.values()))
        return pd.DataFrame(aggregated_results, columns=['Gene ID', 'EC#', 'HMMER', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results'])

    def run(self):
        loaded_data = self.load_data()
        gene_ids = set()
        for dataset in loaded_data.values():
            gene_ids.update(dataset['Target Name'].unique() if 'Target Name' in dataset.columns else dataset['Gene ID'].unique())
        aggregated_results = self.aggregate_data(gene_ids, loaded_data)
        output_path = os.path.join(self.output_dir, 'overview.tsv')
        aggregated_results.to_csv(output_path, sep='\t', index=False)
        print("Aggregated results saved to:", output_path)


# if __name__ == "__main__":
#     generator = OverviewGenerator('/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/test')
#     loaded_data = generator.load_data()
#     gene_ids = set()
#     for dataset in loaded_data.values():
#         gene_ids.update(dataset['Target Name'].unique() if 'Target Name' in dataset.columns else dataset['Gene ID'].unique())
#     aggregated_results = generator.aggregate_data(gene_ids, loaded_data)
#     aggregated_results.to_csv(os.path.join(generator.out_path, 'aggregated_results.tsv'), sep='\t', index=False)
#     print("Aggregated results saved.")
