import pandas as pd
import os
import logging
from BCBio import GFF


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class GFFProcessor:
    def __init__(self, config):
        self.config = config
        self.input_total_faa = config.get('input_total_faa')
        self.output_dir = config.get('output_dir')
        self.cazyme_overview = config.get('cazyme_overview')
        self.non_cazyme_faa = config.get('non_cazyme_faa')
        self.cgc_sig_file = config.get('cgc_sig_file')
        self.input_gff = config.get('input_gff')
        self.output_gff = config.get('output_gff')
        self.gff_type = config.get('gff_type')


    def load_cgc_type(self):
        df = pd.read_csv(self.cazyme_overview, sep='\t', usecols=[0, 5, 6], names=['protein_id', '#ofTools', 'CAZyme'])
        df['#ofTools'] = pd.to_numeric(df['#ofTools'], errors='coerce')
        overview_df = df[df['#ofTools'] >= 2]
        overview_df['CGC_annotation'] = 'CAZyme' + '|' + overview_df['CAZyme'].astype(str)


        cgc_sig_df = pd.read_csv(self.cgc_sig_file, sep='\t', usecols=[0, 2, 10], header=None, names=['function_annotation', 'protein_id', 'type'])
        cgc_sig_df['CGC_annotation'] = cgc_sig_df['type'] + '|' + cgc_sig_df['function_annotation']

        combined_df = pd.concat([overview_df[['protein_id', 'CGC_annotation']], cgc_sig_df[['protein_id', 'CGC_annotation']]], ignore_index=True)
        combined_df = combined_df.groupby('protein_id')['CGC_annotation'].apply(lambda x: '+'.join(set(x))).reset_index()
        #print(combined_df)
        return combined_df.set_index('protein_id').to_dict('index')
        
    def process_gff(self):
        cgc_data = self.load_cgc_type()
        with open(self.input_gff) as input_file, open(self.output_gff, 'w') as output_file:
            for record in GFF.parse(input_file):
                for feature in record.features:
                    if feature.type == 'gene':
                        protein_id = 'unknown'
                        cgc_annotation = 'unknown'
                        if self.gff_type == "NCBI":
                            non_mRNA_found = False
                            for sub_feature in feature.sub_features:
                                if 'mRNA' not in sub_feature.type:
                                    protein_id = 'NA'
                                    cgc_annotation = 'Other' + '|' + sub_feature.type
                                    non_mRNA_found = True
                                    break

                            if non_mRNA_found:
                                start, end, strand = feature.location.start + 1, feature.location.end, '+' if feature.location.strand >= 0 else '-'
                                line = f"{record.id}\t.\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\tprotein_id={protein_id};CGC_annotation={cgc_annotation}\n"
                                output_file.write(line)
                                continue

                            for sub_feature in feature.sub_features:
                                if sub_feature.type == 'mRNA':
                                    for sub_sub_feature in sub_feature.sub_features:
                                        if sub_sub_feature.type == 'CDS':
                                            protein_id = sub_sub_feature.qualifiers.get('protein_id', ['unknown'])[0]
                                            break
                                    if protein_id != 'unknown':
                                        break

                        elif self.gff_type == "JGI":
                            protein_id = feature.qualifiers.get("proteinId", ["unknown"])[0]
                        elif self.gff_type == "prodigal":
                            protein_id = feature.qualifiers.get("ID", ["unknown"])[0]

                        cgc_annotation = cgc_data.get(protein_id, {}).get('CGC_annotation', 'Other|null')

                        start, end, strand = feature.location.start + 1, feature.location.end, '+' if feature.location.strand >= 0 else '-'
                        line = f"{record.id}\t.\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\tprotein_id={protein_id};CGC_annotation={cgc_annotation}\n"
                        output_file.write(line)

        logging.info(f"Updated GFF file saved to {self.output_gff}.")







if __name__ == '__main__':
    config = {
    'input_total_faa': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output/uniInput.faa',
    'output_dir': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output',
    'cazyme_overview': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output/overview.tsv',
    'cgc_sig_file': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output/test_total_cgc_info.tsv',
    'input_gff': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output/uniInput.gff',
    'output_gff': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/NCBI_test/output/cgc.gff',
    'gff_type': 'NCBI'  # ncbi, jgi, prodigal
}
processor = GFFProcessor(config)
processor.process_gff()