import pandas as pd


class DBCANProcessor:
    def __init__(self, config):
        self.config = config
        self.output_file = config.get('output_file')
        self.mapping_file = config.get('mapping_file')

    def load_substrate_mapping(self):
        """
        example of mapping file:
        Substrate_high_level	 Substrate_curated	Family	Name	        EC_Number
        lignin	                 lignin	            AA1	    ferroxidase	    1.10.3.2
        chitin                    chitin            CBM14   Long description    NA
        """
        try:
            df = pd.read_csv(self.mapping_file, sep='\t', header=None, skiprows=1, usecols=[2, 4, 0])
            df[4] = df[4].str.strip().fillna('-')
            df['key'] = df.apply(lambda x: (x[2], x[4]) , axis=1)
            return pd.Series(df[0].values, index=pd.MultiIndex.from_tuples(df['key'])).to_dict()
        except FileNotFoundError:
            print(f"can't find file: {self.mapping_file}")
            return {}

    def process_dbcan_sub(self):

        subs_dict = self.load_substrate_mapping()

        try:
            df = pd.read_csv(self.output_file, sep='\t')
            #  extract information from HMM Name : PL25_e0.hmm|PL25:38|PL0:1|3.2.1.122:13
            df['Subfam Name'] = df['HMM Name'].apply(lambda x: '|'.join(p.split('.')[0] for p in x.split('|') if '.hmm' in p))
            df['Subfam Composition'] = df['HMM Name'].apply(lambda x: '|'.join(p for p in x.split('|') if '.hmm' not in p and len(p.split('.')) != 4))
            df['Subfam EC'] = df['HMM Name'].apply(lambda x: '|'.join(p for p in x.split('|') if len(p.split('.')) == 4))
            df['Substrate'] = df['HMM Name'].apply(lambda x: self.get_substrates(x, subs_dict))
            df.drop('HMM Name', axis=1, inplace=True)
            columns = ['Subfam Name', 'Subfam Composition', 'Subfam EC', 'Substrate', 'HMM Length', 'Target Name', 'Target Length', 'i-Evalue', 'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name']
            df = df[columns]
            df.to_csv(self.output_file, sep='\t', index=False)
            print("suceessfully processed dbcan substrate file")
        except Exception as e:
            print(f"error: {e}")

    def get_substrates(self, profile_info, subs_dict):

        parts = profile_info.split('|')
        substrates = set()
        key1 = parts[0].split('.hmm')[0].split("_")[0]
        if key1.startswith('CBM'):
            key2 = '-'
            if (key1, key2) in subs_dict:
                substrates.add(subs_dict[(key1, key2)])

        else:
            for p in parts:
                if ':' in p and '.' in p.split(':')[0] and len(p.split(':')[0].split('.')) == 4:
                    key2 = p.split(':')[0]
                    if (key1, key2) in subs_dict:
                        substrates.add(subs_dict[(key1, key2)])
                        
                if not substrates and not key1.startswith('CBM'):
                    if (key1, '-') in subs_dict:
                        substrates.add(subs_dict[(key1, '-')])

        return ', '.join(substrates) if substrates else '-'


# if __name__ == "__main__":
#         config = {
#         'input_file': 'dbCAN_hmm_results.tsv',
#         'output_file': 'dbCAN-sub.substrate.tsv',
#         'mapping_file': '/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/part1_script/dbCAN_db/fam-substrate-mapping.tsv'
#     }
#         processor = DBCANProcessor(config)

#         processor.process_dbcan_sub()
