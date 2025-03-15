import pandas as pd
import os
import re
import logging
from parameter import OverviewGeneratorConfig
from Bio import SeqIO

class OverviewGenerator:
    """Generate overview of CAZyme annotations using template method pattern"""
    
    def __init__(self, config: OverviewGeneratorConfig):
        """Initialize with configuration"""
        self.config = config
        self._setup_processor()
    
    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        self.output_dir = self.config.output_dir
        self.file_paths = self._derive_file_paths()
        self.column_names = self._derive_column_names()
        self.overview_columns = ['Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results']
        self.overlap_threshold = 0.5
        
        # For non-CAZyme FAA generation
        self.cazyme_overview = os.path.join(self.output_dir, 'overview.tsv')
        self.input_total_faa = self._derive_input_total_faa()
        
        # Validate required attributes
        self._validate_attributes()
    
    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
        
        # Check that at least one result file exists
        found_files = False
        for key, file_path in self.file_paths.items():
            if os.path.exists(file_path):
                found_files = True
                logging.info(f"Found {key} results at {file_path}")
            else:
                logging.warning(f"{key} results not found at {file_path}")
        
        if not found_files:
            logging.warning("No CAZyme annotation results found. Overview will be empty.")
    
    def _derive_file_paths(self):
        """Derive file paths for annotation results"""
        return {
            'diamond': os.path.join(self.output_dir, 'diamond.out'),
            'dbcan_sub': os.path.join(self.output_dir, 'dbCANsub_hmm_results.tsv'),
            'dbcan_hmm': os.path.join(self.output_dir, 'dbCAN_hmm_results.tsv')
        }
    
    def _derive_column_names(self):
        """Derive column names for annotation results"""
        return {
            'diamond': ['Gene ID', 'CAZy ID'],
            'dbcan_sub': ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue'],
            'dbcan_hmm': ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']
        }
    
    def _derive_input_total_faa(self):
        """Derive path to total input FAA file"""
        return os.path.join(self.output_dir, 'uniInput.faa')

    def load_data(self):
        """Load data from annotation result files"""
        data = {}
        for key, file_path in self.file_paths.items():
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep='\t')
                    
                    # Ensure all required columns exist
                    if not all(col in df.columns for col in self.column_names[key]):
                        logging.warning(f"Missing columns in {file_path}. Expected: {self.column_names[key]}, Found: {df.columns}")
                        continue
                    
                    # Filter to only required columns
                    df = df[self.column_names[key]]
                    
                    # Process specific to each result type
                    if key == 'diamond':
                        df['CAZy ID'] = df['CAZy ID'].apply(self.extract_cazy_id)
                    elif key in ['dbcan_hmm', 'dbcan_sub']:
                        # Handle HMM name formatting
                        hmm_col = 'HMM Name' if key == 'dbcan_hmm' else 'Subfam Name'
                        df[hmm_col] = df[hmm_col].apply(
                            lambda x: x.split('.hmm')[0] if isinstance(x, str) and '.hmm' in x else x
                        )
                    
                    data[key] = df
                    logging.info(f"Loaded {len(df)} rows from {key} results")
                except Exception as e:
                    logging.error(f"Error loading {key} results: {e}")
        
        return data

    @staticmethod
    def extract_cazy_id(cazy_id):
        """Extract CAZy ID from DIAMOND result"""
        if not isinstance(cazy_id, str):
            return cazy_id
        
        parts = cazy_id.split('|')
        for part in parts:
            if re.match(r"^(GH|GT|CBM|AA|CE|PL)", part):
                return '+'.join(parts[parts.index(part):])
        return cazy_id

    def calculate_overlap(self, start1, end1, start2, end2):
        """Calculate overlap between two regions"""
        start_max = max(start1, start2)
        end_min = min(end1, end2)
        overlap = max(0, end_min - start_max + 1)
        length1 = end1 - start1 + 1
        length2 = end2 - start2 + 1
        return overlap / min(length1, length2) > self.overlap_threshold

    def determine_best_result(self, gene_id, data):
        """Determine best result for a gene"""
        results = {'EC#': '-', 'dbCAN_hmm': '-', 'dbCAN_sub': '-', 'DIAMOND': '-', '#ofTools': 0, 'Recommend Results': '-'}

        # Process HMMER results
        if 'dbcan_hmm' in data and not data['dbcan_hmm'].empty:
            hmm_results = data['dbcan_hmm'][data['dbcan_hmm']['Target Name'] == gene_id]
            if not hmm_results.empty:
                results['dbCAN_hmm'] = '+'.join([f"{row['HMM Name']}({row['Target From']}-{row['Target To']})" for _, row in hmm_results.iterrows()])
                results['#ofTools'] += 1

        # Process dbCAN-sub results
        if 'dbcan_sub' in data and not data['dbcan_sub'].empty:
            sub_results = data['dbcan_sub'][data['dbcan_sub']['Target Name'] == gene_id]
            if not sub_results.empty:
                results['dbCAN_sub'] = '+'.join([f"{row['Subfam Name']}({row['Target From']}-{row['Target To']})" for _, row in sub_results.iterrows()])
                results['EC#'] = '|'.join([str(ec) if ec is not None else '-' for ec in sub_results['Subfam EC'].fillna('-').tolist()])
                results['#ofTools'] += 1

        # Process DIAMOND results
        if 'diamond' in data and not data['diamond'].empty:
            diamond_results = data['diamond'][data['diamond']['Gene ID'] == gene_id]
            if not diamond_results.empty:
                results['DIAMOND'] = '+'.join(diamond_results['CAZy ID'].tolist())
                results['#ofTools'] += 1

        # Only add Recommend Results if at least 2 tools detected the gene
        if results['#ofTools'] >= 2:
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
                results['Recommend Results'] = '|'.join([name.split('(')[0] for name in results['dbCAN_hmm'].split('+')])
            elif results['dbCAN_sub'] != '-':
                results['Recommend Results'] = '|'.join([name.split('(')[0] for name in results['dbCAN_sub'].split('+')])

        return results

    def aggregate_data(self, gene_ids, data):
        """Aggregate data for all genes"""
        aggregated_results = []
        for gene_id in sorted(gene_ids):
            result = self.determine_best_result(gene_id, data)
            aggregated_results.append([gene_id] + list(result.values()))
        return pd.DataFrame(aggregated_results, columns=self.overview_columns)

    def generate_non_cazyme_faa(self):
        """Generate FAA file with non-CAZyme sequences"""
        try:
            # Only generate if input and output can be determined
            if not hasattr(self, 'input_total_faa') or self.input_total_faa is None:
                logging.error("Cannot generate non-CAZyme FAA: input_total_faa not set")
                return
                
            if not os.path.exists(self.input_total_faa):
                logging.error(f"Cannot generate non-CAZyme FAA: input file not found: {self.input_total_faa}")
                return
                
            if not os.path.exists(self.cazyme_overview):
                logging.error(f"Cannot generate non-CAZyme FAA: overview file not found: {self.cazyme_overview}")
                return
                
            # Read overview and get CAZyme IDs with >= 2 tools
            df = pd.read_csv(self.cazyme_overview, sep='\t')
            filtered_df = df[df['#ofTools'] >= 2]
            cazyme_ids = set(filtered_df['Gene ID'].tolist())
            
            # Write non-CAZymes to output
            output_path = os.path.join(self.output_dir, 'non_CAZyme.faa')
            count = 0
            
            with open(self.input_total_faa, 'r') as infile, open(output_path, 'w') as outfile:
                for record in SeqIO.parse(infile, 'fasta'):
                    header_id = record.id.split()[0]
                    if header_id not in cazyme_ids:
                        SeqIO.write(record, outfile, 'fasta')
                        count += 1
                        
            logging.info(f"Non-CAZyme FAA file generated with {count} sequences at {output_path}")
            
        except Exception as e:
            logging.error(f"Failed to generate non-CAZyme FAA: {str(e)}")

    def run(self):
        """Run overview generation"""
        try:
            # Load data from result files
            loaded_data = self.load_data()
            
            # If no data was loaded, create empty overview
            if not loaded_data:
                logging.warning("No annotation results found. Creating empty overview.")
                empty_df = pd.DataFrame(columns=self.overview_columns)
                output_path = os.path.join(self.output_dir, 'overview.tsv')
                empty_df.to_csv(output_path, sep='\t', index=False)
                print(f"Empty overview saved to: {output_path}")
                return
            
            # Collect all gene IDs from all datasets
            gene_ids = set()
            for key, dataset in loaded_data.items():
                id_col = 'Target Name' if key in ['dbcan_hmm', 'dbcan_sub'] else 'Gene ID'
                if id_col in dataset.columns:
                    gene_ids.update(dataset[id_col].unique())
            
            # Aggregate data for all genes
            aggregated_results = self.aggregate_data(gene_ids, loaded_data)
            
            # Save overview to file
            output_path = os.path.join(self.output_dir, 'overview.tsv')
            aggregated_results.to_csv(output_path, sep='\t', index=False)
            print(f"Aggregated results saved to: {output_path}")
            
            # Generate non-CAZyme FAA file
            self.generate_non_cazyme_faa()
            
        except Exception as e:
            logging.error(f"Error generating overview: {str(e)}")
            import traceback
            traceback.print_exc()
