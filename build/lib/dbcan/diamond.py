import subprocess
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class DiamondProcessor:
    def __init__(self, config):
        self.config = config
        self.diamond_db = config.get('diamond_db')
        self.input_faa = config.get('input_faa')
        self.output_file = config.get('output_file')
        self.e_value_threshold = config.get('e_value_threshold')
        self.threads = config.get('threads')
        self.verbose_option = config.get('verbose_option', False)
        self.coverage_threshold_tc = config.get('coverage_threshold_tc')
        if self.verbose_option:
            self.verbose_option = '-v'
        else:
            self.verbose_option = '--quiet'


    def run_diamond(self):
        """
        Run DIAMOND BLASTp with the provided configuration.
        """
        cmd = [
            'diamond', 'blastp',
            '--db', self.diamond_db,
            '--query', self.input_faa,
            '--out', self.output_file,
            '--outfmt', '6',
            '--evalue', str(self.e_value_threshold),
            '--max-target-seqs', '1',
            '--threads', str(self.threads), #default use all cpus available
            self.verbose_option
        ]
        logging.info("Running DIAMOND BLASTp...")
        try:
            subprocess.run(cmd, check=True)
            logging.info("DIAMOND BLASTp completed.")
        except subprocess.CalledProcessError as e:
            logging.error(f"DIAMOND BLASTp failed :{e}")

    def format_results(self):
        """
        reformat DIAMOND output and filter based on e-value threshold.
        """
        filtered_df = pd.read_csv(self.output_file, sep='\t', header=None, names=[
            'Gene ID', 
            'CAZy ID', 
            '% Identical', 
            'Length', 
            'Mismatches', 
            'Gap Open',
            'Gene Start', 
            'Gene End', 
            'CAZy Start', 
            'CAZy End', 
            'E Value', 
            'Bit Score'
        ])
        filtered_df.to_csv(self.output_file, sep='\t', index=False)

    def run_tcdb_diamond(self):
        """
        Run DIAMOND BLASTp with the provided configuration.
        """
        cmd = [
            'diamond', 'blastp',
            '--db', self.diamond_db,
            '--query', self.input_faa,
            '--out', self.output_file,
            '--outfmt', '6', 'sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp',
            '--evalue', str(self.e_value_threshold),
            '--max-target-seqs', '1',
            '--query-cover', str(self.coverage_threshold_tc),
            '--threads', str(self.threads), #default use all cpus available
            self.verbose_option
        ]
        logging.info("Running DIAMOND BLASTp TCDB...")
        try:
            subprocess.run(cmd, check=True)
            logging.info("DIAMOND BLASTp TCDB completed.")
        except subprocess.CalledProcessError as e:
            logging.error(f"DIAMOND BLASTp TCDB failed :{e}")

    def format_results_tcdb(self):
        """
        reformat DIAMOND output and filter based on e-value threshold.
        """
        filtered_df = pd.read_csv(self.output_file, sep='\t', header=None, names=[
            'TCDB ID', 
            'TCDB Length', 
            'Target ID', 
            'Target Length', 
            'EVALUE', 
            'TCDB START',
            'TCDB END', 
            'QSTART', 
            'QEND', 
            'COVERAGE'
        ])
        filtered_df['TCDB ID'] = filtered_df['TCDB ID'].apply(lambda x: x.split(' ')[0].split('|')[-1])
        filtered_df['Database'] = 'TC'
        filtered_df.to_csv(self.output_file, sep='\t', index=False)

# # Example usage
# if __name__ == "__main__":
#     config = {
#         'diamond_db': './dbCAN_databases/CAZy.dmnd',
#         'input_faa': 'test.faa',
#         'output_file': './dbCAN_output/diamond_results.tsv',
#         'output_dir': 'dbCAN_output',
#         'e_value_threshold': 1e-102,
#         'threads': os.cpu_count(),
#         'verbose_option': False
#     }
#     processor = DiamondProcessor(config)
#     processor.run_diamond()
#     filtered_results = processor.format_results()
