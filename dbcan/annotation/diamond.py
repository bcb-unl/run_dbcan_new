import os
import subprocess
import logging
import pandas as pd

from dbcan.parameter import DiamondConfig, DiamondTCConfig
from dbcan.constants import CAZY_COLUMN_NAMES, TCDB_COLUMN_NAMES

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class DiamondProcessor:
    """Base Diamond processor class using template method pattern"""
    
    def __init__(self, config):
        """Initialize with configuration"""
        self.config = config
        self._setup_processor()
    
    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        self.diamond_db = self._derive_diamond_db()
        self.input_faa = self._derive_input_faa()
        self.output_file = self._derive_output_file() 
        self.e_value_threshold = self._derive_e_value_threshold()
        self.threads = self._derive_threads()
        self.verbose_option = self._derive_verbose_option()
        
        # Validate required attributes
        self._validate_attributes()
    
    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        required_attrs = ['diamond_db', 'input_faa', 'output_file', 
                            'e_value_threshold', 'threads']
        
        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")
                
        # Also validate file existence
        if not os.path.exists(self.diamond_db):
            raise FileNotFoundError(f"Database file not found: {self.diamond_db}")
            
        if not os.path.exists(self.input_faa):
            raise FileNotFoundError(f"Input file not found: {self.input_faa}")
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(self.output_file)), exist_ok=True)
    
    def _derive_diamond_db(self):
        """Derive DIAMOND database path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_diamond_db()")
    
    def _derive_input_faa(self):
        """Derive input protein sequence file path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_input_faa()")
    
    def _derive_output_file(self):
        """Derive output file path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_output_file()")
    
    def _derive_e_value_threshold(self):
        """Derive E-value threshold - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_e_value_threshold()")
    
    def _derive_threads(self):
        """Derive number of threads to use"""
        return self.config.threads
        
    def _derive_verbose_option(self):
        """Derive verbose option flag"""
        return getattr(self.config, 'verbose_option', False)

    def run_diamond(self, outfmt='6', extra_args=None):
        """Run DIAMOND BLASTP"""
        cmd = [
            'diamond', 'blastp',
            '--db', self.diamond_db,
            '--query', self.input_faa,
            '--out', self.output_file,
            '--outfmt', outfmt,
            '--evalue', str(self.e_value_threshold),
            '--max-target-seqs', '1',
            '--threads', str(self.threads),
            '-v' if self.verbose_option else '--quiet'
        ]
        
        if extra_args:
            cmd.extend(extra_args)

        logging.info(f"Running DIAMOND BLASTp with {os.path.basename(self.diamond_db)}...")
        try:
            subprocess.run(cmd, check=True)
            logging.info("DIAMOND BLASTp completed")
        except subprocess.CalledProcessError as e:
            logging.error(f"DIAMOND BLASTp failed: {e}")
            raise

    def format_results(self, column_names, extra_processing=None):
        """Format results"""
        if not os.path.exists(self.output_file) or os.stat(self.output_file).st_size == 0:
            logging.warning(f"No results to format: {self.output_file} is empty or missing")
            return
            
        try:
            filtered_df = pd.read_csv(self.output_file, sep='\t', header=None, names=column_names)
            
            if extra_processing:
                extra_processing(filtered_df)
                
            filtered_df.to_csv(self.output_file, sep='\t', index=False)
            logging.info(f"Results formatted and saved to {self.output_file}")
        except Exception as e:
            logging.error(f"Error formatting results: {e}")
            raise


class CAZyDiamondProcessor(DiamondProcessor):
    """CAZyme DIAMOND processor"""
    
    def _derive_diamond_db(self):
        """Get CAZyme DIAMOND database path"""
        return os.path.join(self.config.db_dir, "CAZy.dmnd")
    
    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")
    
    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, "diamond.out")
    
    def _derive_e_value_threshold(self):
        """Get E-value threshold for CAZyme searches"""
        # Use the value from config or default to 1e-15
        return getattr(self.config, 'e_value_threshold', 1e-15)
    
    def run(self):
        """Run CAZyme DIAMOND search"""
        self.run_diamond()

    def format_results(self):
        """Format CAZyme DIAMOND results"""
        super().format_results(CAZY_COLUMN_NAMES)


class TCDBDiamondProcessor(DiamondProcessor):
    """TCDB DIAMOND processor"""
    
    def _derive_diamond_db(self):
        """Get TCDB DIAMOND database path"""
        return os.path.join(self.config.db_dir, "TCDB.dmnd")
    
    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")
    
    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, "diamond.out.tc")
    
    def _derive_e_value_threshold(self):
        """Get E-value threshold for TCDB searches"""
        return getattr(self.config, 'e_value_threshold_tc', 1e-5)
    
    def _derive_coverage_threshold(self):
        """Get coverage threshold for TCDB searches"""
        return getattr(self.config, 'coverage_threshold_tc', 0.4)
    
    def run(self):
        """Run TCDB DIAMOND search"""
        # Get coverage threshold
        coverage_threshold = self._derive_coverage_threshold()
        
        # Set additional parameters
        extra_args = [
            '--outfmt', '6', 'sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp',
            '--query-cover', str(coverage_threshold)
        ]
        
        self.run_diamond(outfmt='6', extra_args=extra_args)

    def format_results(self):
        """Format TCDB DIAMOND results"""
        def extra_processing(df):
            """Additional processing for TCDB results"""
            if 'TCDB ID' in df.columns:
                df['TCDB ID'] = df['TCDB ID'].apply(lambda x: x.split(' ')[0].split('|')[-1] if isinstance(x, str) else x)
            df['Database'] = 'TC'
            
        super().format_results(TCDB_COLUMN_NAMES, extra_processing)