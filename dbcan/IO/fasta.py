import os
import logging
from multiprocessing.pool import ThreadPool
import pyrodigal
import contextlib
import gzip
from Bio import SeqIO
from dbcan.parameter import GeneralConfig

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class BaseProcessor:
    """Base processor for sequence data with template method pattern"""
    
    def __init__(self, config: GeneralConfig):
        """Initialize processor with configuration"""
        self.config = config
        self._setup_processor()
        
    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        # Setup common attributes
        self.input_raw_data = self._derive_input_path()
        self.output_dir = self._derive_output_dir()
        
        # Ensure output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        # Additional setup from subclasses
        self._additional_setup()
        
        # Validate attributes
        self._validate_attributes()
    
    def _derive_input_path(self):
        """Derive input file path"""
        if not hasattr(self.config, 'input_raw_data') or not self.config.input_raw_data:
            raise ValueError("Input data path not specified in configuration")
        return self.config.input_raw_data
    
    def _derive_output_dir(self):
        """Derive output directory path"""
        if not hasattr(self.config, 'output_dir') or not self.config.output_dir:
            raise ValueError("Output directory not specified in configuration")
        return self.config.output_dir
    
    def _additional_setup(self):
        """Additional setup steps - can be overridden by subclasses"""
        pass
    
    def _validate_attributes(self):
        """Validate processor attributes"""
        # Check if input file exists
        if not os.path.exists(self.input_raw_data):
            raise FileNotFoundError(f"Input file not found: {self.input_raw_data}")
        
        # Check if input file is empty
        if os.path.getsize(self.input_raw_data) == 0:
            raise ValueError(f"Input file is empty: {self.input_raw_data}")

    def parse(self, path):
        """Generic FASTA parser that handles gzipped files"""
        def zopen(f, mode="r"):
            return gzip.open(f, mode) if f.endswith(".gz") else open(f, mode)

        try:
            with contextlib.ExitStack() as ctx:
                file = ctx.enter_context(zopen(path, "rt"))
                id_, desc, seq = None, None, []
                for line in file:
                    if line.startswith(">"):
                        if id_ is not None:
                            yield (id_, "".join(seq), desc)
                        fields = line[1:].strip().split(maxsplit=1)
                        id_ = fields[0] if fields else ""
                        desc = fields[1] if len(fields) > 1 else ""
                        seq = []
                    else:
                        seq.append(line.strip())
                if id_ is not None:
                    yield (id_, "".join(seq), desc)
        except Exception as e:
            logging.error(f"Error parsing file {path}: {str(e)}")
            raise

    def process_input(self):
        """Process input data - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement process_input method")
        
    def _verify_output(self, output_path, expected_content=True):
        """Verify that output file exists and has expected content"""
        if not os.path.exists(output_path):
            logging.error(f"Output file was not created: {output_path}")
            return False
            
        if expected_content and os.path.getsize(output_path) == 0:
            logging.error(f"Output file is empty: {output_path}")
            return False
            
        return True


class ProkProcessor(BaseProcessor):
    """Processor for prokaryotic genomic data"""
    
    def _additional_setup(self):
        """Set up attributes specific to prokaryotic processing"""
        self.output_faa = self._derive_output_faa()
        self.output_gff = self._derive_output_gff()
        self.threads = self._derive_threads()
        
    def _derive_output_faa(self):
        """Derive path for protein output file"""
        return os.path.join(self.output_dir, 'uniInput.faa')
        
    def _derive_output_gff(self):
        """Derive path for GFF output file"""
        return os.path.join(self.output_dir, 'uniInput.gff')
        
    def _derive_threads(self):
        """Derive number of threads to use"""
        return getattr(self.config, 'threads', os.cpu_count())

    def process_input(self):
        """Process prokaryotic genome input"""
        return self.process_fna(False)

    def process_fna(self, is_meta):
        """Process nucleotide sequence data with Pyrodigal"""
        processor_type = "metagenomic" if is_meta else "prokaryotic"
        logging.info(f'Processing {processor_type} genome with Pyrodigal')
        
        try:
            # Initialize gene finder
            gene_finder = pyrodigal.GeneFinder(meta=is_meta)
            
            # Parse input sequences
            logging.info(f"Reading sequences from {self.input_raw_data}")
            sequence_data = []
            for record in self.parse(self.input_raw_data):
                sequence_data.append((record[0], bytes(record[1], 'utf-8')))
            
            if not sequence_data:
                logging.error(f"No valid sequences found in {self.input_raw_data}")
                return None, None
                
            logging.info(f"Found {len(sequence_data)} sequences")
            
            # Train model (except for metagenome mode)
            if not is_meta:
                logging.info("Training Pyrodigal on input sequences")
                gene_finder.train(*(seq[1] for seq in sequence_data))
            
            # Find genes with thread pool
            logging.info(f"Finding genes using {self.threads} threads")
            with ThreadPool(self.threads) as pool:
                results = pool.map(gene_finder.find_genes, [seq[1] for seq in sequence_data])
            
            # Write output files
            logging.info(f"Writing protein translations to {self.output_faa}")
            logging.info(f"Writing gene annotations to {self.output_gff}")
            
            with open(self.output_faa, 'w') as prot_file, open(self.output_gff, 'w') as gff_file:
                genes_found = 0
                for (ori_seq_id, _), genes in zip(sequence_data, results):
                    genes.write_gff(gff_file, sequence_id=ori_seq_id)
                    genes.write_translations(prot_file, sequence_id=ori_seq_id)
                    genes_found += len(genes)
            
            logging.info(f"Processed {len(sequence_data)} sequences, found {genes_found} genes")
            
            # Verify output
            if not self._verify_output(self.output_faa):
                return None, None
                
            if not self._verify_output(self.output_gff):
                return None, None
                
            return self.output_faa, self.output_gff
            
        except Exception as e:
            logging.error(f"Error processing {processor_type} genome: {str(e)}")
            import traceback
            traceback.print_exc()
            return None, None


class MetaProcessor(ProkProcessor):
    """Processor for metagenomic data"""
    
    def process_input(self):
        """Process metagenomic input"""
        return self.process_fna(True)


class ProteinProcessor(BaseProcessor):
    """Processor for protein sequence data"""
    
    def _additional_setup(self):
        """Set up attributes specific to protein processing"""
        self.output_faa = self._derive_output_faa()
        
    def _derive_output_faa(self):
        """Derive path for protein output file"""
        return os.path.join(self.output_dir, 'uniInput.faa')

    def process_input(self):
        """Process protein sequence input"""
        logging.info('Processing protein sequences')
        
        try:
            # Determine if file is compressed
            input_opener = gzip.open if self.input_raw_data.endswith('.gz') else open
            
            record_count = 0
            
            # Process sequences
            logging.info(f"Reading protein sequences from {self.input_raw_data}")
            with input_opener(self.input_raw_data, "rt") as input_handle, open(self.output_faa, "w") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    # Clean up sequence ID (take only first part before any spaces)
                    new_id = record.id.split()[0]
                    record.id = new_id
                    record.description = ''
                    
                    # Write to output file
                    SeqIO.write(record, output_handle, "fasta")
                    record_count += 1
            
            if record_count == 0:
                logging.warning(f"No valid protein sequences found in {self.input_raw_data}")
                return None
                
            logging.info(f"Processed {record_count} protein sequences to {self.output_faa}")
            
            # Verify output
            if not self._verify_output(self.output_faa):
                return None
                
            return self.output_faa
            
        except Exception as e:
            logging.error(f"Error processing protein sequences: {str(e)}")
            import traceback
            traceback.print_exc()
            return None


def get_processor(config: GeneralConfig):
    """Factory function to get the appropriate processor based on mode"""
    mode = config.mode
    
    try:
        if mode == 'prok':
            return ProkProcessor(config)
        elif mode == 'meta':
            return MetaProcessor(config)
        elif mode == 'protein':
            return ProteinProcessor(config)
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    except Exception as e:
        logging.error(f"Error creating processor for mode '{mode}': {str(e)}")
        raise
