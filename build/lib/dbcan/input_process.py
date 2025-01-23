import os
import logging
from multiprocessing.pool import ThreadPool
import pyrodigal
import contextlib
import gzip
from Bio import SeqIO




logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class InputProcessor:

    #use config to initialize the processor
    def __init__(self, config):
        self.config = config
        self.mode = config.get('mode') 
        self.input_format = config.get('input_format')
        self.input_raw_data = config.get('input_raw_data')
        self.output_dir = config.get('output_dir')
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    #parse the input fasta file
    def parse(self, path):
        def zopen(f, mode="r"):
            return gzip.open(f, mode) if f.endswith(".gz") else open(f, mode)
        
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

    #process the input file based on the mode
    def process_input(self):
        if self.mode == 'prok':
            faa, gff = self.process_fna(False)
        elif self.mode == 'meta':
            faa, gff = self.process_fna(True)
        elif self.mode == 'protein':
            faa = self.process_protein()

    def process_fna(self, is_meta):
        logging.info(f'Processing {"metagenomic" if is_meta else "prokaryotic"} genome with Pyrodigal')
        gene_finder = pyrodigal.GeneFinder(meta=is_meta)
        faa_path = os.path.join(self.output_dir, 'uniInput.faa')
        gff_path = os.path.join(self.output_dir, 'uniInput.gff')
        sequence_data = [(record[0], bytes(record[1], 'utf-8')) for record in self.parse(self.input_raw_data)]
        #single mode need to be trained before finding genes
        if not is_meta:
            gene_finder.train(*(seq[1] for seq in sequence_data))

        #use threads to speed up the process
        with ThreadPool(os.cpu_count()) as pool:
            results = pool.map(gene_finder.find_genes, [seq[1] for seq in sequence_data])

        with open(faa_path, 'w') as prot_file, open(gff_path, 'w') as out_file:
            for (ori_seq_id, _), genes in zip(sequence_data, results):
                genes.write_gff(out_file, sequence_id=ori_seq_id)
                genes.write_translations(prot_file, sequence_id=ori_seq_id)

        return faa_path, gff_path
        #directly copy input files file to output directory
    def process_protein(self):
        logging.info('Processing protein sequences based on input format: {}'.format(self.input_format))
        faa_path = os.path.join(self.output_dir, 'uniInput.faa')

        with open(self.input_raw_data, "r") as input_handle, open(faa_path, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                if self.input_format == "NCBI":
                    new_id = record.id.split()[0]  
                elif self.input_format == "JGI":
                    new_id = record.id.split('|')[2]  
                else:
                    logging.warning(f"Unrecognized input format: {self.input_format}")
                    new_id = record.id

                record.id = new_id
                record.description = '' 
                SeqIO.write(record, output_handle, "fasta")

        return faa_path

# if __name__ == '__main__':
#     config = {
#         'mode': 'prok',  # or 'meta' or 'protein'
#         'input_file': 'input/EscheriaColiK12MG1655.fna',
#         'output_dir': 'output'
#     }
#     processor = InputProcessor(config)
#     processor.process_input()