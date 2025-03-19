import pandas as pd
import numpy as np
import os
import logging
from dbcan.parameter import CGCFinderConfig
from dbcan.constants import GFF_COLUMNS

class CGCFinder:
    """CGCFinder"""
    
    def __init__(self, config: CGCFinderConfig):
        """Initialize the CGCFinder with configuration."""
        self.config = config
        self._setup_processor()

    def _setup_processor(self):
        """setup the processor with derived attributes"""
        # 基本属性
        self.output_dir = self._derive_output_dir()
        self.filename = self._derive_filename()
        
        # CGC查找参数
        self.num_null_gene = self._derive_num_null_gene()
        self.base_pair_distance = self._derive_base_pair_distance()
        self.use_null_genes = self._derive_use_null_genes()
        self.use_distance = self._derive_use_distance()
        self.additional_genes = self._derive_additional_genes()
        
        # 验证必需的属性
        self._validate_attributes()
    
    def _validate_attributes(self):
        """check if required attributes are set and output directory exists"""
        required_attrs = ['output_dir', 'filename']
        
        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")
        
        # ensure output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)
            logging.info(f"Created output directory: {self.output_dir}")
    
    def _derive_output_dir(self):
        """generate output directory from config"""
        return self.config.output_dir
    
    def _derive_filename(self):
        """generate filename from config"""
        return os.path.join(self.output_dir, 'cgc.gff')
    
    def _derive_num_null_gene(self):
        """generate the maximum number of null genes allowed between signature genes"""
        return getattr(self.config, 'num_null_gene', 2)
    
    def _derive_base_pair_distance(self):
        """generate the maximum base pair distance between signature genes"""
        return getattr(self.config, 'base_pair_distance', 15000)
    
    def _derive_use_null_genes(self):
        """control whether to consider null genes in the distance calculation"""
        return getattr(self.config, 'use_null_genes', True)
    
    def _derive_use_distance(self):
        """consider distance between signature genes when identifying clusters"""
        return getattr(self.config, 'use_distance', False)
    
    def _derive_additional_genes(self):
        """generate additional genes to be considered as signature genes"""
        return getattr(self.config, 'additional_genes', ['TC', 'TF', 'STP'])

    def read_gff(self):
        """read GFF file and extract relevant information"""
        try:
            if not os.path.exists(self.filename):
                logging.error(f"GFF file not found: {self.filename}")
                return False
                
            logging.info(f"Reading GFF file: {self.filename}")
            self.df = pd.read_csv(self.filename, sep='\t', names=GFF_COLUMNS, comment='#')
            
            # extract relevant columns
            self.df['CGC_annotation'] = self.df['attributes'].apply(
                lambda x: dict(item.split('=') for item in x.split(';') if '=' in item).get('CGC_annotation', '')
            )
            self.df['Protein_ID'] = self.df['attributes'].apply(
                lambda x: dict(item.split('=') for item in x.split(';') if '=' in item).get('protein_id', '')
            )
            self.df = self.df[['Contig ID', 'start', 'end', 'strand', 'CGC_annotation', 'Protein_ID']]
            logging.info(f"Loaded {len(self.df)} records from GFF file")
            return True
        except Exception as e:
            logging.error(f"Error reading GFF file: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def mark_signature_genes(self):
        """annotate signature genes in the dataframe"""
        try:
            if not hasattr(self, 'df') or self.df.empty:
                logging.error("No GFF data loaded. Run read_gff() first.")
                return False
                
            core_sig_types = ['CAZyme']
            self.df['is_core'] = self.df['CGC_annotation'].str.contains('|'.join(core_sig_types), na=False)
            self.df['is_additional'] = self.df['CGC_annotation'].str.contains('|'.join(self.additional_genes), na=False)
            self.df['is_signature'] = self.df['is_core'] | self.df['is_additional']
            
            sig_gene_count = self.df['is_signature'].sum()
            logging.info(f"Marked {sig_gene_count} signature genes ({self.df['is_core'].sum()} core, {self.df['is_additional'].sum()} additional)")
            return True
        except Exception as e:
            logging.error(f"Error marking signature genes: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def find_cgc_clusters(self):
        """identify CGC clusters based on the defined criteria"""
        try:
            if not hasattr(self, 'df') or self.df.empty:
                logging.error("No GFF data loaded or no signature genes marked.")
                return []
                
            clusters = []
            cgc_id = 1
            
            logging.info(f"Finding CGC clusters using {self.num_null_gene} max null genes, {self.base_pair_distance} bp distance")
            logging.info(f"Use null genes: {self.use_null_genes}, Use distance: {self.use_distance}")

            for contig, contig_df in self.df.groupby('Contig ID'):
                sig_indices = contig_df[contig_df['is_signature']].index.to_numpy()
                
                if len(sig_indices) < 2:
                    continue  # need at least 2 signature genes to form a cluster
                    
                starts = contig_df.loc[sig_indices, 'start'].to_numpy()
                ends = contig_df.loc[sig_indices, 'end'].to_numpy()

                last_index = None
                start_index = None

                for i, sig_index in enumerate(sig_indices):
                    if last_index is None:
                        start_index = last_index = sig_index
                        continue

                    distance_valid = (starts[i] - ends[i - 1] <= self.base_pair_distance) if self.use_distance else True
                    null_gene_count = (sig_index - last_index - 1)
                    null_gene_valid = (null_gene_count <= self.num_null_gene) if self.use_null_genes else True

                    if distance_valid and null_gene_valid:
                        last_index = sig_index
                    else:
                        cluster_df = contig_df.loc[start_index:last_index]
                        if self.validate_cluster(cluster_df):
                            clusters.append(self.process_cluster(cluster_df, cgc_id))
                            cgc_id += 1
                        start_index = last_index = sig_index

                # process the last cluster if it exists
                if last_index is not None and start_index is not None:
                    cluster_df = contig_df.loc[start_index:last_index]
                    if self.validate_cluster(cluster_df):
                        clusters.append(self.process_cluster(cluster_df, cgc_id))
                        cgc_id += 1

            logging.info(f"Found {len(clusters)} CGC clusters")
            return clusters
        except Exception as e:
            logging.error(f"Error finding CGC clusters: {str(e)}")
            import traceback
            traceback.print_exc()
            return []

    def validate_cluster(self, cluster_df):
        """validate if a cluster meets the criteria for being a CGC"""
        has_core = cluster_df['is_core'].any()
        has_additional = cluster_df['is_additional'].any()
        return (has_core and has_additional) or (has_core and cluster_df['is_core'].sum() > 1)

    def process_cluster(self, cluster_df, cgc_id):
        """format a cluster for output"""
        return [{
            'CGC#': f'CGC{cgc_id}',
            'Gene Type': gene['CGC_annotation'].split('|')[0] if '|' in gene['CGC_annotation'] else 'null',
            'Contig ID': gene['Contig ID'],
            'Protein ID': gene['Protein_ID'],
            'Gene Start': gene['start'],
            'Gene Stop': gene['end'],
            'Gene Strand': gene['strand'],
            'Gene Annotation': gene['CGC_annotation']
        } for _, gene in cluster_df.iterrows()]

    def output_clusters(self, clusters):
        """export identified CGC clusters to a TSV file"""
        try:
            if not clusters:
                logging.warning("No CGC clusters found to output")
                # generate empty file
                empty_df = pd.DataFrame(columns=['CGC#', 'Gene Type', 'Contig ID', 'Protein ID', 
                                               'Gene Start', 'Gene Stop', 'Gene Strand', 'Gene Annotation'])
                output_path = os.path.join(self.output_dir, 'cgc_standard_out.tsv')
                empty_df.to_csv(output_path, sep='\t', index=False)
                logging.info(f"Empty CGC output file created at {output_path}")
                return
                
            rows = []
            for cluster in clusters:
                rows.extend(cluster)
                
            df_output = pd.DataFrame(rows)
            output_path = os.path.join(self.output_dir, 'cgc_standard_out.tsv')
            df_output.to_csv(output_path, sep='\t', index=False)
            logging.info(f"CGC clusters have been written to {output_path}")
        except Exception as e:
            logging.error(f"Error outputting CGC clusters: {str(e)}")
            import traceback
            traceback.print_exc()

    def run(self):
        """run the CGCFinder"""
        if not self.read_gff():
            return False
        if not self.mark_signature_genes():
            return False
        clusters = self.find_cgc_clusters()
        self.output_clusters(clusters)
        logging.info("CGCFinder run completed")
        return True
