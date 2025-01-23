import pandas as pd
import os

class CGCFinder:
    def __init__(self,config):
        self.config = config
        self.filename = config.get('filename')
        self.num_null_gene = config.get('num_null_gene', 2)
        self.base_pair_distance = config.get('base_pair_distance', 15000)
        self.use_null_genes = config.get('use_null_genes', True)
        self.use_distance = config.get('use_distance', False)
        self.additional_genes = config.get('additional_genes')
        self.output_dir = config.get('output_dir')

    def read_gff(self):
        """Read GFF file using Pandas and extract required information."""
        cols = ['Contig ID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        self.df = pd.read_csv(self.filename, sep='\t', names=cols, comment='#')
        self.df['CGC_annotation'] = self.df['attributes'].apply(lambda x: dict(item.split('=') for item in x.split(';')).get('CGC_annotation', ''))
        self.df['Protein_ID'] = self.df['attributes'].apply(lambda x: dict(item.split('=') for item in x.split(';')).get('protein_id', ''))
        self.df = self.df[['Contig ID', 'start', 'end', 'strand', 'CGC_annotation', 'Protein_ID']]

    def mark_signature_genes(self):
        """Mark signature genes based on their annotations."""
        core_sig_types = ['CAZyme']
        self.df['is_core'] = self.df['CGC_annotation'].str.contains('|'.join(core_sig_types))
        self.df['is_additional'] = self.df['CGC_annotation'].str.contains('|'.join(self.additional_genes))
        self.df['is_signature'] = self.df['is_core'] | self.df['is_additional']

    def find_cgc_clusters(self):
        """Identify CGC clusters using vectorized operations within the same contig."""
        clusters = []
        cgc_id = 1
        
        for contig, contig_df in self.df.groupby('Contig ID'):
            sig_indices = contig_df[contig_df['is_signature']].index
            last_index = None
            start_index = None

            for i in sig_indices:
                if last_index is None:
                    start_index = last_index = i
                    continue

                distance_valid = (contig_df.loc[i, 'start'] - contig_df.loc[last_index, 'end'] <= self.base_pair_distance) if self.use_distance else True
                null_gene_count = (i - last_index - 1)
                null_gene_valid = (null_gene_count <= self.num_null_gene) if self.use_null_genes else True

                if distance_valid and null_gene_valid:
                    last_index = i
                else:
                    cluster_df = contig_df.loc[start_index:last_index]
                    if self.validate_cluster(cluster_df):
                        clusters.append(self.process_cluster(cluster_df, cgc_id))
                        cgc_id += 1
                    start_index = last_index = i

            cluster_df = contig_df.loc[start_index:last_index]
            if self.validate_cluster(cluster_df):
                clusters.append(self.process_cluster(cluster_df, cgc_id))
                cgc_id += 1

        return clusters

    def validate_cluster(self, cluster_df):
        """Validate if the cluster meets the defined CGC criteria."""
        has_core = cluster_df['is_core'].any()
        has_additional = cluster_df['is_additional'].any()
        return (has_core and has_additional) or (has_core and cluster_df['is_core'].sum() > 1)

    def process_cluster(self, cluster_df, cgc_id):
        """Format cluster data for output."""
        return [{
            'CGC#': f'CGC{cgc_id}',
            'Gene Type': gene['CGC_annotation'].split('|')[0],
            'Contig ID': gene['Contig ID'],
            'Protein ID': gene['Protein_ID'],
            'Gene Start': gene['start'],
            'Gene Stop': gene['end'],
            'Gene Strand': gene['strand'],
            'Gene Annotation': gene['CGC_annotation']
        } for _, gene in cluster_df.iterrows()]

    def output_clusters(self, clusters):
        """Output CGC clusters to a TSV file."""
        rows = []
        for cluster in clusters:
            rows.extend(cluster)
        df_output = pd.DataFrame(rows)
        df_output.to_csv(os.path.join(self.output_dir, 'cgc_standard_out.tsv'), sep='\t', index=False)
        print(f"CGC clusters have been written to {self.output_dir}cgc_standard_out.tsv")

# if  __name__ == '__main__':
#     config = {
#         'filename': 'cgc.gff',
#         'num_null_gene': 2,
#         'base_pair_distance': 15000,
#         'use_null_genes': True,
#         'use_distance': False,
#         'additional_genes': ['tcDoms', 'tf-1']  # 自定义附加基因类型
#     }
#     cgc_finder = CGCFinder(**config)
#     cgc_finder.read_gff()
#     cgc_finder.mark_signature_genes()
#     clusters = cgc_finder.find_cgc_clusters()
#     cgc_finder.output_clusters(clusters)
