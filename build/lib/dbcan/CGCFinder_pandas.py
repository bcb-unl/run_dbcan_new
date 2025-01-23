import pandas as pd

def read_gff(filename):
    """Read GFF file using Pandas and extract required information."""
    cols = ['Contig ID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = pd.read_csv(filename, sep='\t', names=cols, comment='#')
    df['CGC_annotation'] = df['attributes'].apply(lambda x: dict(item.split('=') for item in x.split(';')).get('CGC_annotation', ''))
    df['Protein_ID'] = df['attributes'].apply(lambda x: dict(item.split('=') for item in x.split(';')).get('protein_id', ''))
    return df[['Contig ID', 'start', 'end', 'strand', 'CGC_annotation', 'Protein_ID']]

def mark_signature_genes(df):
    """Mark signature genes based on their annotations."""
    core_sig_types = ['CAZyme']
    additional_sig_types = ['TC', 'TF', 'STP']
    df['is_core'] = df['CGC_annotation'].str.contains('|'.join(core_sig_types))
    df['is_additional'] = df['CGC_annotation'].str.contains('|'.join(additional_sig_types))
    df['is_signature'] = df['is_core'] | df['is_additional']
    return df

def find_cgc_clusters(df, num_null_gene=2, base_pair_distance=15000, use_null_genes=True, use_distance=False):
    """Identify CGC clusters using vectorized operations within the same contig."""
    clusters = []
    cgc_id = 1
    
    for contig, contig_df in df.groupby('Contig ID'):
        sig_indices = contig_df[contig_df['is_signature']].index
        last_index = None
        start_index = None

        for i in sig_indices:
            if last_index is None:
                start_index = last_index = i
                continue

            # Apply distance and null gene conditions
            distance_valid = (contig_df.loc[i, 'start'] - contig_df.loc[last_index, 'end'] <= base_pair_distance) if use_distance else True
            null_gene_count = (i - last_index - 1)
            null_gene_valid = (null_gene_count <= num_null_gene) if use_null_genes else True

            if distance_valid and null_gene_valid:
                last_index = i
            else:
                # Check if the current cluster meets the criteria
                cluster_df = contig_df.loc[start_index:last_index]
                if validate_cluster(cluster_df):
                    clusters.append(process_cluster(cluster_df, cgc_id))
                    cgc_id += 1
                start_index = last_index = i

        # Check the last cluster
        cluster_df = contig_df.loc[start_index:last_index]
        if validate_cluster(cluster_df):
            clusters.append(process_cluster(cluster_df, cgc_id))
            cgc_id += 1

    return clusters

def validate_cluster(cluster_df):
    """Validate if the cluster meets the defined CGC criteria."""
    has_core = cluster_df['is_core'].any()
    has_additional = cluster_df['is_additional'].any()
    return (has_core and has_additional) or (has_core and cluster_df['is_core'].sum() > 1)

def process_cluster(cluster_df, cgc_id):
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

def output_clusters(clusters):
    """Output CGC clusters to a TSV file."""
    rows = []
    for cluster in clusters:
        rows.extend(cluster)
    df_output = pd.DataFrame(rows)
    df_output.to_csv('output.tsv', sep='\t', index=False)

# Usage example
filename = 'cgc.gff'
df = read_gff(filename)
df = mark_signature_genes(df)
clusters = find_cgc_clusters(df)
output_clusters(clusters)
