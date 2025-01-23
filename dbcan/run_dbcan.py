import argparse
import os
import psutil
import pandas as pd
from dbcan.dbCAN_database import DBDownloader
from dbcan.input_process import InputProcessor
from dbcan.diamond import DiamondProcessor
from dbcan.pyhmmer_search import PyHMMERProcessor
from dbcan.process_dbcan_sub import DBCANProcessor
from dbcan.OverviewGenerator import OverviewGenerator
from dbcan.generate_cgc_gff import GFFProcessor
from dbcan.CGCFinder_pandas_class import CGCFinder
from dbcan.cgc_substrate_prediction import cgc_substrate_prediction
from dbcan.syntenic_plot import  syntenic_plot_allpairs


def run_dbCAN_database(args):

    dbCAN_database_config = {'db_dir': args.db_dir}
    downloader = DBDownloader(dbCAN_database_config)
    downloader.download_file()
    downloader.extract_tar_file()

def run_dbCAN_input_process(args):
    dbCAN_input_process_config = {
                'input_raw_data': args.input_raw_data, 
                'mode': args.mode, 
                'output_dir': args.output_dir,
                'input_format': args.input_format}
    input_processor = InputProcessor(dbCAN_input_process_config)
    input_processor.process_input()

def run_dbCAN_CAZyme_annotation(args):
    method  = args.methods
    diamond_config = {
            'diamond_db': os.path.join(args.db_dir, 'CAZy.dmnd'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_file': os.path.join(args.output_dir, 'diamond_results.tsv'),
            'e_value_threshold': args.diamond_evalue,
            'threads': args.diamond_threads,
            'verbose_option': args.diamond_verbose_option
        }
    hmm_config = {
            'hmm_file': os.path.join(args.db_dir, 'dbCAN.hmm'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_file': os.path.join(args.output_dir, 'dbCAN_hmm_results.tsv'),
            'e_value_threshold': args.dbcan_hmm_evalue,
            'coverage_threshold': args.dbcan_hmm_coverage,
            'hmmer_cpu': args.hmmer_cpu
        }
    dbcan_sub_config = {
            'hmm_file': os.path.join(args.db_dir, 'dbCAN_sub.hmm'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'e_value_threshold': args.dbcansub_hmm_evalue,
            'coverage_threshold': args.dbcansub_hmm_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_file': os.path.join(args.output_dir, 'dbCAN-sub.substrate.tsv'),
            'mapping_file': os.path.join(args.db_dir, 'fam-substrate-mapping.tsv')
        }
    generate_overview_config = {'output_dir': args.output_dir}

    if 'diamond' in method:
        diamond_processor = DiamondProcessor(diamond_config)
        diamond_processor.run_diamond()
        diamond_processor.format_results()
    if 'hmm' in method:
        hmm_processor = PyHMMERProcessor(hmm_config)
        hmm_processor.hmmsearch()
    if 'dbCANsub' in method:
        dbcan_sub_processor = PyHMMERProcessor(dbcan_sub_config)
        dbcan_sub_processor.hmmsearch()
        parser_dbcan_sub = DBCANProcessor(dbcan_sub_config)
        parser_dbcan_sub.process_dbcan_sub()
    overview_generator = OverviewGenerator(generate_overview_config)
    overview_generator.run()



def run_CGC_annotation_preprocess(args):
    CGC_info_config = { 
            'input_total_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_dir': args.output_dir,
            'cazyme_overview': os.path.join(args.output_dir, 'overview.tsv'),
            'cgc_sig_file': os.path.join(args.output_dir, 'total_cgc_info.tsv'),
            'input_gff': args.input_gff,
            'output_gff': os.path.join(args.output_dir, 'cgc.gff'),
            'gff_type': args.input_gff_format
        }

    tc_config = {
            'diamond_db': os.path.join(args.db_dir, 'tcdb.dmnd'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'TC_results.tsv'),
            'e_value_threshold': args.tc_evalue,
            'coverage_threshold_tc': args.tc_coverage,
            'threads': args.diamond_threads,
            'verbose_option': args.diamond_verbose_option
        }

    tf_config = {
            'hmm_file': os.path.join(args.db_dir, 'TF.hmm'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'TF_results.tsv'),
            'e_value_threshold': args.tf_evalue,
            'coverage_threshold': args.tf_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_dir': args.output_dir
        }

    stp_config = {
            'hmm_file': os.path.join(args.db_dir, 'STP.hmm'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'STP_results.tsv'),
            'e_value_threshold': args.stp_evalue,
            'coverage_threshold': args.stp_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_dir': args.output_dir
        }

    cgc_infoProcessor = GFFProcessor(CGC_info_config)
    cgc_infoProcessor.generate_non_cazyme_faa()
#    tc_hmm_processor   = PyHMMERProcessor(tc_config)
#    tc_hmm_processor.hmmsearch()   #used to process tcDoms but now use diamond

    tc_diamond_processor = DiamondProcessor(tc_config)
    tf_hmm_processor   = PyHMMERProcessor(tf_config)
    stp_hmm_processor  = PyHMMERProcessor(stp_config)
    tc_diamond_processor.run_tcdb_diamond()
    tc_diamond_processor.format_results_tcdb()
    tf_hmm_processor.hmmsearch()
    stp_hmm_processor.hmmsearch()

    columns = ['Annotate Name', 'Annotate Length', 'Target Name', 'Target Length', 'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To', 'Coverage', 'Annotate File Name']
    tc_df = pd.read_csv(tc_config['output_file'], names=columns, header=0, sep='\t')
    tf_df = pd.read_csv(tf_config['output_file'], names=columns, header=0, sep='\t')
    stp_df = pd.read_csv(stp_config['output_file'], names=columns, header=0, sep='\t')
    total_function_annotation_df = pd.concat([tc_df, tf_df, stp_df], ignore_index=True)
    tf_hmm_processor.filter_overlaps(total_function_annotation_df).to_csv(os.path.join(tf_config["output_dir"], 'total_cgc_info.tsv'), index=False, sep='\t')
    cgc_infoProcessor.process_gff()

def run_CGC_annotation(args):

    cgc_config = {
            'output_dir': args.output_dir,
            'filename': os.path.join(args.output_dir, 'cgc.gff'),
            'num_null_gene': args.num_null_gene,
            'base_pair_distance': args.base_pair_distance,
            'use_null_genes': args.use_null_genes,
            'use_distance': args.use_distance,
            'additional_genes': args.additional_genes
        }

    cgc_finder = CGCFinder(cgc_config)
    cgc_finder.read_gff()
    cgc_finder.mark_signature_genes()
    clusters = cgc_finder.find_cgc_clusters()   
    cgc_finder.output_clusters(clusters)


def main():
    parser = argparse.ArgumentParser(description='CAZyme analysis workflow management.')

    subparsers = parser.add_subparsers(dest='command', help='The function to run: \n'
'database - Downloads and prepares the dbCAN database files.\n'
'input_process - Processes the input fasta files for annotation.\n'
'CAZyme_annotation - Runs the CAZyme annotation using specified methods.\n'
'CGC_info - Prepares the input files for CGC annotation.\n'
'CGC_annotation - Runs the CGC annotation using specified methods.\n'
'CGC_substrate_prediction - Predicts the substrates of CGCs.\n'
'CGC_substrate_plot - Generates syntenic plots for CGC substrate prediction. \n'
'easy_CAZyme  - Runs the easy CAZyme annotation pipeline.\n'
'easy_CGC  - Runs the easy CGC annotation pipeline.\n'
'easy_substrate - Runs the easy CGC substrate prediction pipeline.',  required=True)
    
    #subparser for database
    parser_database = subparsers.add_parser('database')
    parser_database.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.')

    #subparser for input_process
    parser_input_process = subparsers.add_parser('input_process')
    parser_input_process.add_argument('--input_raw_data', help='Specify the input fasta file for data preprocessing.', required=True)
    parser_input_process.add_argument('--mode', help='Check the mode of the input file: prok, meta, or protein',required=True)
    parser_input_process.add_argument('--input_format', default='NCBI', choices=['NCBI', 'JGI'], help='Specify the input format for protein sequences, only needed when mode is protein')
    parser_input_process.add_argument('--output_dir', default='./output', help='Output folder.',required=True)
    parser_input_process.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.',required=True)

    #subparser for CAZyme_annotation
    parser_caz_annotation = subparsers.add_parser('CAZyme_annotation')
    parser_caz_annotation.add_argument('--methods', nargs='+', choices=['diamond', 'hmm', 'dbCANsub'],
        default=['diamond', 'hmm', 'dbCANsub'],help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.')
    parser_caz_annotation.add_argument('--diamond_evalue', type=float, default=1e-102, help='E-value threshold for diamond annotation.')
    parser_caz_annotation.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser_caz_annotation.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')
    parser_caz_annotation.add_argument('--dbcan_hmm_evalue', type=float, default=1e-15, help='E-value threshold for HMM annotation.')
    parser_caz_annotation.add_argument('--dbcan_hmm_coverage', type=float, default=0.35, help='Coverage threshold for HMM annotation.')
    parser_caz_annotation.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')
    parser_caz_annotation.add_argument('--dbcansub_hmm_evalue', type=float, default=1e-15, help='E-value threshold for dbCANsub annotation.')
    parser_caz_annotation.add_argument('--dbcansub_hmm_coverage', type=float, default=0.35, help='Coverage threshold for dbCANsub annotation.')
    parser_caz_annotation.add_argument('--output_dir', default='./output', help='Output folder.',required=True)
    parser_caz_annotation.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.',required=True)

    #subparser for CGC_info
    parser_cgc_info = subparsers.add_parser('CGC_info')
    parser_cgc_info.add_argument('--input_gff', help='Specify the input gff file for CGC annotation.', required=True)
    parser_cgc_info.add_argument('--input_gff_format', default='NCBI_prok', choices=['NCBI_euk', 'JGI', 'NCBI_prok', 'prodigal'], help='Specify the input format for protein sequences.',required=True)
    parser_cgc_info.add_argument('--tc_evalue', type=float, default=1e-15, help='E-value threshold for TC annotation.')
    parser_cgc_info.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser_cgc_info.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')
    parser_cgc_info.add_argument('--tc_coverage', type=float, default=0.35, help='Coverage threshold for TC annotation.')
    parser_cgc_info.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')
    parser_cgc_info.add_argument('--tf_evalue', type=float, default=1e-15, help='E-value threshold for TF annotation.')
    parser_cgc_info.add_argument('--tf_coverage', type=float, default=0.35, help='Coverage threshold for TF annotation.')
    parser_cgc_info.add_argument('--stp_evalue', type=float, default=1e-15, help='E-value threshold for STP annotation.')
    parser_cgc_info.add_argument('--stp_coverage', type=float, default=0.35, help='Coverage threshold for STP annotation.')
    parser_cgc_info.add_argument('--output_dir', default='./output', help='Output folder.',required=True)
    parser_cgc_info.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.',required=True)

    #subparser for CGC_annotation
    parser_cgc_annotation = subparsers.add_parser('CGC_annotation')
    parser_cgc_annotation.add_argument('--additional_genes', nargs='+', default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')
    parser_cgc_annotation.add_argument('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')
    parser_cgc_annotation.add_argument('--base_pair_distance', type=int, default=15000, help='Base pair distance of sig genes for CGC annotation.')
    parser_cgc_annotation.add_argument('--use_null_genes', action='store_true', default=True, help='Use null genes in CGC annotation.')
    parser_cgc_annotation.add_argument('--use_distance', action='store_true', default=False, help='Use base pair distance in CGC annotation.')
    parser_cgc_annotation.add_argument('--output_dir', default='./output', help='Output folder.',required=True)


    #subparser for CGC_substrate_prediction
    parser_cgc_substrate = subparsers.add_parser('CGC_substrate_prediction')
    group = parser_cgc_substrate.add_argument_group('general optional arguments')
    group.add_argument('-i','--input',help="input file: dbCAN3 output folder", required=True)
    group.add_argument('--pul',help="dbCAN-PUL PUL.faa")
    group.add_argument('-o','--out',default="substrate.out")
    group.add_argument('-w','--workdir',type=str,default=".")
    group.add_argument('-rerun','--rerun',type=bool,default=False)
    group.add_argument('-env','--env',type=str,default="local")
    group.add_argument('-odbcan_sub','--odbcan_sub', help="output dbcan_sub prediction intermediate result?")
    group.add_argument('-odbcanpul','--odbcanpul',type=bool,default=True,help="output dbCAN-PUL prediction intermediate result?")
    group.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.',required=True)
    group.add_argument('--output_dir', default='./output', help='Output folder.',required=True)

    ### paramters to identify a homologous PUL
    ### including blastp evalue,number of CAZyme pair, number of pairs, extra pair, bitscore_cutoff, uniq query cgc gene hits.
    ### uniq PUL gene hits. identity cutoff. query coverage cutoff.
    group1 = parser_cgc_substrate.add_argument_group('dbCAN-PUL homologous conditons', 'how to define homologous gene hits and PUL hits')
    group1.add_argument('-upghn','--uniq_pul_gene_hit_num',default = 2,type=int)
    group1.add_argument('-uqcgn','--uniq_query_cgc_gene_num',default = 2,type=int)
    group1.add_argument('-cpn','--CAZyme_pair_num',default = 1,type=int)
    group1.add_argument('-tpn','--total_pair_num',default = 2,type=int)
    group1.add_argument('-ept','--extra_pair_type',default = None,type=str,help="None[TC-TC,STP-STP]. Some like sigunature hits")
    group1.add_argument('-eptn','--extra_pair_type_num',default ="0",type=str,help="specify signature pair cutoff.1,2")
    group1.add_argument('-iden','--identity_cutoff',default = 0.,type=float,help="identity to identify a homologous hit")
    group1.add_argument('-cov','--coverage_cutoff',default = 0.,type=float,help="query coverage cutoff to identify a homologous hit")
    group1.add_argument('-bsc','--bitscore_cutoff',default = 50,type=float,help="bitscore cutoff to identify a homologous hit")
    group1.add_argument('-evalue','--evalue_cutoff',default = 0.01,type=float,help="evalue cutoff to identify a homologous hit")

    group2 = parser_cgc_substrate.add_argument_group('dbCAN-sub conditons', 'how to define dbsub hits and dbCAN-sub subfamily substrate')
    group2.add_argument('-hmmcov','--hmmcov',default = 0.,type=float)
    group2.add_argument('-hmmevalue','--hmmevalue',default = 0.01,type=float)
    group2.add_argument('-ndsc','--num_of_domains_substrate_cutoff',default = 2,type=int,help="define how many domains share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-npsc','--num_of_protein_substrate_cutoff',default = 2,type=int,help="define how many sequences share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-subs','--substrate_scors',default = 2,type=int,help="each cgc contains with substrate must more than this value")

    #subparser for syntenic plot
    parser_syntenic = subparsers.add_parser('CGC_substrate_plot')
    parser_syntenic.add_argument('-input_sub_out',help='substrate out ')
    parser_syntenic.add_argument('-b','--blastp',help='blastp result for cgc')
    parser_syntenic.add_argument('--cgc', help='cgc_finder output')
    parser_syntenic.add_argument('--db_dir', default="db", help='Database directory',required=True)
    parser_syntenic.add_argument('--output_dir', default='./output', help='Output folder.',required=True)

    #subparser for easy_CAZyme
    parser_easy_caz = subparsers.add_parser('easy_CAZyme')
    parser_easy_caz.add_argument('--input_raw_data', help='Specify the input fasta file for data preprocessing.', required=True)
    parser_easy_caz.add_argument('--mode', help='Check the mode of the input file: prok, meta, or protein', required=True)
    parser_easy_caz.add_argument('--input_format', default='NCBI', choices=['NCBI', 'JGI'], help='Specify the input format for protein sequences.')
    parser_easy_caz.add_argument('--output_dir', default='./output', help='Output folder.',required=True)
    parser_easy_caz.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.',required=True)
    parser_easy_caz.add_argument('--methods', nargs='+', choices=['diamond', 'hmm', 'dbCANsub'],
        default=['diamond', 'hmm', 'dbCANsub'],help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.')
    parser_easy_caz.add_argument('--diamond_evalue', type=float, default=1e-102, help='E-value threshold for diamond annotation.')
    parser_easy_caz.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser_easy_caz.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')
    parser_easy_caz.add_argument('--dbcan_hmm_evalue', type=float, default=1e-15, help='E-value threshold for HMM annotation.')
    parser_easy_caz.add_argument('--dbcan_hmm_coverage', type=float, default=0.35, help='Coverage threshold for HMM annotation.')
    parser_easy_caz.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')
    parser_easy_caz.add_argument('--dbcansub_hmm_evalue', type=float, default=1e-15, help='E-value threshold for dbCANsub annotation.')
    parser_easy_caz.add_argument('--dbcansub_hmm_coverage', type=float, default=0.35, help='Coverage threshold for dbCANsub annotation.')

    #subparser for easy_CGC
    parser_easy_cgc = subparsers.add_parser('easy_CGC')
    parser_easy_cgc.add_argument('--input_raw_data', help='Specify the input fasta file for data preprocessing.', required=True)
    parser_easy_cgc.add_argument('--mode', help='Check the mode of the input file: prok, meta, or protein', required=True)
    parser_easy_cgc.add_argument('--input_format', default='NCBI', choices=['NCBI', 'JGI'], help='Specify the input format for protein sequences.')
    parser_easy_cgc.add_argument('--output_dir', default='./output', help='Output folder.', required=True)
    parser_easy_cgc.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.', required=True)
    parser_easy_cgc.add_argument('--methods', nargs='+', choices=['diamond', 'hmm', 'dbCANsub'],
        default=['diamond', 'hmm', 'dbCANsub'],help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.')
    parser_easy_cgc.add_argument('--diamond_evalue', type=float, default=1e-102, help='E-value threshold for diamond annotation.')
    parser_easy_cgc.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser_easy_cgc.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')
    parser_easy_cgc.add_argument('--dbcan_hmm_evalue', type=float, default=1e-15, help='E-value threshold for HMM annotation.')
    parser_easy_cgc.add_argument('--dbcan_hmm_coverage', type=float, default=0.35, help='Coverage threshold for HMM annotation.')
    parser_easy_cgc.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')
    parser_easy_cgc.add_argument('--dbcansub_hmm_evalue', type=float, default=1e-15, help='E-value threshold for dbCANsub annotation.')
    parser_easy_cgc.add_argument('--dbcansub_hmm_coverage', type=float, default=0.35, help='Coverage threshold for dbCANsub annotation.')
    parser_easy_cgc.add_argument('--input_gff', help='Specify the input gff file for CGC annotation.', required=True)
    parser_easy_cgc.add_argument('--input_gff_format', default='NCBI_prok', choices=['NCBI_euk', 'JGI', 'NCBI_prok', 'prodigal'], help='Specify the input format for protein sequences.', required=True)
    parser_easy_cgc.add_argument('--tc_evalue', type=float, default=1e-15, help='E-value threshold for TC annotation.')
    parser_easy_cgc.add_argument('--tc_coverage', type=float, default=0.35, help='Coverage threshold for TC annotation.')
    parser_easy_cgc.add_argument('--tf_evalue', type=float, default=1e-15, help='E-value threshold for TF annotation.')
    parser_easy_cgc.add_argument('--tf_coverage', type=float, default=0.35, help='Coverage threshold for TF annotation.')
    parser_easy_cgc.add_argument('--stp_evalue', type=float, default=1e-15, help='E-value threshold for STP annotation.')
    parser_easy_cgc.add_argument('--stp_coverage', type=float, default=0.35, help='Coverage threshold for STP annotation.')
    parser_easy_cgc.add_argument('--additional_genes', nargs='+', default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')
    parser_easy_cgc.add_argument('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')
    parser_easy_cgc.add_argument('--base_pair_distance', type=int, default=15000, help='Base pair distance of sig genes for CGC annotation.')
    parser_easy_cgc.add_argument('--use_null_genes', action='store_true', default=True, help='Use null genes in CGC annotation.')
    parser_easy_cgc.add_argument('--use_distance', action='store_true', default=False, help='Use base pair distance in CGC annotation.')

    #subparser for easy_substrate
    parser_easy_sub = subparsers.add_parser('easy_substrate')
    parser_easy_sub.add_argument('--input_raw_data', help='Specify the input fasta file for data preprocessing.', required=True)
    parser_easy_sub.add_argument('--mode', help='Check the mode of the input file: prok, meta, or protein', required=True)
    parser_easy_sub.add_argument('--input_format', default='NCBI', choices=['NCBI', 'JGI'], help='Specify the input format for protein sequences.')
    parser_easy_sub.add_argument('--output_dir', default='./output', help='Output folder.', required=True)
    parser_easy_sub.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.', required=True)
    parser_easy_sub.add_argument('--methods', nargs='+', choices=['diamond', 'hmm', 'dbCANsub'],
        default=['diamond', 'hmm', 'dbCANsub'],help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.')
    parser_easy_sub.add_argument('--diamond_evalue', type=float, default=1e-102, help='E-value threshold for diamond annotation.')
    parser_easy_sub.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser_easy_sub.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')
    parser_easy_sub.add_argument('--dbcan_hmm_evalue', type=float, default=1e-15, help='E-value threshold for HMM annotation.')
    parser_easy_sub.add_argument('--dbcan_hmm_coverage', type=float, default=0.35, help='Coverage threshold for HMM annotation.')
    parser_easy_sub.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')
    parser_easy_sub.add_argument('--dbcansub_hmm_evalue', type=float, default=1e-15, help='E-value threshold for dbCANsub annotation.')
    parser_easy_sub.add_argument('--dbcansub_hmm_coverage', type=float, default=0.35, help='Coverage threshold for dbCANsub annotation.')
    parser_easy_sub.add_argument('--input_gff', help='Specify the input gff file for CGC annotation.', required=True)
    parser_easy_sub.add_argument('--input_gff_format', default='NCBI_prok', choices=['NCBI_euk', 'JGI', 'NCBI_prok', 'prodigal'], help='Specify the input format for protein sequences.', required=True)
    parser_easy_sub.add_argument('--tc_evalue', type=float, default=1e-15, help='E-value threshold for TC annotation.')
    parser_easy_sub.add_argument('--tc_coverage', type=float, default=0.35, help='Coverage threshold for TC annotation.')
    parser_easy_sub.add_argument('--tf_evalue', type=float, default=1e-15, help='E-value threshold for TF annotation.')
    parser_easy_sub.add_argument('--tf_coverage', type=float, default=0.35, help='Coverage threshold for TF annotation.')
    parser_easy_sub.add_argument('--stp_evalue', type=float, default=1e-15, help='E-value threshold for STP annotation.')
    parser_easy_sub.add_argument('--stp_coverage', type=float, default=0.35, help='Coverage threshold for STP annotation.')
    parser_easy_sub.add_argument('--additional_genes', nargs='+', default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')
    parser_easy_sub.add_argument('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')
    parser_easy_sub.add_argument('--base_pair_distance', type=int, default=15000, help='Base pair distance of sig genes for CGC annotation.')
    parser_easy_sub.add_argument('--use_null_genes', action='store_true', default=True, help='Use null genes in CGC annotation.')
    parser_easy_sub.add_argument('--use_distance', action='store_true', default=False, help='Use base pair distance in CGC annotation.')
    group = parser_easy_sub.add_argument_group('general optional arguments')
    group.add_argument('-i','--input',help="input file: dbCAN3 output folder")
    group.add_argument('--pul',help="dbCAN-PUL PUL.faa")
    group.add_argument('-o','--out',default="substrate.out")
    group.add_argument('-w','--workdir',type=str,default=".")
    group.add_argument('-rerun','--rerun',type=bool,default=False)
    group.add_argument('-env','--env',type=str,default="local")
    group.add_argument('-odbcan_sub','--odbcan_sub', help="output dbcan_sub prediction intermediate result?")
    group.add_argument('-odbcanpul','--odbcanpul',type=bool,default=True,help="output dbCAN-PUL prediction intermediate result?")

    ### paramters to identify a homologous PUL
    ### including blastp evalue,number of CAZyme pair, number of pairs, extra pair, bitscore_cutoff, uniq query cgc gene hits.
    ### uniq PUL gene hits. identity cutoff. query coverage cutoff.
    group1 = parser_easy_sub.add_argument_group('dbCAN-PUL homologous conditons', 'how to define homologous gene hits and PUL hits')
    group1.add_argument('-upghn','--uniq_pul_gene_hit_num',default = 2,type=int)
    group1.add_argument('-uqcgn','--uniq_query_cgc_gene_num',default = 2,type=int)
    group1.add_argument('-cpn','--CAZyme_pair_num',default = 1,type=int)
    group1.add_argument('-tpn','--total_pair_num',default = 2,type=int)
    group1.add_argument('-ept','--extra_pair_type',default = None,type=str,help="None[TC-TC,STP-STP]. Some like sigunature hits")
    group1.add_argument('-eptn','--extra_pair_type_num',default ="0",type=str,help="specify signature pair cutoff.1,2")
    group1.add_argument('-iden','--identity_cutoff',default = 0.,type=float,help="identity to identify a homologous hit")
    group1.add_argument('-cov','--coverage_cutoff',default = 0.,type=float,help="query coverage cutoff to identify a homologous hit")
    group1.add_argument('-bsc','--bitscore_cutoff',default = 50,type=float,help="bitscore cutoff to identify a homologous hit")
    group1.add_argument('-evalue','--evalue_cutoff',default = 0.01,type=float,help="evalue cutoff to identify a homologous hit")

    group2 = parser_easy_sub.add_argument_group('dbCAN-sub conditons', 'how to define dbsub hits and dbCAN-sub subfamily substrate')
    group2.add_argument('-hmmcov','--hmmcov',default = 0.,type=float)
    group2.add_argument('-hmmevalue','--hmmevalue',default = 0.01,type=float)
    group2.add_argument('-ndsc','--num_of_domains_substrate_cutoff',default = 2,type=int,help="define how many domains share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-npsc','--num_of_protein_substrate_cutoff',default = 2,type=int,help="define how many sequences share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-subs','--substrate_scors',default = 2,type=int,help="each cgc contains with substrate must more than this value")

    #subparser for CGC_substrate_prediction
    parser_easy_sub.add_argument('-input_sub_out',help='substrate out ')
    parser_easy_sub.add_argument('-b','--blastp',help='blastp result for cgc')
    parser_easy_sub.add_argument('--cgc', help='cgc_finder output')


    args = parser.parse_args()


    if args.command == 'database':
        run_dbCAN_database(args)

    if args.command == 'input_process':
        run_dbCAN_input_process(args)

    if args.command == 'CAZyme_annotation' and args.methods:
        run_dbCAN_CAZyme_annotation(args)

    if args.command == 'CGC_info':
        run_CGC_annotation_preprocess(args)


    if args.command == 'CGC_annotation':
        run_CGC_annotation(args)

    if args.command == 'CGC_substrate_prediction':
        args.input = args.output_dir
        cgc_substrate_prediction(args)

    if args.command == 'CGC_substrate_plot':
        args.input_sub_out = os.path.join(args.output_dir, 'substrate.out')
        args.blastp = os.path.join(args.output_dir, 'PUL_blast.out')
        args.cgc = os.path.join(args.output_dir, 'cgc_standard_out.tsv')
        syntenic_plot_allpairs(args)


    if args.command == 'easy_CAZyme':
        run_dbCAN_input_process(args)
        run_dbCAN_CAZyme_annotation(args)

    if args.command == 'easy_CGC':
        run_dbCAN_input_process(args)
        run_dbCAN_CAZyme_annotation(args)      
        run_CGC_annotation_preprocess(args)
        run_CGC_annotation(args)

    if args.command == 'easy_substrate':
        run_dbCAN_input_process(args)
        run_dbCAN_CAZyme_annotation(args)      
        run_CGC_annotation_preprocess(args)
        run_CGC_annotation(args)

        args.input = args.output_dir
        cgc_substrate_prediction(args)

        args.input_sub_out = os.path.join(args.output_dir, 'substrate.out')
        args.blastp = os.path.join(args.output_dir, 'PUL_blast.out')
        args.cgc = os.path.join(args.output_dir, 'cgc_standard_out.tsv')
        syntenic_plot_allpairs(args)


if __name__ == '__main__':
    main()