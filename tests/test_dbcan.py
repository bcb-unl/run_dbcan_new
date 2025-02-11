# run_dbcan testing file

import os
from pathlib import Path

from dbcan.cgc_substrate_prediction import cgc_substrate_prediction
from dbcan.run_dbcan import (
    run_CGC_annotation,
    run_CGC_annotation_preprocess,
    run_dbCAN_CAZyme_annotation,
    run_dbCAN_database,
    run_dbCAN_input_process,
)
from dbcan.syntenic_plot import syntenic_plot_allpairs

TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")


class Test_dbCAN:

    def test_run_dbCAN_database():
        args = type('', (), {})()
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        run_dbCAN_database(args)
        assert os.path.exists(os.path.join(args.db_dir, 'CAZy.dmnd'))
        assert os.path.exists(os.path.join(args.db_dir, 'dbCAN.hmm'))
        assert os.path.exists(os.path.join(args.db_dir, 'dbCAN_sub.hmm'))


    def test_easy_substrate_prok_fna():
        args = type('', (), {})()
        args.input_raw_data = os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.fna')
        args.mode = 'prok'
        args.output_dir = os.path.join(DATA_ROOT, 'test_dbCAN_prok_fna')
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        args.input_gff_format = 'prodigal'
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


        assert os.path.exists.join(args.output, 'uniInput.faa')
        assert os.path.exists.join(args.output, 'overview.tsv')
        assert os.path.exists.join(args.output, 'cgc_standard_out.tsv')
        assert os.path.exists.join(args.output, 'substrate.out')
        assert os.path.exists.join(args.output, 'synteny_pdf')

    def test_easy_substrate_prok_faa():
        args = type('', (), {})()
        args.input_raw_data = os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.faa')
        args.mode = 'protein'
        args.output_dir = os.path.join(DATA_ROOT, 'test_dbCAN_prok_faa')
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        args.input_format = "NCBI"
        args.input_gff_format = 'NCBI_prok'
        args.input_gff = os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.gff')

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


        assert os.path.exists.join(args.output, 'uniInput.faa')
        assert os.path.exists.join(args.output, 'overview.tsv')
        assert os.path.exists.join(args.output, 'cgc_standard_out.tsv')
        assert os.path.exists.join(args.output, 'substrate.out')
        assert os.path.exists.join(args.output, 'synteny_pdf')

    def test_easy_substrate_meta_fna():
        args = type('', (), {})()
        args.input_raw_data = os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.fna')
        args.mode = 'meta'
        args.output_dir = os.path.join(DATA_ROOT, 'test_dbCAN_meta_fna')
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        args.input_gff_format = 'prodigal'
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


        assert os.path.exists.join(args.output, 'uniInput.faa')
        assert os.path.exists.join(args.output, 'overview.tsv')
        assert os.path.exists.join(args.output, 'cgc_standard_out.tsv')
        assert os.path.exists.join(args.output, 'substrate.out')
        assert os.path.exists.join(args.output, 'synteny_pdf')


    def test_easy_substrate_NCBI_euk_faa():
        args = type('', (), {})()
        args.input_raw_data = os.path.join(DATA_ROOT, 'Xylona_heveae_TC161.faa')
        args.mode = 'protein'
        args.output_dir = os.path.join(DATA_ROOT, 'test_dbCAN_NCBI_euk_faa')
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        args.input_format = "NCBI"
        args.input_gff = os.path.join(DATA_ROOT, 'Xylona_heveae_TC161.gff')
        args.input_gff_format = 'NCBI_euk'

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


        assert os.path.exists.join(args.output, 'uniInput.faa')
        assert os.path.exists.join(args.output, 'overview.tsv')
        assert os.path.exists.join(args.output, 'cgc_standard_out.tsv')
        assert os.path.exists.join(args.output, 'substrate.out')
        assert os.path.exists.join(args.output, 'synteny_pdf')



    def test_easy_substrate_JGI_euk():

        args = type('', (), {})()
        args.input_raw_data = os.path.join(DATA_ROOT, 'Xylhe1_GeneCatalog_proteins_20130827.aa.fasta')
        args.mode = 'protein'
        args.output_dir = os.path.join(DATA_ROOT, 'test_dbCAN_JGI_euk_faa')
        args.db_dir = os.path.join(DATA_ROOT, 'test_dbCAN_database')
        args.input_format = "JGI"
        args.input_gff_format = 'JGI'
        args.input_gff = os.path.join(DATA_ROOT, 'Xylhe1_GeneCatalog_proteins_20130827.gff')

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


        assert os.path.exists.join(args.output, 'uniInput.faa')
        assert os.path.exists.join(args.output, 'overview.tsv')
        assert os.path.exists.join(args.output, 'cgc_standard_out.tsv')
        assert os.path.exists.join(args.output, 'substrate.out')
        assert os.path.exists.join(args.output, 'synteny_pdf')