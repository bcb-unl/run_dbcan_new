import rich_click as click
from parameter import (
    create_config, GeneralConfig, DBDownloaderConfig, DiamondConfig, DiamondTCConfig,
    PyHMMERConfig, DBCANSUBProcessorConfig, PyHMMERTFConfig, PyHMMERSTPConfig,OverviewGeneratorConfig, GFFConfig, CGCFinderConfig, CGCSubstrateConfig, 
    SynPlotConfig, CGCPlotConfig,
    general_options, database_options, output_dir_option, methods_option, threads_option, diamond_options, diamond_tc_options, 
    pyhmmer_dbcan_options, dbcansub_options ,pyhmmer_tf, pyhmmer_stp, cgc_gff_option, cgc_options, cgc_sub_options, syn_plot_options, 
    cgc_circle_plot_options, cgc_substrate_base_options, cgc_substrate_homology_params_options, cgc_substrate_dbcan_sub_param_options
)
from core import (
    run_dbCAN_database, run_dbCAN_input_process, run_dbCAN_CAZyme_annotation,
    run_dbCAN_CGCFinder_preprocess, run_dbCAN_CGCFinder,
    run_dbCAN_CGCFinder_substrate, run_dbcan_syn_plot, run_dbCAN_cgc_circle
)

@click.group()
def cli():
    """use dbCAN tools to annotate and analyze CAZymes and CGCs."""
    pass

@cli.command('database')
@database_options
@click.pass_context
def database_cmd(ctx, **kwargs):
    """download dbCAN databases."""
    config = create_config(DBDownloaderConfig, **kwargs)
    run_dbCAN_database(config)

@cli.command('CAZyme_annotation')
@general_options
@database_options
@output_dir_option
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@click.pass_context
def cazyme_annotation_cmd(ctx, **kwargs):
    """process input data."""
    config = create_config(GeneralConfig, **kwargs)
    run_dbCAN_input_process(config)
    diamond_config = create_config(DiamondConfig,  **kwargs)
    pyhmmer_config = create_config(PyHMMERConfig, **kwargs)
    dbcansubconfig = create_config(DBCANSUBProcessorConfig, **kwargs)
    overviewconfig = create_config(OverviewGeneratorConfig, **kwargs)
    methods_option = kwargs.get('methods')
    run_dbCAN_CAZyme_annotation(diamond_config, pyhmmer_config, dbcansubconfig, overviewconfig, methods_option)




@cli.command('gff_process')
@database_options
@output_dir_option
@threads_option
@pyhmmer_stp
@pyhmmer_tf
@diamond_tc_options
@cgc_gff_option
@click.pass_context
def gff_process_cmd(ctx, **kwargs):
    diamond_tc_config = create_config(DiamondTCConfig, **kwargs)
    pyhmmer_tf_config = create_config(PyHMMERTFConfig, **kwargs)
    pyhmmer_stp_config = create_config(PyHMMERSTPConfig, **kwargs)
    gff_config = create_config(GFFConfig, **kwargs)
    run_dbCAN_CGCFinder_preprocess(diamond_tc_config, pyhmmer_tf_config, pyhmmer_stp_config, gff_config)



@cli.command('cgc_finder')
@cgc_options
@click.pass_context
def cgc_finder_cmd(ctx, **kwargs):
    """identify CAZyme Gene Clusters(CGCs)"""
    config = create_config(CGCFinderConfig, **kwargs)
    run_dbCAN_CGCFinder(config)




@cli.command('substrate_prediction')
@cgc_substrate_base_options
@cgc_substrate_homology_params_options
@cgc_substrate_dbcan_sub_param_options
@click.pass_context
def substrate_prediction_cmd(ctx, **kwargs):
    """predict substrate specificities of CAZyme Gene Clusters(CGCs)."""
    config = create_config(CGCSubstrateConfig, **kwargs)
    run_dbCAN_CGCFinder_substrate(config)


@cli.command('syntenic_plot')
@syn_plot_options
@click.pass_context
def syntenic_plot_cmd(ctx, **kwargs):
    """generate syntenic plots for CAZyme Gene Clusters(CGCs)."""
    config = create_config(SynPlotConfig, **kwargs)
    run_dbcan_syn_plot(config)

@cli.command('cgc_circle_plot')
@cgc_circle_plot_options
@click.pass_context
def cgc_circle_plot_cmd(ctx, **kwargs):
    """generate circular plots for CAZyme Gene Clusters(CGCs)."""
    config = create_config(CGCPlotConfig, **kwargs)
    run_dbCAN_cgc_circle(config)

if __name__ == "__main__":
    cli()