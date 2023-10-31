#########The structure of the main script is modified from bin3C########
from Script.detect_marker_gene import detect_marker_gene
from Script.imputation import ImputeMatrix
from Script.exceptions import ApplicationException
from Script.pre_common import get_bin_dirs
from Script.pre_profile import Profile
from Script.utility import save_object, load_object, make_dir, gen_bins
from Script.final_cluster import FinalCluster
from Script.pre_clustering import Clust4CheckM, PreCluster

#######Import python packages
import subprocess
import argparse
import warnings
import logging
import shutil
import sys
import os
'''
'''
##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.1.0, released at 08/2023'

if __name__ == '__main__':
    
    def mk_version():
        return 'ImputeCC v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'gene_cov': 0.9,
        'rwr_rp': 0.5,
        'rwr_thres': 80,
        'max_markers': 8000,
        'intra': 50,
        'inter': 0,
        'min_binsize':100000,
        'cont_weight': 2,
        'min_comp': 50.0,
        'max_cont': 10.0,
        'report_quality': 10.0
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/MetaCC.log]')


    parser = argparse.ArgumentParser(description='ImputeCC: a metagenomic Hi-C-based binning method')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    #Impute the Hi-C matrix for contigs with marker genes
    cmd_impute = subparsers.add_parser('impute', parents=[global_parser],
                                      description='Imputation step.')
    
    #Clustering and ensemble                              
    cmd_cl = subparsers.add_parser('cluster', parents=[global_parser],
                                      description='Clustering step.')

    cmd_pipe = subparsers.add_parser('pipeline', parents=[global_parser],
                                      description='Run the whole ImputeCC pipeline.')
    cmd_test = subparsers.add_parser('test', parents=[global_parser],
                                      description='Test the whole ImputeCC pipeline using the test data.')
    
    '''
    Run the whole ImputeCC pipeline
    '''
    cmd_pipe.add_argument('--gene-cov', type=float, help='gene coverage used in detecting marker genes, default 0.9')
    cmd_pipe.add_argument('--rwr-rp', type=float, help='random walk restart probability, default 0.5')
    cmd_pipe.add_argument('--rwr-thres', type=int, help='cut-off to maintain sparsity in each random walk step, default 80')
    cmd_pipe.add_argument('--max-markers', type=int, help='maximum number of contigs with marker genes, default 8000')
    cmd_pipe.add_argument('--intra', type=int, help='percentile threshold to assign the contigs to preliminary bins in pre-clustering step, default 50')
    cmd_pipe.add_argument('--inter', type=int, help='percentile threshold to assign the contigs to new bins in pre-clustering step, default 0')
    cmd_pipe.add_argument('--cont-weight', type=float, help='coefficient of completeness - cont_weight * completeness, default 2')
    cmd_pipe.add_argument('--min-comp', type=float, help='minimum completeness of bin to consider during bin selection process, default 50')
    cmd_pipe.add_argument('--max-cont' , type=float, help='maximum contamination of bin to consider during bin selection process, default 10')
    cmd_pipe.add_argument('--report-quality', type=float, help="minimum quality of bin to report, default 10")
    cmd_pipe.add_argument('--min-binsize', type=int, help='Minimum bin size used in output [100000]')
    cmd_pipe.add_argument('-t', '--threads', default=30, type=int, help="the number of threads. default is 30.")
    
    cmd_pipe.add_argument('FASTA', help='Reference fasta sequence')
    cmd_pipe.add_argument('CONTIG_INFO', help='contig information file in csv format') # generate by NormCC
    cmd_pipe.add_argument('HIC_MATRIX', help='normalized Hi-C contact matrix in npz format') # generate by NormCC
    cmd_pipe.add_argument('OUTDIR', help='Output directory of sub bins')
    
    '''
    Imputation 
    '''
    #Parameter
    cmd_impute.add_argument('--gene-cov', type=float, help='gene coverage used in detecting marker genes, default 0.9')
    cmd_impute.add_argument('--rwr-rp', type=float, help='random walk restart probability, default 0.5')
    cmd_impute.add_argument('--rwr-thres', type=int, help='cut-off to maintain sparsity in each random walk step, default 80')
    cmd_impute.add_argument('--max-markers', type=int, help='maximum number of contigs with marker genes, default 8000')
    cmd_impute.add_argument('-t', '--threads', default=30, type=int, help="the number of threads. default is 30.")
    #Input files
    cmd_impute.add_argument('FASTA', help='contig sequence file in fasta format')
    cmd_impute.add_argument('CONTIG_INFO', help='contig information file in csv format') # generate by NormCC
    cmd_impute.add_argument('HIC_MATRIX', help='normalized Hi-C contact matrix in npz format') # generate by NormCC
    cmd_impute.add_argument('OUTDIR', help='Output directory')
        

    '''
    Clutering subsparser input
    '''
    cmd_cl.add_argument('--intra', type=int, help='percentile threshold to assign the contigs to preliminary bins in pre-clustering step, default 50')
    cmd_cl.add_argument('--inter', type=int, help='percentile threshold to assign the contigs to new bins in pre-clustering step, default 0')
    cmd_cl.add_argument('--cont-weight', type=float, help='coefficient of completeness - cont_weight * completeness, default 2')
    cmd_cl.add_argument('--min-comp', type=float, help='minimum completeness of bin to consider during bin selection process, default 50')
    cmd_cl.add_argument('--max-cont' , type=float, help='maximum contamination of bin to consider during bin selection process, default 10')
    cmd_cl.add_argument('--report-quality', type=float, help="minimum quality of bin to report, default 10")
    cmd_cl.add_argument('--min-binsize', type=int, help='Minimum bin size used in output [100000]')
    cmd_cl.add_argument('-t', '--threads', default=30, type=int, help="the number of threads. default is 30.")
    #cmd_cl.add_argument('--seed', type=int, default=None,
    #                           help='Random seed')
    cmd_cl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_cl.add_argument('CONTIG_INFO', help='contig information file in csv format') # generate by NormCC
    cmd_cl.add_argument('HIC_MATRIX', help='normalized Hi-C contact matrix in npz format') # generate by NormCC
    cmd_cl.add_argument('OUTDIR', help='Output directory of sub bins')
    
    
    '''
    Testing of NormCC software
    
    cmd_test.add_argument('--OUTDIR', type=str, default='Test/out_test', help='Output directory of testing results')
    '''
    
    args = parser.parse_args()

    if args.command == 'test':
        imputecc_dir = os.path.abspath(os.path.dirname(__file__))
        script_path = os.path.join(imputecc_dir, 'ImputeCC.py')
        contig_path = os.path.join(imputecc_dir, 'test_data/test_contigs.fa')
        info_path = os.path.join(imputecc_dir, 'test_data/test_info.csv')
        matrix_path = os.path.join(imputecc_dir, 'test_data/test_normalized_contact_matrix.npz')
        output_dir = os.path.join(imputecc_dir, 'test_result/')
        cmd = ' '.join(['python', script_path, 'pipeline -t 60', contig_path, info_path,  matrix_path, output_dir, '--cover -v'])
        print('The pipeline was tested by the command: ' + cmd)
        subprocess.run(cmd, shell=True)
        shutil.rmtree(output_dir)
        print('Pipeline test successfully finished...')
        sys.exit(0)

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    
    # Create temp folder
    # Temp folder is deleted for each step
    temp_folder = os.path.join(args.OUTDIR , 'tmp')
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    else:
        shutil.rmtree(temp_folder)           
        os.mkdir(temp_folder)
    
    # Create Intermediate folder
    # Intermediate folder is not deleted
    interm_folder = os.path.join(args.OUTDIR , 'intermediate')
    if not os.path.exists(interm_folder):
        os.mkdir(interm_folder)
           
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'ImputeCC.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        if args.command == 'pipeline':
            logger.info('The ImputeCC binning pipeline starts...')
            marker_file = os.path.join(interm_folder , 'contigs_marker_genes.hmmout')
            if os.path.exists(marker_file):
                logger.info('Existing single-copy marker gene file is detected.')
            else:
                logger.info('Begin detecting single-copy marker genes from assembled contigs...')
                detect_marker_gene(args.FASTA, args.threads, interm_folder)
            
            if not os.path.exists(marker_file):
                raise ModuleNotFoundError('Single-copy marker gene detection failed!')
            else:
                logger.info('Single-copy marker gene detection finished.')
            
            logger.info('Begin CRWR Imputation...')
            imp = ImputeMatrix(args.CONTIG_INFO,
                                args.HIC_MATRIX,
                                marker_file,
                                gene_cov = ifelse(args.gene_cov, runtime_defaults['gene_cov']),
                                rwr_rp = ifelse(args.rwr_rp, runtime_defaults['rwr_rp']),
                                rwr_thres= ifelse(args.rwr_thres, runtime_defaults['rwr_thres']),
                                max_markers= ifelse(args.max_markers, runtime_defaults['max_markers'])                                                  
                                )
            
            save_object(os.path.join(interm_folder, 'ImputeCC_storage'), imp)
            logger.info('CRWR Imputation finished.')
            
            #######Construct preliminary bins#####
            logger.info('Begin preclustering marker-gene-containing contigs...')
            bins, bin_of_contigs = PreCluster(imp.marker_contig_counts, imp.marker_contigs, imp.contig_markers,
                                             imp.imputed_matrix, imp.dict_contigRevLocal,
                                             intra= ifelse(args.intra , runtime_defaults['intra']),
                                             inter = ifelse(args.inter , runtime_defaults['inter']))
            logger.info('Preclustering finished with {} preliminary bins established.'.format(len(bins)))
            
            
            logger.info('Begin final binning for all contigs utilizing the information of preliminary bins...')            
            checkm_bac_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_bac_gene_table.tsv')
            checkm_ar_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_ar_gene_table.tsv')
            if os.path.exists(checkm_ar_gene_table) and os.path.exists(checkm_bac_gene_table):
                logger.info('CheckM lineage-specific gene tables detected.')
            else:
                logger.info('Detect lineage-specific genes...')
                Clust4CheckM(args.FASTA, args.CONTIG_INFO, args.HIC_MATRIX, args.OUTDIR)
                output_dir = os.path.join(temp_folder, 'out_checkm_qa')
                bin_file = os.path.join(temp_folder , 'dir4checkm.tsv')
                cpus = args.threads
                bin_dirs = get_bin_dirs(bin_file)
                profile = Profile(cpus)
                profile.run(bin_dirs,
                            output_dir) 
                os.mkdir(os.path.join(interm_folder, 'checkm_gene_table'))

                mv1Cmd = 'mv ' + os.path.join(temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                            'checkm_bac', 'marker_gene_table.tsv') + ' ' + checkm_bac_gene_table
                
                mv2Cmd = 'mv ' + os.path.join(temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                            'checkm_ar', 'marker_gene_table.tsv') + ' ' + checkm_ar_gene_table
                os.system(mv1Cmd)
                os.system(mv2Cmd)
                logger.info('Lineage-specific gene detection finished.')
                
            cluster_process = FinalCluster(imp.contig_info, imp.contig_local, 
                 imp.dict_contigRev, imp.dict_contigRevLocal, imp.dict_contig_len, 
                 imp.contig_markers, bins, bin_of_contigs,
                 imp.normcc_matrix, imp.imputed_matrix, 
                 checkm_bac_gene_table, checkm_ar_gene_table, args.OUTDIR,
                 intra= ifelse(args.intra , runtime_defaults['intra']),
                 inter = ifelse(args.inter , runtime_defaults['inter']),
                 cont_weight = ifelse(args.cont_weight , runtime_defaults['cont_weight']),
                 min_comp = ifelse(args.min_comp , runtime_defaults['min_comp']), 
                 max_cont = ifelse(args.max_cont , runtime_defaults['max_cont']), 
                 report_quality = ifelse(args.report_quality , runtime_defaults['report_quality']),
                 min_binsize = ifelse(args.min_binsize , runtime_defaults['min_binsize']))
            
            logger.info('Final binning finished!')
            logger.info('Writing the final bins...')
            gen_bins(args.FASTA , os.path.join(temp_folder , 'cluster_imputecc.txt') , os.path.join(args.OUTDIR ,'FINAL_BIN'))
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            logger.info('The ImputeCC binning pipeline fininshed.')
                
            
            
        if args.command == 'impute':
            logger.info('Run the CRWR imputation step...')
            marker_file = os.path.join(interm_folder , 'contigs_marker_genes.hmmout')
            if os.path.exists(marker_file):
                logger.info('Existing single-copy marker gene file is detected.')
            else:
                logger.info('Begin detecting single-copy marker genes from assembled contigs...')
                detect_marker_gene(args.FASTA, args.threads, interm_folder)
            
            if not os.path.exists(marker_file):
                raise ModuleNotFoundError('Single-copy marker gene detection failed!')
            else:
                logger.info('Single-copy marker gene detection finished.')
            
            logger.info('Begin CRWR Imputation...')
            imp = ImputeMatrix(args.CONTIG_INFO,
                                args.HIC_MATRIX,
                                marker_file,
                                gene_cov = ifelse(args.gene_cov, runtime_defaults['gene_cov']),
                                rwr_rp = ifelse(args.rwr_rp, runtime_defaults['rwr_rp']),
                                rwr_thres= ifelse(args.rwr_thres, runtime_defaults['rwr_thres']),
                                max_markers= ifelse(args.max_markers, runtime_defaults['max_markers'])                                                  
                                )
            
            save_object(os.path.join(interm_folder, 'ImputeCC_storage'), imp)
            logger.info('CRWR Imputation finished.')
            shutil.rmtree(temp_folder) ######Remove temp folder, intermediate folder should not be deleted#######
            
           
           

        if args.command == 'cluster':
            logger.info('Run the clusterings step...')
            if not os.path.exists(os.path.join(interm_folder , 'ImputeCC_storage.gz')):
                raise IOError('Please run the Imputation step before clustering step')
                 
            ###########Load the normalization instance to get access to the normalized Hi-C contact maps##########
            logger.info('Loading imputed contact maps by ImputeCC from: {}'.format(os.path.join(interm_folder, 'ImputeCC_storage.gz')))
            imp = load_object(os.path.join(interm_folder , 'ImputeCC_storage.gz'))
            
            #######Construct preliminary bins#####
            logger.info('Begin preclustering marker-gene-containing contigs...')
            bins, bin_of_contigs = PreCluster(imp.marker_contig_counts, imp.marker_contigs, imp.contig_markers,
                                             imp.imputed_matrix, imp.dict_contigRevLocal,
                                             intra= ifelse(args.intra , runtime_defaults['intra']),
                                             inter = ifelse(args.inter , runtime_defaults['inter']))
            logger.info('Preclustering finished with {} preliminary bins established.'.format(len(bins)))
            
            
            logger.info('Begin final binning of all contigs utilizing the information of preliminary bins...')            
            checkm_bac_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_bac_gene_table.tsv')
            checkm_ar_gene_table = os.path.join(interm_folder, 'checkm_gene_table' , 'checkm_ar_gene_table.tsv')
            if os.path.exists(checkm_ar_gene_table) and os.path.exists(checkm_bac_gene_table):
                logger.info('CheckM lineage-specific gene tables detected.')
            else:
                logger.info('Detect CheckM lineage-specific genes...')
                Clust4CheckM(args.FASTA, args.CONTIG_INFO, args.HIC_MATRIX, args.OUTDIR)
                output_dir = os.path.join(temp_folder, 'out_checkm_qa')
                bin_file = os.path.join(temp_folder , 'dir4checkm.tsv')
                cpus = args.threads
                bin_dirs = get_bin_dirs(bin_file)
                profile = Profile(cpus)
                profile.run(bin_dirs,
                            output_dir) 
                os.mkdir(os.path.join(interm_folder, 'checkm_gene_table'))

                mv1Cmd = 'mv ' + os.path.join(temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                            'checkm_bac', 'marker_gene_table.tsv') + ' ' + checkm_bac_gene_table
                
                mv2Cmd = 'mv ' + os.path.join(temp_folder, 'out_checkm_qa', 'binning_methods', 'INITIAL_BIN', 
                                            'checkm_ar', 'marker_gene_table.tsv') + ' ' + checkm_ar_gene_table
                os.system(mv1Cmd)
                os.system(mv2Cmd)
                logger.info('CheckM lineage-specific genes detection finished.')
                
            
            #if not args.seed:
            #    args.seed = make_random_seed()   
            #logger.info('The random seed for clustering is {}'.format(args.seed))
            
            cluster_process = FinalCluster(imp.contig_info, imp.contig_local, 
                 imp.dict_contigRev, imp.dict_contigRevLocal, imp.dict_contig_len, 
                 imp.contig_markers, bins, bin_of_contigs,
                 imp.normcc_matrix, imp.imputed_matrix, 
                 checkm_bac_gene_table, checkm_ar_gene_table, args.OUTDIR,
                 intra= ifelse(args.intra , runtime_defaults['intra']),
                 inter = ifelse(args.inter , runtime_defaults['inter']),
                 cont_weight = ifelse(args.cont_weight , runtime_defaults['cont_weight']),
                 min_comp = ifelse(args.min_comp , runtime_defaults['min_comp']), 
                 max_cont = ifelse(args.max_cont , runtime_defaults['max_cont']), 
                 report_quality = ifelse(args.report_quality , runtime_defaults['report_quality']),
                 min_binsize = ifelse(args.min_binsize , runtime_defaults['min_binsize']))
            
            logger.info('Final binning finished.')
            logger.info('Writing the final bins...')
            gen_bins(args.FASTA , os.path.join(temp_folder , 'cluster_imputecc.txt') , os.path.join(args.OUTDIR ,'FINAL_BIN'))
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            logger.info('The ImputeCC clustering step fininshed.')
   

        '''  
        if args.command == 'test':
            logger.info('Begin to test MetaCC...')
            ENZ = 'HindIII'
            FASTA = 'Test/final.contigs.fa'
            BAM = 'Test/MAP_SORTED.bam'
            OUT = args.OUTDIR
            logger.info('Begin to test the contact map construction section...')
            
            shutil.rmtree(OUT, ignore_errors=True)
            logger.info('Testing finished!')
        ''' 


    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
        
        
