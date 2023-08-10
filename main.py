#########The structure of the main script is modified from bin3C########
from imputecc import ImputeMatrix
from utility import save_object, load_object
from utils import make_dir, gen_bins
from exceptions import ApplicationException
from biolib.common import make_sure_path_exists
from unitem_common import get_bin_dirs
from unitem_profile import Profile
import argparse
import warnings
import logging
import shutil
import sys
import os

#import scipy.sparse as scisp

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
        #'min_len': 1000,
        #'min_signal': 2,
        #'min_mapq': 30,
        #'min_match': 30,
        #'thres': 0.05,
        'gene_cov': 0.7,
        'rwr_rp': 0.3,
        'rwr_thres': 70,
        'intra': 50,
        'inter': 5,
        'min_binsize':150000
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/MetaCC.log]')


    parser = argparse.ArgumentParser(description='MetaCC: a scalable and integrative analysis framework for both short-read and long-read metagenomic Hi-C datasets')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    #Impute the Hi-C matrix for contigs with marker genes
    cmd_impute = subparsers.add_parser('impute', parents=[global_parser],
                                      description='Normalize contacts.')
    
    cmd_unitem = subparsers.add_parser('unitem', parents=[global_parser],
                                      description='unitem.')
    #Clustering and ensemble                              
    cmd_cl = subparsers.add_parser('cluster', parents=[global_parser],
                                      description='Do the binning.')

    cmd_ensemble = subparsers.add_parser('ensemble', parents=[global_parser],
                                        description='ensemble.')

    '''
    Imputation subparser input
    '''
    #cmd_norm.add_argument('--min-len', type=int,
    #                       help='Minimum acceptable contig length [1000]')
    #cmd_norm.add_argument('--min-signal', type=int,
    #                       help='Minimum acceptable Hi-C signal [2]')
    #cmd_norm.add_argument('--min-mapq', type=int,
    #                       help='Minimum acceptable mapping quality [30]')
    #cmd_norm.add_argument('--min-match', type=int,
    #                       help='Accepted alignments must being N matches [30]')
    #cmd_norm.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
    #                       help='Case-sensitive enzyme name. Use multiple times for multiple enzymes')
    #cmd_norm.add_argument('--thres', type=float,
    #                       help='the fraction of discarded NormCC-normalized Hi-C contacts [0.05]')
    #Parameter
    cmd_impute.add_argument('--gene-cov', type=float, help='gene coverage used in detecting marker genes')
    cmd_impute.add_argument('--rwr-rp', type=float, help='random walk restart probability')
    cmd_impute.add_argument('--rwr-thres', type=int, help='cut-off to maintain sparsity in each random walk step')
    cmd_impute.add_argument('--inter', type=int, help='inter-cluster threshold in pre-clustering step')
    cmd_impute.add_argument('--intra', type=int, help='intra-cluster threshold in pre-clustering step')
    
    #Input files
    cmd_impute.add_argument('FASTA', help='final.contigs.fa')
    cmd_impute.add_argument('CONTIG', help='contig_info.csv')
    cmd_impute.add_argument('MARKER', help='marker.hmmout')
    cmd_impute.add_argument('MATRIX', help='hic_matrix.npz')
    cmd_impute.add_argument('OUTDIR', help='Output directory')
    
    '''
    Unitem
    '''
    cmd_unitem.add_argument('--threads', default=20, type=int, help="the number of threads. default is 20.")
    cmd_unitem.add_argument('OUTDIR', help='Output directory for ImputeCC')
    
    
    '''
    Clutering subsparser input
    '''
    cmd_cl.add_argument('--min-binsize', type=int,
                               help='Minimum bin size used in output [150000]')
    cmd_cl.add_argument('--marker-gene', type=str,
                               help='Number of maker genes detected, automatically detected if not input')
    cmd_cl.add_argument('--seed', type=int, default=None,
                               help='Random seed')
    cmd_cl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_cl.add_argument('OUTDIR', help='Output directory of sub bins')
    

    '''
    Testing of NormCC software
    
    cmd_test.add_argument('--OUTDIR', type=str, default='Test/out_test', help='Output directory of testing results')
    '''
    
    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    
    # Create temp folder
    temp_folder = os.path.join(args.OUTDIR , 'tmp')
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    #else:
    #    shutil.rmtree(temp_folder)           
    #    os.mkdir(temp_folder)
           
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
        if args.command == 'impute':
            logger.info('Begin Imputation')
            imp = ImputeMatrix(args.CONTIG,
                                args.MARKER,
                                args.MATRIX,
                                args.OUTDIR,
                                gene_cov = ifelse(args.gene_cov, runtime_defaults['gene_cov']),
                                rwr_rp = ifelse(args.rwr_rp, runtime_defaults['rwr_rp']),
                                rwr_thres= ifelse(args.rwr_thres, runtime_defaults['rwr_thres']),
                                intra= ifelse(args.intra , runtime_defaults['intra']),
                                inter= ifelse(args.inter, runtime_defaults['inter'])                                                         
                                )
            
            save_object(os.path.join(args.OUTDIR, 'ImputeCC_storage'), imp)
            gen_bins(args.FASTA , os.path.join(temp_folder , 'cluster4unitem.txt') , os.path.join(temp_folder ,'BIN4unitem'))
            logger.info('Imputation finished!')
            
        
        if args.command == 'unitem':
            logger.info('run unitem...')
            
            output_dir = os.path.join(temp_folder, 'out_unitem')
            bin_file = os.path.join(temp_folder , 'dir4unitem.tsv')
            cpus = args.threads

            bin_dirs = get_bin_dirs(bin_file)
            make_sure_path_exists(output_dir)

            profile = Profile(cpus)
            profile.run(bin_dirs,
                        output_dir) 
           
           
        ''' 
        if args.command == 'cluster':
            #if not os.path.exists(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')):
            #    raise IOError('Please run the NormCC normalization step before binning')
            
            ###########Load the normalization instance to get access to the normalized Hi-C contact maps##########
            logger.info('Loading imputed contact maps by ImputeCC from: {}'.format(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')))
            hzmap = load_object(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz'))
            
            #########Scan the marker gene to determine the hyperparameter in the Leiden clustering#########
            if args.marker_gene is None:
                logger.info('Begin scanning marker genes...')
                
                marker_gene_folder = os.path.join(args.OUTDIR , 'marker_gene')
                if not os.path.exists(marker_gene_folder):
                    os.mkdir(marker_gene_folder)
        
                num_gene = gen_bestk(args.OUTDIR , args.FASTA)
                args.marker_gene = os.path.join(marker_gene_folder, 'contigs.hmmout')
                logger.warning('Finish detecting marker genes from the assembled contigs!')
            
            if not args.seed:
                args.seed = make_random_seed()
                
            logger.info('The random seed for clustering is {}'.format(args.seed))
                
            cluster_process = ClusterBin(args.OUTDIR , hzmap.name , hzmap.len , hzmap.seq_map ,
                                            ifelse(args.min_binsize, runtime_defaults['min_binsize']), args.marker_gene, args.seed)
            logger.info('Writing bins...')
            gen_bins(args.FASTA , os.path.join(temp_folder , 'cluster.txt') , os.path.join(args.OUTDIR ,'BIN'))
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            logger.info('ImputeCC binning fininshed.')
        '''     

        '''  
        if args.command == 'test':
            logger.info('Begin to test MetaCC...')
            ENZ = 'HindIII'
            FASTA = 'Test/final.contigs.fa'
            BAM = 'Test/MAP_SORTED.bam'
            OUT = args.OUTDIR
            logger.info('Begin to test the contact map construction section...')
            cm = ContactMatrix(BAM,
                                ENZ,
                                FASTA,
                                OUT,
                                min_mapq=runtime_defaults['min_mapq'],
                                min_len=runtime_defaults['min_len'],
                                min_match=runtime_defaults['min_match'],
                                min_signal=0)

            logger.info('Contact map construction section works!')
            
            logger.info('Begin to test the NormCC normalization module...')
            logger.info('Normalizing raw contacts by NormCC...')
            
            from rpy2 import robjects
            r = robjects.r
            r.source('NormCC/normcc.R')
            contig_file = 'Test/contig_info_test.csv'
            norm_result = r.normcc(contig_file)
            
            ######Construct normalized matrix of Hi-C interaction maps#############
            hzmap = NormCCMap(OUT,
                            cm.seq_info,
                            cm.seq_map,
                            norm_result,
                            thres = runtime_defaults['thres'])
                            
            logger.info('NormCC Normalization module works!')


            logger.info('Begin to test the MetaCC binning module...')
            logger.info('Begin scanning marker genes...')
            logger.info('Leiden clustering starts...')
            cluster_process = ClusterBin(OUT , hzmap.name , hzmap.len , hzmap.seq_map ,
                                            0, 0, 0)
                                            
            logger.info('MetaCC binning module works!')
            
            shutil.rmtree(OUT, ignore_errors=True)
            logger.info('Testing finished!')
        ''' 


    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
        
        

        '''
        Semi-supervised Leiden
        
        is_membership_fixed = []
        new_membership = []
        index = len(bins)+1

        for i in contig_info[: , 0]:
            if i in bin_of_contig:
                new_membership.append(bin_of_contig[i])
                is_membership_fixed.append(True)
            else:
                new_membership.append(index)
                index += 1
                is_membership_fixed.append(False)
                
                
        map_lil = normcc_matrix
        contig_temp = contig_info

        map_del = map_lil.tocoo()
        vcount = map_del.shape[0]
        sources = map_del.row
        targets = map_del.col
        wei = map_del.data
        index = sources<targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        edgelist = list(zip(sources, targets))
        g = ig.Graph(vcount, edgelist)
        
        for rp in reso_para_list:

            new_partition = leidenalg.RBConfigurationVertexPartition(g, new_membership, weights=wei , resolution_parameter = rp)

            #new_partition = leidenalg.CPMVertexPartition(g, new_membership, weights=wei)

            optimiser = leidenalg.Optimiser()
            diff = optimiser.optimise_partition(new_partition, is_membership_fixed=is_membership_fixed, n_iterations = -1)
            mem = new_partition.membership
            
            part = [[] for i in range(len(set(mem)))]
            for i, j in enumerate(mem):
                part[j].append(i)
                
            dist_cluster = {}
            numnode = 0
            rang = []
            for ci in range(len(part)):
                if sum(contig_temp[part[ci] , 2]) >= 150000:
                    rang.append(ci)
                    numnode = numnode+len(part[ci])
                    for id in part[ci]:
                        dist_cluster[contig_info[id, 0]] = 'group' + str(ci)
                    
            with open('cluster_'+ str(rwr_rp) + '_' + str(rwr_perc) + '_' + str(w_intra) + '_' + str(w_inter) + '_' + str(rp) +  '.txt','w') as out:
                for key , value in dist_cluster.items():
                    out.write(str(key)+ '\t' +str(value))
                    out.write('\n')
        '''
