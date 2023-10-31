import os
import sys
import warnings
import logging
##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
# package logger
logger = logging.getLogger(__name__)

def detect_marker_gene(contig_file, threads, output_dir):
    # Execution of the software
    fragScanURL = 'run_FragGeneScan.pl'
    hmmExeURL = 'hmmsearch'
    markerURL = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'marker.hmm')
    
    # Output file for fraggenescan and hmmer 
    fragResultURL = os.path.join(output_dir , 'contigs.frag.faa')
    hmmResultURL = os.path.join(output_dir, 'contigs_marker_genes.hmmout')

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + os.path.join(output_dir , 'contigs.frag') + " -complete=0 -train=complete -thread=" + str(
            threads) + " 1>" + os.path.join(output_dir , 'contigs.frag.out') + " 2>" + os.path.join(output_dir , 'contigs.frag.err')
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu " + str(
                threads) + " " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            os.system(hmmCmd)
            
            if not (os.path.exists(hmmResultURL)):
                logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
                sys.exit()
            
            rm1Cmd = 'rm' + ' -f ' + os.path.join(output_dir , 'contigs.frag.*') 
            os.system(rm1Cmd)
            
            rm2Cmd = 'rm' + ' -f ' + os.path.join(output_dir , 'contigs_marker_genes.hmmout.*')
            os.system(rm2Cmd)

    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()
        

