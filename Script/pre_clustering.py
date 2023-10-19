import igraph as ig
import leidenalg
import os
import scipy.sparse as scisp
import pandas as pd
import numpy as np
import operator
from Script.utility import gen_bins, match_contigs

def PreCluster(marker_contig_counts, marker_contigs, contig_markers, imputed_matrix, dict_contigRevLocal, intra, inter):
        _my_gene_counts = list(marker_contig_counts.values())
        _my_gene_counts.sort(reverse=True)      

        smg_iteration = {}
        n = 0
        _unique_my_gene_counts = list(set(_my_gene_counts))
        _unique_my_gene_counts.sort(reverse=True)

        # Get contigs for each iteration of single-copy marker gene
        for _g_count in _unique_my_gene_counts:
            # Get the single-copy marker genes with maximum count of contigs and
            # sort them in the descending order of the total marker genes contained

            total_contig_mgs = {}
            
            ####item here is the single copy marker genes
            #If multiple genes correspond to the same number of contigs, we will sort the genes according to the number of genes on contigs
            for item in marker_contig_counts:
                if marker_contig_counts[item] == _g_count:
                    total_contig_lengths = 0
                    
                    #Marker contigs are dict gene: contigs###
                    for contig in marker_contigs[item]:
                        contig_mg_counts = len(contig_markers[contig])
                        total_contig_lengths += contig_mg_counts

                    total_contig_mgs[item] = total_contig_lengths

            total_contig_mgs_sorted = sorted(
                total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True
            )

            for item in total_contig_mgs_sorted:
                smg_iteration[n] = marker_contigs[item[0]]
                n += 1
                
        bins, bin_of_contig, _, _, _ = match_contigs(
            smg_iteration,
            contig_markers,
            imputed_matrix.tolil(),
            dict_contigRevLocal,
            np.percentile(imputed_matrix.tocoo().data , intra),
            np.percentile(imputed_matrix.tocoo().data , inter)
        )
        return bins, bin_of_contig



def Clust4CheckM(fasta_file, contig_info_file, normcc_matrix_file, path):
    _map_del = scisp.load_npz(normcc_matrix_file).tocoo()
    contig_info = pd.read_csv(contig_info_file , sep = ',' , header = 0).values
    
    #########Use Leiden Algorithm to do clustering########
    _vcount = _map_del.shape[0]
    _sources = _map_del.row
    _targets = _map_del.col
    _wei = _map_del.data
    _index = _sources<_targets
    _sources = _sources[_index]
    _targets = _targets[_index]
    _wei = _wei[_index]
    _edgelist = list(zip(_sources, _targets))
    g = ig.Graph(_vcount, _edgelist)
    
    #############determine the resolution parameter###########
    part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=_wei)
    part = list(part)
    _dict_cluster = {}
    # dict of communities

    for ci in range(len(part)):
        for id in part[ci]:
            _dict_cluster[contig_info[id , 0]] = 'group'+str(ci)
            
    with open(os.path.join(path , 'tmp' , 'cluster4checkm.txt'),'w') as out:
        for key , value in _dict_cluster.items():
            out.write(str(key)+ '\t' +str(value)+ '\n')
            
    with open(os.path.join(path , 'tmp' , 'dir4checkm.tsv'),'w') as out:
        out.write('INITIAL_BIN' + '\t' + os.path.join(path , 'tmp' ,'BIN4checkm'))
    
    gen_bins(fasta_file , os.path.join(path , 'tmp' , 'cluster4checkm.txt') , os.path.join(path , 'tmp' ,'BIN4checkm'))
    
    
            
