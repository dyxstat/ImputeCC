#!/usr/bin/env python
# coding: utf-8

import logging
import warnings
import numpy as np
import pandas as pd
from utility import get_contigs_with_marker_genes
import operator
import scipy.sparse as scisp
import igraph as ig
import leidenalg
from utility import random_walk_cpu, get_contigs_with_marker_genes, match_contigs
import os


##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
# package logger
logger = logging.getLogger(__name__)


'''
Cut off selection
'''
class ImputeMatrix:
    def __init__(self, contig_file, marker_file, matrix_file, path, gene_cov, rwr_rp, rwr_thres, intra, inter):
        self.thres_mg = gene_cov
        self.rwr_rp = rwr_rp
        self.rwr_perc = rwr_thres
        self.w_intra = intra
        self.w_inter = inter
        self.path = path

        '''
        Information of assembled contigs
        dict_contig_len:    dict[contig_name] = contig_length
        dict_contigRev:     dict[contig_name] = row number in contig info np array
        '''
        ######Input the contig information
        #####input the blastn result for the contigs on reference genome###########

        self.contig_info = pd.read_csv(contig_file , sep = ',' , header = 0).values

        self.dict_contig_len = {}
        for i in self.contig_info:
            self.dict_contig_len[i[0]] = i[2]
            
        self.dict_contigRev = {} #dict[contig_name] = row number in contig info np array
        for i , j in enumerate(self.contig_info[: , 0]):
            self.dict_contigRev[j] = i
            
        ######Input normalized Hi-C matrix#######
        self.normcc_matrix = scisp.load_npz(matrix_file)
        self.normcc_matrix.tolil()
        
        '''
        Detect Marker genes
        '''
        ####Input marker gene information
        self.marker_contigs, self.marker_contig_counts, self.contig_markers = get_contigs_with_marker_genes(marker_file , self.thres_mg, self.dict_contigRev)

        self.contig_local, self.dict_contigRevLocal, self.imputed_matrix = self._imputation()
        
        self.bins, self.bin_of_contigs = self._pre_clustering()
        self._unitem_clustering()

    def  _imputation(self):
        '''
        Impute normcc-normalized Hi-C contact matrix using randow walk with restart
        '''
        ###Only impute part of matrix instead of the whole matrix
        _contig_with_marker = list(self.contig_markers.keys())
        _num_marker_contig = len(_contig_with_marker)
        logger.info('There are {} contigs containing marker genes'.format(_num_marker_contig))
        
        _contig_with_marker_index = list(np.sort([self.dict_contigRev[i] for i in _contig_with_marker]))
        _contig_local = self.contig_info[_contig_with_marker_index, :]
        
        _dict_contigRevLocal = {}
        for i, j in enumerate(_contig_local[: , 0]):
            _dict_contigRevLocal[j] = i
        
        for i in range(self.contig_info.shape[0]):
            if i not in _contig_with_marker_index:
                _contig_with_marker_index.append(i)
                
        _normcc_shuffle = self.normcc_matrix.copy().tocsr()[_contig_with_marker_index, :]
        _normcc_shuffle = self.normcc_shuffle.tocsc()[: , _contig_with_marker_index]

        A = _normcc_shuffle.copy()
        A = A - scisp.diags(A.diagonal())
        B = A + scisp.diags((A.sum(axis=0).A.ravel() == 0).astype(int))
        d = scisp.diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B).astype(np.float32)
        Q = random_walk_cpu(P , self.rwr_rp , self.rwr_perc , 0.01 , _num_marker_contig)

        E = Q.copy()
        E += E.T
        d = E.sum(axis=0).A.ravel()
        d[d == 0] = 1
        b = scisp.diags(1 / np.sqrt(d))
        E = b.dot(E).dot(b)
        
        return _contig_local, _dict_contigRevLocal, E.tolil()
        
        
    def _pre_clustering(self):  
        
        _my_gene_counts = list(self.marker_contig_counts.values())
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
            for item in self.marker_contig_counts:
                if self.marker_contig_counts[item] == _g_count:
                    total_contig_lengths = 0
                    
                    #Marker contigs are dict gene: contigs###
                    for contig in self.marker_contigs[item]:
                        contig_mg_counts = len(self.contig_markers[contig])
                        total_contig_lengths += contig_mg_counts

                    total_contig_mgs[item] = total_contig_lengths

            total_contig_mgs_sorted = sorted(
                total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True
            )

            for item in total_contig_mgs_sorted:
                smg_iteration[n] = self.marker_contigs[item[0]]
                n += 1
                
        bins, bin_of_contig, _, _, _ = match_contigs(
            smg_iteration,
            self.contig_markers,
            self.imputed_matrix.tolil(),
            self.dict_contigRevLocal,
            np.percentile(self.imputed_matrix.tocoo().data , self.w_intra),
            np.percentile(self.imputed_matrix.tocoo().data , self.w_inter)
        )
        return bins, bin_of_contig
    
    
    def _unitem_clustering(self):
        #########Use Leiden Algorithm to do clustering########
        _map_del = self.normcc_matrix.tocoo()
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
        part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=_wei, n_iterations = -1)
        part = list(part)
        _dict_cluster = {}
        # dict of communities

        for ci in range(len(part)):
            for id in part[ci]:
                _dict_cluster[self.contig_info[id , 0]] = 'group'+str(ci)
                
        with open(os.path.join(self.path , 'tmp' , 'cluster4unitem.txt'),'w') as out:
            for key , value in _dict_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                
        with open(os.path.join(self.path , 'tmp' , 'dir4unitem.tsv'),'w') as out:
            out.write('INITIAL_BIN' + '\t' + os.path.join(self.path , 'tmp' ,'BIN4unitem'))
                
                
        


