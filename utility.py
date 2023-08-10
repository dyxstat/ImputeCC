#!/usr/bin/env python
# coding: utf-8
import warnings
##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")
import logging
import sys
import os
import numpy as np
import pandas as pd
import operator
import time
import scipy.sparse as scisp
from scipy.sparse import csr_matrix, diags, eye
from scipy.sparse.linalg import norm
import igraph as ig
import leidenalg
import networkx as nx
import copy
from sklearn import metrics
import bz2
import pickle
import gzip
import io
import subprocess

logger = logging.getLogger(__name__)


'''
Function Section
'''
def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with open_output(file_name, compress='gzip') as out_h:
        pickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialzied object
    """
    with open_input(file_name) as in_h:
        return pickle.load(in_h)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognising gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.BZ2File(file_name, 'r')
    elif suffix == 'gz':
        return gzip.GzipFile(file_name, 'r')
    else:
        return open(file_name, 'r')


def open_output(file_name, append=False, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param file_name: file name of output
    :param append: append to any existing file
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """

    mode = 'w' if not append else 'w+'

    if compress == 'bzip2':
        if not file_name.endswith('.bz2'):
            file_name += '.bz2'
        # bz2 missing method to be wrapped by BufferedWriter. Just directly
        # supply a buffer size
        return bz2.BZ2File(file_name, mode, buffering=65536)
    elif compress == 'gzip':
        if not file_name.endswith('.gz'):
            file_name += '.gz'
        return io.BufferedWriter(gzip.GzipFile(file_name, mode, compresslevel=gzlevel))
    else:
        return io.BufferedWriter(io.FileIO(file_name, mode))


def calc_sparsity(matrix):
    row, col = matrix.shape
    sparsity = matrix.nnz / row / col
    return sparsity
    

def random_walk_cpu(P, rp, perc, tol, num_marker_contig):
    #logger.info('the original sparsity of normcc-normalized matrix is {}'.format(calc_sparsity(P)))
    _record = 0
    _start_time = time.time()
    n_genes = P.shape[0]
    row = np.arange(num_marker_contig)
    col = np.arange(num_marker_contig)
    data = np.ones(num_marker_contig)
    I = scisp.coo_matrix((data, (row, col)), shape=(num_marker_contig, P.shape[0]), dtype=np.float32)
    Q = I.copy()
    delta_previous = 0
    for i in range(500):
        Q_new = (1 - rp) * Q.dot(P) + rp * I
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        
        if i >= 1:
            min_cutoff = np.percentile(Q.tocoo().data, perc)
            if min_cutoff > 0:
                #s_before = calc_sparsity(Q)
                Q = Q.multiply(Q > min_cutoff)
                #s_after = calc_sparsity(Q)
            
        if delta < tol:
            _end_time = time.time()
            logger.info('End because of tolerance condition')
            break

        if np.abs(delta - delta_previous) < 0.001:
            _record += 1
            
        delta_previous = delta
        if _record == 5:
            logger.info('End because of early stop')
            break
    _end_time = time.time()
    Q = Q.tocsr()[np.arange(num_marker_contig) , :]
    Q = Q.tocsc()[: , np.arange(num_marker_contig)]
    logger.info('Inputaion takes {} time and {} steps; loss is {}; sparsity of imputed matrix is {}'.format(_end_time-_start_time, i+1, delta, calc_sparsity(Q)))
    return Q
    
    


# Get contigs containing marker genes
def get_contigs_with_marker_genes(
    contigs_file, mg_length_threshold, dict_contigRev
):
    '''
    contig_file: hmm file end with '.hmmout'
    mg_length_thres: thres of marker gene
    
    '''
    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(contigs_file + ".hmmout", "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = contig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # Contig name
                contig_name = "_".join(name_strings)
                
                if contig_name in dict_contigRev:

                    if mapped_marker_length > marker_gene_length * mg_length_threshold:
                        marker_repeated_in_contig = False

                        # Get marker genes in each contig
                        if contig_name not in contig_markers:
                            contig_markers[contig_name] = [marker_gene]
                        else:
                            if marker_gene not in contig_markers[contig_name]:
                                contig_markers[contig_name].append(marker_gene)

                        # Get contigs containing each marker gene
                        if marker_gene not in marker_contigs:
                            marker_contigs[marker_gene] = [contig_name]
                        else:
                            if contig_name not in marker_contigs[marker_gene]:
                                marker_contigs[marker_gene].append(contig_name)
                            else:
                                marker_repeated_in_contig = True

                        # Get contig counts for each marker
                        if marker_gene not in marker_contig_counts:
                            marker_contig_counts[marker_gene] = 1
                        else:
                            if not marker_repeated_in_contig:
                                marker_contig_counts[marker_gene] += 1

    return marker_contigs, marker_contig_counts, contig_markers




def count_contigs_with_marker_genes(marker_contig_counts):
    marker_frequencies = {}

    for marker in marker_contig_counts:
        if marker_contig_counts[marker] not in marker_frequencies:
            marker_frequencies[marker_contig_counts[marker]] = 1
        else:
            marker_frequencies[marker_contig_counts[marker]] += 1

    return marker_frequencies
    
    
#smg_iteration: a list of integers representing contig IDs that contain seed marker genes
#bins: a list of lists, where each inner list contains contig IDs that are currently assigned to the same bin
#n_bins: an integer representing the number of bins currently in use
#bin_of_contig: a dictionary mapping each contig ID to its current bin
#binned_contigs_with_markers: a list of contig IDs that have already been assigned to a bin and contain seed marker genes
#bin_markers: a dictionary mapping each bin to a list of seed marker genes present in its member contigs
#contig_markers: a dictionary mapping each contig to a list of seed marker genes present in it
#contig_lengths: a dictionary mapping each contig ID to its length
#hic_matrix: hic_matrix sparse lil format
#w_intra: a float representing the maximum allowed weight for intra-bin edges
#w_inter: a float representing the minimum required weight for inter-bin edges
#d_limit: an integer representing the maximum allowed shortest path length between contigs in the same bin
    
    
def match_contigs(
    smg_iteration,
    contig_markers,
    normcc_matrix,
    dict_contigRev,
    w_intra,
    w_inter
):
    
    edge_weights_per_iteration = {}
    #smg_iterations = len(smg_iteration) ####number of iterations

    for i in range(len(smg_iteration)):
        if i == 0:
        # Initialise bins with the contigs having the first single-copy
        # marker gene according to the ordering
            #logger.info("Initialize" + ": " + "construct " + str(len(smg_iteration[i])) + " initial bins")
            bins = {} #dict[bin index] = [contig name]
            bin_of_contig = {} #dict[contig name]: bin index
            bin_markers = {} #dict[bin index] = [marker genes existing in bins already]

            binned_contigs_with_markers = []

            for i in range(len(smg_iteration[0])):
                binned_contigs_with_markers.append(smg_iteration[0][i])
                contig_num = smg_iteration[0][i]

                bins[i] = [contig_num]
                bin_of_contig[contig_num] = i

                bin_markers[i] = contig_markers[contig_num]
    
        else:
            B = nx.Graph()
            #print("Iteration " + str(i) + ": " + str(len(smg_iteration[i])) + " contig(s) with seed marker genes")
            common = set(binned_contigs_with_markers).intersection(set(smg_iteration[i])) #common show contigs that are not binned in the iteration
            to_bin = list(set(smg_iteration[i]) - common) #remove the already binned contigs from smg_iteration[i]
            #print(str(len(to_bin)) + " contig(s) to bin in the iteration")
            n_bins = len(bins)
            bottom_nodes = []

            for n in range(n_bins):
                contigid = bins[n][0]
                if contigid not in bottom_nodes:
                    bottom_nodes.append(contigid) #找到每个bin的第一个contig存储在buttom node list中用于代表这个bin

            top_nodes = []
            edges = []

            binned_count = 0

            if len(to_bin) != 0: #to_bin shows the unbinner contigs in the ith iteration
                for contig in to_bin:
                    contigid = contig

                    if contigid not in top_nodes:
                        top_nodes.append(contigid) ###top nodes shows the nodes belong to the top level of biparate graph

                    for b in range(n_bins):
                        hic_sum = 0
                        n_contigs = len(bins[b]) #bins[b]: the bth bin; and there is n_contigs contigs in the bth contig
                        
                        for j in range(n_contigs):
                            hic_sum += normcc_matrix[dict_contigRev[contigid] , dict_contigRev[bins[b][j]]]
                            
                        edges.append((bins[b][0], contigid, -hic_sum / n_contigs))

                B.add_nodes_from(top_nodes, bipartite=0)
                B.add_nodes_from(bottom_nodes, bipartite=1)

                edge_weights = {}

                # Add edges only between nodes of opposite node sets
                for edge in edges:
                    edge_weights[(edge[0], edge[1])] = edge[2]
                    B.add_edge(edge[0], edge[1], weight=edge[2])

                edge_weights_per_iteration[i] = edge_weights

                top_nodes = {n for n, d in B.nodes(data=True) if d["bipartite"] == 0}
                bottom_nodes = set(B) - top_nodes

                if len(top_nodes) > 0:
                    my_matching = (
                        nx.algorithms.bipartite.matching.minimum_weight_full_matching(
                            B, top_nodes, "weight"
                        )
                    )   #my matching is double the direction but repetitive (i.e., a to b and b to a)

                    not_binned = {}

                    for l in my_matching:
                        if l in bin_of_contig:
                            b = bin_of_contig[l]

                            if (
                                my_matching[l] not in bins[b]
                                and (l, my_matching[l]) in edge_weights
                            ):

                                if ( -edge_weights[(l, my_matching[l])] >= w_intra):
                                    can_assign = False

                                    common_mgs = set(bin_markers[b]).intersection(
                                        set(contig_markers[my_matching[l]])
                                    )

                                    if len(common_mgs) == 0:
                                        can_assign = True

                                    if can_assign:
                                        bins[b].append(my_matching[l])
                                        bin_of_contig[my_matching[l]] = b
                                        binned_contigs_with_markers.append(
                                            my_matching[l]
                                        )
                                        binned_count += 1

                                        # # logger.debug("To assign contig " + my_matching[l] + " to bin "+str(
                                        #     b) + " based on contig " + str(l) + " weight="+str(edge_weights[(l, my_matching[l])]))

                                        bin_markers[b] = list(
                                            set(
                                                bin_markers[b]
                                                + contig_markers[my_matching[l]]
                                            )
                                        )

                                    else:
                                        not_binned[my_matching[l]] = (l, b)

                                else:
                                    not_binned[my_matching[l]] = (l, b)   #not binner contigs include contigs that are excluded though in the binary mapping
                    
                    #####The following codes used to check whether create a new bin or not#######
                    for nb in not_binned:
                        new_bin = True
                        for b in bottom_nodes:
                            if -edge_weights_per_iteration[i][(b , nb)] > w_inter:
                                new_bin = False
                                
                        if new_bin:
                            #print("Creating new bin..." + str(nb) + " to bin " + str(n_bins+1))
                    
                            bins[n_bins] = [nb]
                            bin_of_contig[nb] = n_bins
                            binned_count += 1
                    
                            bin_markers[n_bins] = contig_markers[nb]
                            n_bins += 1
                            binned_contigs_with_markers.append(nb)

            #print(str(binned_count) + " contig(s) binned in the iteration")

    if len(smg_iteration) > 0:
        del edge_weights_per_iteration
        del B
        del my_matching
        del not_binned
        del edge_weights
        del common
        del to_bin
        del top_nodes
        del bottom_nodes
        del edges

    return bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers

