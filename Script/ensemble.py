#!/usr/bin/env python

#modified from https://github.com/dparks1134/UniteM/blob/master/unitem/ensemble.py
import random
from collections import defaultdict


def calculateN50(seq_lens):
    '''
    bin_len: a list of contig length belong to one bin
    '''    
    thresholdN50 = sum(seq_lens) / 2
    seq_lens.sort(reverse=True)
    
    test_sum = 0
    N50 = 0
    for seq_len in seq_lens:
        test_sum += seq_len
        if test_sum >= thresholdN50:
            N50 = seq_len
            break

    return N50


class Ensemble():
    """
    Bin selection:
    consensus -> Consensus clustering across multiple binning methods #
    greedy    -> Greedy bin selection across multiple binning methods （DAS_tool）
    unanimous -> Unanimous bin filtering across multiple binning methods (Binning_refiner )
    """

    def __init__(self):
        """

        """

    def _matched_bin_sets(self, bins, contig_lens, bin_quality, min_quality, min_comp, max_cont, no_bin_matching, min_perc_common,
                          overlap_strategy='unitem_overlap_strategy', sort_quality_strategy='sort_matches_sum_q',
                          quality_filter=True):
        """Determine all sets of matched bins.

        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contig_lens : d[cid] -> length of contig
          Length of contigs.
        bin_quality : list of tuples
          Bin quality information. Must be sorted by quality!
        min_quality : float
          Minimum quality of bin to consider during bin matching.
        no_bin_matching : boolean
          Flag indicating bin matching should be skipped and all bins treated independently.

        Return
        ------
          Matched bin sets.
        """

        matched_sets = []
        processed_bins = defaultdict(set)
        for cur_bm, cur_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
            if cur_bid in processed_bins[cur_bm]:
                continue  # bin has already been considered

            if not bins[cur_bm][cur_bid]:
                continue  # null bin
            matched_bins = []
            matched_bins.append((cur_bm, cur_bid, quality, N50, gs))
            if quality_filter:
                quality_flag = (quality >= min_quality and comp >= min_comp and cont <= max_cont)
            else:
                quality_flag = True
            if not no_bin_matching and quality_flag:
                binned_bm = []
                binned_bm.append(cur_bm)
                for test_bm, test_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
                    if test_bm in binned_bm:
                        # can't group with a bin from the same binning method
                        continue

                    if not bins[test_bm][test_bid]:
                        continue  # null bin

                    if test_bid in processed_bins[test_bm]:
                        continue  # bin has already been considered

                    if quality_filter:
                        if quality < min_quality or comp < min_comp or cont > max_cont:
                            continue

                    bp_in_common, total_bp1, total_bp2 = self._bases_in_common(bins[cur_bm][cur_bid],
                                                                               bins[test_bm][test_bid],
                                                                               contig_lens)

                    if overlap_strategy == 'unitem_overlap_strategy':
                        per_bp_in_common = bp_in_common / max(total_bp1, total_bp2)
                        if per_bp_in_common > min_perc_common:
                            matched_bins.append((test_bm, test_bid, quality, N50, gs))
                            binned_bm.append(test_bm)
                    elif overlap_strategy == 'metawrap_overlap_strategy':
                        per_bp_in_common = bp_in_common / min(total_bp1, total_bp2)
                        if per_bp_in_common > min_perc_common:
                            matched_bins.append((test_bm, test_bid, quality, N50, gs))
                            binned_bm.append(test_bm)

            # removed matched bins
            for bm, bid, _q, _n50, _gs in matched_bins:
                processed_bins[bm].add(bid)

            matched_sets.append(tuple(matched_bins))

        # sort by total quality of bins in matched sets
        if sort_quality_strategy == 'sort_matches_max_q':
            matched_sets.sort(key=lambda ms: (max([x[2] for x in ms]),
                                              max([x[3] for x in ms]),
                                              max([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        elif sort_quality_strategy == 'sort_matches_sum_q':
            matched_sets.sort(key=lambda ms: (sum([x[2] for x in ms]),
                                              sum([x[3] for x in ms]),
                                              sum([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        elif sort_quality_strategy == 'sort_matches_maxaddsum_q':
            matched_sets.sort(key=lambda ms: (max([x[2] for x in ms]),
                                              sum([x[2] for x in ms]),
                                              max([x[3] for x in ms]),
                                              max([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        return matched_sets

    def _bases_in_common(self, contig_ids1, contig_ids2, contig_lens):
        """Calculate number of base pairs in common between two bins.

        Parameters
        ----------
        contig_ids1 : iterable
          Contigs in bin.
        contig_ids2 : iterable
          Contigs in bin.
        contig_lens : d[seq ID] -> length
          Length of contigs.

        Returns
        -------
          Base pairs in common, total bases in bin 1, total bases in bin 2
        """

        bp_in_common = sum([contig_lens[seq_id]
                            for seq_id in contig_ids2
                            if seq_id in contig_ids1])

        total_bp1 = sum([contig_lens[seq_id] for seq_id in contig_ids1])
        total_bp2 = sum([contig_lens[seq_id] for seq_id in contig_ids2])

        return bp_in_common, total_bp1, total_bp2

    # from unitem Ensemble class
    # update for metabinner
    def _bin_quality(self, markers, bins, contig_lens, methods_sorted, quality_weight):
        """Determine estimated completeness, contamination, and quality of bins.

        Parameters
        ----------
        bins : d[binning method][bin ID] -> set(cid1, cid2, ... cidN)
          Contigs for bins across all binning methods.
        contig_len : d[cid] -> contig length
          Contigs across all bins.
        quality_weight : float
          Weight given to contamination when assessing genome quality.

        Return
        ------
          List with bin metadata sorted by quality, then N50, the genome size.
        """

        # markers = Markers()

        q = []

        for method_id in methods_sorted:
            for bin_id in bins[method_id]:
                domain, comp, cont = markers.bin_quality(bins[method_id][bin_id])
                bin_seqs = [contig_lens[cid] for cid in bins[method_id][bin_id]]
                n50 = calculateN50(bin_seqs)
                genome_size = sum(bin_seqs)
                q.append((method_id, bin_id, domain, comp, cont, comp - quality_weight * cont, n50, genome_size))

        # sort bins by quality follwed by N50 followed by genome size, and
        # break remaining ties randomly
        q.sort(key=lambda x: (x[5], x[6], x[7], random.random()),
               reverse=True)

        return q

    def _reconstruct_match_sets(self, matched_set,
                                bins,
                                contig_lens,
                                em_mode='greedy'):
        """
        changed from unitem ensemble.py def _resolve_matched_set()
        Select contigs to cluster from matched bin set.

        Parameters
        ----------
        matched_set : iterable with bin metadata (binning method, bin ID, quality)
        Matched set of bins. Highest quality bin must be first.
        bins : d[binning method][bin ID]
        Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
        Contigs across all bins.
        bin_quality : list of tuples
        Bin quality information. Must be sorted by quality!!!!!
        remove_perc : float
        Minimum percentage of bins from other binning methods require to remove contig in highest quality bin.
        add_perc : float
        Minimum percentage of matched bins required to add contig to highest quality bin.
        add_matches : float
        Minimum number of matched bins required to 'add' contigs.
        sel_min_quality : float
        Minimum quality of bins to consider for filtering contigs.
        """
        # identify contig count for bins in matched in set

        primary_bm, primary_bid, primary_q, primary_n50, primary_gs = matched_set[0]
        primary_contigs = bins[primary_bm][primary_bid]

        if primary_q == 100:
            new_bin = {}
            removed_contigs = set()
            added_contigs = set()
            for cid in primary_contigs:
                new_bin[cid] = contig_lens[cid]

        else:
            matched_contigs = defaultdict(int)
            matched_bin_ids = set()
            matched_methods = set()
            for bm, bid, _q, _n50, _gs in matched_set:
                matched_bin_ids.add(bm + bid)
                matched_methods.add(bm)
                for cid in bins[bm][bid]:
                    matched_contigs[cid] += 1

            new_bin = {}
            removed_contigs = set()
            added_contigs = set()

            if em_mode == 'greedy':
                for cid in primary_contigs:
                    new_bin[cid] = contig_lens[cid]

        return new_bin, removed_contigs, added_contigs

    def _update_bins(self, bins, cids_to_remove):
        """Remove specified contigs from all bins.

            Parameters
            ----------
            bins : d[binning method][bin ID]
              Contigs for bins across all binning methods.
            cids_to_remove : iterable
              Contigs to remove from bins.
        """

        cids_to_remove = set(cids_to_remove)

        # remove scaffolds in highest quality bin from all other bins
        for method_id in bins:
            for bin_id in bins[method_id]:
                bins[method_id][bin_id] = list(set(bins[method_id][bin_id]) - cids_to_remove)