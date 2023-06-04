import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import sys, os
from polyleven import levenshtein
from sklearn.cluster import DBSCAN
from utils import *
from hdbscan import HDBSCAN
import seaborn as sns




if __name__ == '__main__':


    class args:
        # i = '/Users/steveabecassis/PycharmProjects/AF_Cluster/data_sep2022/01_KaiB_eps_scan/2QKEE_colabfold.a3m'
        i = {'/Users/steveabecassis/PycharmProjects/AF_Cluster/inputs/AF_cluster_73815.a3m'}
        o = '/Users/steveabecassis/PycharmProjects/AF_Cluster/output'
        keyword = '2QKEE'
        resample = False
        gap_cutoff = 0.25
        min_eps = 3
        max_eps = 20
        eps_step = .5
        verbose = True
        n_controls = 10
        scan = None
        eps_val = None
        min_samples = 3




    os.makedirs(args.o, exist_ok=True)
    f = open("%s.log" % args.keyword, 'w')
    IDs, seqs = load_fasta(args.i)

    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment

    df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs})

    query_ = df.iloc[:1]
    df = df.iloc[1:]

    if args.resample:
        df = df.sample(frac=1)

    L = len(df.sequence.iloc[0])
    N = len(df)

    df['frac_gaps'] = [x.count('-') / L for x in df['sequence']]

    former_len = len(df)
    df = df.loc[df.frac_gaps < args.gap_cutoff]

    new_len = len(df)
    lprint(args.keyword, f)
    lprint("%d seqs removed for containing more than %d%% gaps, %d remaining" % (former_len - new_len, int(args.gap_cutoff * 100), new_len), f)
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)

    n_clusters = []
    eps_test_vals = np.arange(args.min_eps, args.max_eps + args.eps_step, args.eps_step)

    if args.eps_val is None:  # performing scan
        lprint('eps\tn_clusters\tn_not_clustered', f)

    # perform actual clustering

    # clustering = HDBSCAN(min_cluster_size=50, min_samples=args.min_samples).fit_predict(ohe_seqs)
    clustering = HDBSCAN(min_cluster_size = 15, min_samples=args.min_samples)
    clustering.fit_predict(ohe_seqs)
    clusters = np.unique(clustering.labels_)
    _ = clustering.condensed_tree_.plot(select_clusters=True,selection_palette=sns.color_palette("deep", np.unique(clusters).shape[0]), label_clusters=True)
    lprint("%d total seqs" % len(df), f)

    df['dbscan_label'] = clustering.labels_

    clusters = [x for x in df.dbscan_label.unique() if x >= 0]
    unclustered = len(df.loc[df.dbscan_label == -1])

    lprint('%d clusters, %d of %d not clustered (%.2f)' % (len(clusters), unclustered, len(df), unclustered / len(df)),f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in df.loc[df.dbscan_label == -1]['sequence'].tolist()])
    lprint('avg identity to query of unclustered: %.2f' % avg_dist_to_query, f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in
                                 df.loc[df.dbscan_label != -1]['sequence'].tolist()])
    lprint('avg identity to query of clustered: %.2f' % avg_dist_to_query, f)

    cluster_metadata = []
    for clust in clusters:
        tmp = df.loc[df.dbscan_label == clust]

        cs = consensusVoting(tmp.sequence.tolist())

        avg_dist_to_cs = np.mean([1 - levenshtein(x, cs) / L for x in tmp.sequence.tolist()])
        avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in tmp.sequence.tolist()])

        if args.verbose:
            print('Cluster %d consensus seq, %d seqs:' % (clust, len(tmp)))
            print(cs)
            print('#########################################')
            for _, row in tmp.iterrows():
                print(row['SequenceName'], row['sequence'])
            print('#########################################')

        tmp = pd.concat([query_, tmp], axis=0)

        cluster_metadata.append({'cluster_ind': clust, 'consensusSeq': cs, 'avg_lev_dist': '%.3f' % avg_dist_to_cs,
                                 'avg_dist_to_query': '%.3f' % avg_dist_to_query, 'size': len(tmp)})

        write_fasta(tmp.SequenceName.tolist(), tmp.sequence.tolist(),outfile=args.o + '/' + args.keyword + '_' + "%03d" % clust + '.a3m')


    outfile = args.o + "/" + args.keyword + '_clustering_assignments.tsv'
    lprint('wrote clustering data to %s' % outfile, f)
    df.to_csv(outfile, index=False, sep='\t')

    metad_outfile = args.o + "/" + args.keyword + '_cluster_metadata.tsv'
    lprint('wrote cluster metadata to %s' % metad_outfile, f)
    metad_df = pd.DataFrame.from_records(cluster_metadata)
    metad_df.to_csv(metad_outfile, index=False, sep='\t')

    print('Saved this output to %s.log' % args.keyword)
    f.close()

