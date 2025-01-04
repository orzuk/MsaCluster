import argparse
import pandas as pd
import os
from polyleven import levenshtein
from hdbscan import HDBSCAN
import numpy as np
import pickle
from glob import glob
from scripts.msa_utils import *




def consensusVoting(seqs):
    ## Find the consensus sequence
    consensus = ""
    residues = "ACDEFGHIKLMNPQRSTVWY-"
    n_chars = len(seqs[0])
    for i in range(n_chars):
        baseArray = [x[i] for x in seqs]
        baseCount = np.array([baseArray.count(a) for a in list(residues)])
        vote = np.argmax(baseCount)
        consensus += residues[vote]
    return consensus


def encode_seqs(seqs, max_len=108, alphabet=None):
    if alphabet is None:
        alphabet = "ACDEFGHIKLMNPQRSTVY-"

    arr = np.zeros([len(seqs), max_len, len(alphabet)])
    for j, seq in enumerate(seqs):
        for i, char in enumerate(seq):
            for k, res in enumerate(alphabet):
                if char == res:
                    arr[j, i, k] += 1
    return arr.reshape([len(seqs), max_len * len(alphabet)])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=
                                """
                                Cluster sequences in a MSA using DBSCAN algorithm and write .a3m file for each cluster.
                                Assumes first sequence in fasta is the query sequence.
                            
                                H Wayment-Steele, 2022
                                """)

    p.add_argument("--keyword", action="store", help="Keyword to call all generated MSAs.")
    p.add_argument("-i", action='store', help='fasta/a3m file of original alignment.')
    p.add_argument("-o", action="store", help='name of output directory to write MSAs to.')
    p.add_argument("--n_controls", action="store", default=10, type=int,help='Number of control msas to generate (Default 10)')
    p.add_argument('--verbose', action='store_true', help='Print cluster info as they are generated.')

    # p.add_argument('--scan', action='store_true', help='Select eps value on 1/4 of data, shuffled.')
    # p.add_argument('--eps_val', action='store', type=float, help="Use single value for eps instead of scanning.")
    p.add_argument('--resample', action='store_true',help='If included, will resample the original MSA with replacement before writing.')
    p.add_argument("--gap_cutoff", action='store', type=float, default=0.25, help='Remove sequences with gaps representing more than this frac of seq.')
    p.add_argument('--min_eps', action='store', default=3, help='Min epsilon value to scan for DBSCAN (Default 3).')
    p.add_argument('--max_eps', action='store', default=20, help='Max epsilon value to scan for DBSCAN (Default 20).')
    p.add_argument('--eps_step', action='store', default=.5, help='step for epsilon scan for DBSCAN (Default 0.5).')
    p.add_argument('--min_samples', action='store', default=3,help='Default min_samples for DBSCAN (Default 3, recommended no lower than that).')
    p.add_argument('--min_cluster_size', action='store', default=15,help='HDBSCAN minimum cluster size.')


    # p.add_argument('--run_PCA', action='store_true', help='Run PCA on one-hot embedding of sequences and store in output_cluster_metadata.tsv')
    # p.add_argument('--run_TSNE', action='store_true', help='Run TSNE on one-hot embedding of sequences and store in output_cluster_metadata.tsv')

    args = p.parse_args()

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
    print("Seqs len: " + str(L))

    df['frac_gaps'] = [x.count('-') / L for x in df['sequence']]

    former_len = len(df)
    df = df.loc[df.frac_gaps < args.gap_cutoff]

    new_len = len(df)
    lprint(args.keyword, f)
    lprint("%d seqs removed for containing more than %d%% gaps, %d remaining" % (
    former_len - new_len, int(args.gap_cutoff * 100), new_len), f)
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)

    n_clusters = []
    eps_test_vals = np.arange(args.min_eps, args.max_eps + args.eps_step, args.eps_step)

    # Save inputs:
    with open('hdbscan_input.pkl', 'wb') as f_pickle:  # Python 3: open(..., 'wb')
        pickle.dump([ohe_seqs, df, IDs, seqs, args], f_pickle)

    clustering = HDBSCAN(min_cluster_size=10, min_samples=args.min_samples, cluster_selection_method='leaf')
    clustering.fit_predict(ohe_seqs)
    clusters = np.unique(clustering.labels_)
    if len(clusters) > 35:
            clustering = HDBSCAN(min_cluster_size=15, min_samples=args.min_samples, cluster_selection_method='eom')
            clustering.fit_predict(ohe_seqs)
            clusters = np.unique(clustering.labels_)

    # _ = clustering.condensed_tree_.plot(select_clusters=True,selection_palette=sns.color_palette("deep", np.unique(clusters).shape[0]),label_clusters=True)
    lprint("%d total seqs" % len(df), f)
    df['dbscan_label'] = clustering.labels_

    clusters = [x for x in df.dbscan_label.unique() if x >= 0]
    unclustered = len(df.loc[df.dbscan_label == -1])

    lprint('%d clusters, %d of %d not clustered (%.2f)' % (len(clusters), unclustered, len(df), unclustered / len(df)),f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in
                                 df.loc[df.dbscan_label == -1]['sequence'].tolist()])
    lprint('avg identity to query of unclustered: %.2f' % avg_dist_to_query, f)

    avg_dist_to_query = np.mean([1 - levenshtein(x, query_['sequence'].iloc[0]) / L for x in df.loc[df.dbscan_label != -1]['sequence'].tolist()])
    lprint('avg identity to query of clustered: %.2f' % avg_dist_to_query, f)

    print("Delete old MSA files: " + args.o + '/' + args.keyword + '*.a3m')
    for f_old in glob(args.o + '/*.a3m'):
        os.remove(f_old)
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

        write_fasta(tmp.SequenceName.tolist(), tmp.sequence.tolist(), outfile=args.o + '/' + args.keyword + '_' + "%03d" % clust + '.a3m')

    print('Saved this output to %s.log' % args.keyword)
    f.close()

