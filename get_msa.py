from pathlib import Path
import os
from argparse import ArgumentParser

from colabfold.batch import get_msa_and_templates, msa_to_str, get_queries

print("Get MSA, Current dir: ", os.getcwd())

def _ungap(s: str) -> str:
    return (s or "").replace("-", "").replace(".", "").strip()

def main():
    parser = ArgumentParser()
    parser.add_argument("input", help="Seed A3M/FASTA with the two chains aligned")
    parser.add_argument("results", help="Directory to write the results to")
    parser.add_argument("--name", default="DeepMsa", help="MSA base name (default: DeepMsa)")
    args = parser.parse_args()

    result_dir = Path(args.results)
    result_dir.mkdir(parents=True, exist_ok=True)

    # Parse queries (works for A3M/FASTA); use first sequence as query
    queries, is_complex = get_queries(args.input, 'random')
    print("queries:"); print(queries)
    sequence = _ungap(queries[0][1])
    print("sequence:"); print(sequence)
    print("results dir:"); print(result_dir)
    print("args.name="); print(args.name)

    # Read raw A3M lines; some ColabFold builds require passing them explicitly
    with open(args.input, "r") as fh:
        a3m_lines = fh.read().splitlines()

    # Try signature that accepts a3m_lines; fall back if not supported
    try:
        tup = get_msa_and_templates(
            query_sequences=sequence,
            jobname=args.name,
            result_dir=result_dir,
            use_templates=False,
            custom_template_path='',
            pair_mode='',
            msa_mode='mmseqs2_uniref_env',
            a3m_lines=a3m_lines
        )
    except TypeError:
        tup = get_msa_and_templates(
            query_sequences=sequence,
            jobname=args.name,
            result_dir=result_dir,
            use_templates=False,
            custom_template_path='',
            pair_mode='',
            msa_mode='mmseqs2_uniref_env'
        )

    (unpaired_msa, paired_msa, query_seqs_unique,
     query_seqs_cardinality, template_features) = tup

    print("unpaired msa:"); print(unpaired_msa)
    print("paired msa:"); print(paired_msa)
    print("query seqs unique:"); print(query_seqs_unique)

    msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
    result_dir.joinpath(f"{args.name}.a3m").write_text(msa)

if __name__ == "__main__":
    main()
