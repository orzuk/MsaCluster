#!/bin/bash
#SBATCH --job-name=gen_pair_html
#SBATCH --account=course-52017-25
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=gen_pair_html.out
# Rebuild all per-pair HTML (embeds regenerated 9-method figs + new 3D viewer
# superposition + section.fig <h2> layout).
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
python3 -u TableResults/gen_pair_html.py --pairs ALL --mode inline
