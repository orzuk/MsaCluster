import sys
from pathlib import Path
import pandas as pd

# --- repo imports (top-level only) ---
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import *            # expects MAIN_DIR, etc.


def dataframe_to_latex_longtable(df):
    # Start the LaTeX table
    latex_str = r"\begin{longtable}{| " + " | ".join(
        ["l" if i == 0 else "c" for i in range(len(df.columns))]) + r" |}" + "\n"
    latex_str += r"\hline" + "\n"

    # Add column headers
    headers = " & ".join([r"\textbf{" + col + "}" for col in df.columns]) + r" \\ \hline" + "\n"
    latex_str += headers

    # Iterate over each row in the DataFrame and add to the LaTeX table
    for _, row in df.iterrows():
        row_str = " & ".join(str(value) if value != '-' else '-' for value in row) + r" \\ \hline" + "\n"
        latex_str += row_str

    # Add caption and label
    latex_str += r"\caption[Results Summary]{Results Summary}" + "\n"
    latex_str += r"\label{tab:best_comparison}" + "\n"
    latex_str += r"\end{longtable}"

    return latex_str

if __name__ == '__main__':
    df = pd.read_parquet(TABLES_RES + '/final_res_df_2510.parq')
    df = df[['fold_pair','BEST_AF1','BEST_AF2','BEST_RECALL_FOLD1','BEST_RECALL_FOLD2','BEST_ESM1_AF','BEST_ESM2_AF']]
    df.columns = [ 'fold_pair', 'AF\_F1', 'AF\_F2', 'MSAT\_F1','MSAT\_F2','ESM\_AF1', 'ESM\_AF2']
    df.sort_values(by='fold_pair',inplace=True)
    #df['fold_pair'] = df['fold_pair'].apply(lambda x: x.split('_')[0] + '\_' + x.split('_')[1])
    df['url'] = df['fold_pair'].apply(lambda x:'\href{https://steveabecassis.github.io/MsaCluster/HTML/'+x+'.html}'+'{'+ x.split('_')[0] + '\_' + x.split('_')[1] +'}')
    # Generate LaTeX table string
    df_ = df[[ 'url', 'AF\_F1', 'AF\_F2', 'MSAT\_F1','MSAT\_F2','ESM\_AF1', 'ESM\_AF2']]
    latex_table = dataframe_to_latex_longtable(df_)
    print(latex_table)