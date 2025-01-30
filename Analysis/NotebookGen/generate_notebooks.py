import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor
import os
from tqdm import tqdm
import pandas as pd
from config import *

def run_notebook(notebook_path, input_string, output_html_path):
    # Load the notebook
    with open(notebook_path) as f:
        notebook = nbformat.read(f, as_version=4)

    # Inject the input string into the first cell (assuming it's a code cell)
    if notebook['cells'][1]['cell_type'] == 'code':
        notebook['cells'][1]['source'] = f"fold_pair = '{input_string}' "
    if notebook['cells'][2]['cell_type'] == 'code':
        notebook['cells'][2]['source'] = f"plot_tool = PlotTool(folder='{DATA_DIR}',fold_pair=fold_pair)"

    # Execute the notebook
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    ep.preprocess(notebook, {'metadata': {'path': os.path.dirname(notebook_path)}})

    # Convert the notebook to HTML
    html_exporter = HTMLExporter()
    html_exporter.exclude_input = True
    body, resources = html_exporter.from_notebook_node(notebook)

    # Save the HTML output
    with open(output_html_path, 'w') as f:
        f.write(body)

    print(f'Notebook converted to HTML and saved to {output_html_path}')

# Example usage
if __name__ == '__main__':
    '''
    fold_pairs  : List of protein pairs you want to run the analysis
    '''
    OUTPUT_PATH = OUTPUT_PATH_NOTEBOOKS
    # fold_pairs = ['4cmqB_4zt0C', '4n9wA_4nc9C', '1nrjB_2gedB']
    fold_pairs = ['4n9wA_4nc9C', '1nrjB_2gedB']

    errors = []
    no_cluster = []
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    notebook_path = os.path.join(project_root, 'TemplateNotebook.ipynb')
    for fold_pair in fold_pairs:
        if '.sh' in fold_pair:
            continue
        try:
            print('##############################################################################')
            print(fold_pair)
            print('##############################################################################')
            output_html_path = f'{OUTPUT_PATH}/{fold_pair}.html'
            run_notebook(notebook_path,fold_pair,output_html_path)
        except Exception as e:
            print(e)
            errors.append({'fold_pair': fold_pair, 'error':e})
            continue
    print('Finish to run !')
    # pd.DataFrame(errors).to_csv('/Users/steveabecassis/Desktop/Pipeline_res/errors_2110.csv',index=False)
