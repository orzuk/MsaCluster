import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor
import os
from tqdm import tqdm
import pandas as pd

def run_notebook(notebook_path, input_string, output_html_path):
    # Load the notebook
    with open(notebook_path) as f:
        notebook = nbformat.read(f, as_version=4)

    # Inject the input string into the first cell (assuming it's a code cell)
    if notebook['cells'][1]['cell_type'] == 'code':
        notebook['cells'][1]['source'] = f"fold_pair = '{input_string}' "

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
    OUTPUT_PATH : The output you want the html notebook outputs
    fold_pairs  : List of protein pairs you want to run the analysis
    '''
    OUTPUT_PATH = '/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_2110'
    fold_pairs = ['4cmqB_4zt0C', '4n9wA_4nc9C', '1nrjB_2gedB']
    errors = []
    no_cluster = []
    notebook_path = './MsaCluster/Analysis/NotebookGen/TemplateNotebook.ipynb'
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
