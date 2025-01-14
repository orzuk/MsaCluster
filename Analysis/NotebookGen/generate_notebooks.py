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
    errors = []
    no_cluster = []
    notebook_path = '/Users/steveabecassis/Desktop/Pipeline_res/TemplateNotebook.ipynb'
    files_doned = os.listdir('/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_2110')
    fold_pairs = ['4cmqB_4zt0C', '4n9wA_4nc9C', '1nrjB_2gedB']
    for fold_pair in fold_pairs:
        if f'{fold_pair}.html' in files_doned:
            continue
        if '.sh' in fold_pair:
            continue
        try:
            if(len(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/output_msa_cluster'))<2):
                no_cluster.append(fold_pair)
                continue
            print('##############################################################################')
            print(fold_pair)
            print('##############################################################################')
            input_string     = fold_pair
            output_html_path = f'/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_2110/{input_string}.html'
            run_notebook(notebook_path,input_string,output_html_path)
        except Exception as e:
            print(e)
            errors.append({'fold_pair': fold_pair, 'error':e})
            continue
    pd.DataFrame(errors).to_csv('/Users/steveabecassis/Desktop/Pipeline_res/errors_2110.csv',index=False)
