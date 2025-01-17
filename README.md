# MsaCluster

You can find the full results table in this [link](https://steveabecassis.github.io/MsaCluster/table.html).

## Running from python: 
The file "run_foldswitch_pipeline.py" enables to run the entire pipeline or parts of it from python. This program can send jobs as needed. Usage (from command line): 

python3 run_foldswitch_pipeline [operation] [foldpair_id]

For example: 

python3 run_foldswitch_pipeline cluster_msa 1jfkA_2nxqB

will cluster the MSAs of the family that have sequences folding according to pdbids 1fjk (chain A) and 2nxq (chain B). 

The second argument (foldpair_id) is optional. If not given, the program will loop over all pairs of ids


## Entire Pipeline - Example: <br>
The file 'pipeline.sh' implements an entire pipeline of computing cluster-specific contact maps and structures for a protein family <br>
**Input:** A fasta sequence representing a protein chain (1st argument) <br>
Directory name for outputs (2nd argument) <br>
**Output:** Predicted structures and attention maps for each cluster <br>
**Running example:**
sbatch -out "run_pipeline.out" ./pipeline_with_params.sh  ./input/2qke.fasta 2qke <br>

The output will be in the directory according to the second argument. In this case: <br>
Pipeline/2qke


**Steps:**
1. Genreate a Multiple Sequence Alignment (MSA) for the query sequence
2. Cluster the MSA into different sub-MSAs
3. Run AlphaFold and make prediction for the query sequence with each sub-MSAs
4. Run ESMfold on sampled sequences from each cluster
5. Compute attention maps (predicted contact map) for each cluster based on ESM's MSA-transformer model.
6. Plot predictedattention maps va. true ones

The specific steps can be run individually, as shown in the following commands: 

## Get msa (alignement)
To get the full msa run 

```
python3 ./get_msa.py ./FASTA_FILE.fasta ./OUTPUT_DIR -name BIG_MSA
```
where 

- _FASTA_FILE_ is a fasta file in format .fasta and contain the sequence query
- _BIG_MSA_ is the msa output name (.a3m file)
- _OUTPUT_DIR_ is the output directory 


## Run ClusterMSA

To run ClusterMSA 
```
python3 ./ClusterMSA_moriah.py --keyword MSA -i ./output_dir/BIG_MSA.a3m -o ./results/cluster_msa_output 
```
where 
- _MSA_  is the perfix name of each msa cluster (.a3m file)
- _cluster_msa_output_ is your output foler that will contains all the msa cluster (.a3m file)


## Run AlphaFold 

To run AlphaFold with each msa cluster as feature run 
```
python3 ./RunAF_colabfold.py ./results/cluster_msa_output  af_cmsa_output
```
where 
- _cluster_msa_output_  is the folder where you have all your msas (.a3m files)
- _af_cmsa_output_ is your output foler that will contains all the alphafold predictions (.pdb files)


## Run sample Alphafold  

To get Alphafold predictions on sequences sampled from each msa cluster run
 
```
python3 ./get_sample_msa_algn.py ./results/cluster_msa_output 
```
where 

- _./results/cluster_msa_output _ is your folder where you have all your msas clusters (.a3m files)
- _output_ is your output foler that will contains all the msas (.a3m files)


## Run ESM 

To run ESM for contact map predition 
```
python3  ./runESM.py ./msa_file.a3m ./output
``` 

- _msa_file_ is your msa in .a3m format
- _output_ is your output foler that will contains all the msas (.a3m files)

## Pipeline
To run the full pipeline you need a sequence as input in .fasta format and run 
```
mkdir output_pipeline/output_get_msa 
python3 ./get_msa.py ./fasta_files/rcsb_pdb_1EBO.fasta  ./output_pipeline/output_get_msa  -name '1ebo'
mkdir output_pipeline_1ebo/cluster_msa_output 
python3  ./ClusterMSA_moriah.py --keyword cluster -i ./output_pipeline/output_get_msa/1ebo.a3m -o ./output_pipeline/cluster_msa_output
mkdir output_pipeline/sample_msa 
python3 ./get_sample_msa_algn.py ./output_pipeline/cluster_msa_output  ./output_pipeline/sample_msa
``` 

## HURCS cluster

(Note: All the .sh files are in the folder _/sci/home/steveabecassis/colabfold_new_ )

Activate the virtual environement :
```
source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
```

or (depending on shell) 
```
. /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate
```


### Get msa (alignement)

This is the get_msa.sh file
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


python3 /sci/home/steveabecassis/colabfold_new/get_msa.py PATH_TO_YOUR_SEQUENCE_FILE OUTPUT_PATH  -name 'PREFIX_MSA_NAME'
```

After you replace _PATH_TO_YOUR_SEQUENCE_FILE_ by the path of your sequence file (.a3m or .fasta file),_PREFIX_MSA_NAME_ by the prefix you want for the msa and _OUTPUT_PATH_ by the path you want your msa output send the job with: 
```
sbatch get_msa.sh
```

### ClusterMSA

This is the ClusterMsa.sh file
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
python3  /sci/home/steveabecassis/colabfold_new/ClusterMSA_moriah.py --keyword PREFIX_MSAS -i BIG_MSA_INPUT -o OUTPUT_PATH
```

After you replace _PREFIX_MSAS_ ,_BIG_MSA_INPUT_ (.a3m or .fasta file),_PREFIX_MSAS_ and _OUTPUT_PATH_ send the job with: 
```
sbatch ClusterMsa.sh
```

### Run Alphafold

This is the colabfold_dina_gpu.sh file
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

module load cuda/11.1
module load cudnn/8.0.5

source /sci/labs/dina/dina/collabFold_phoenix/bin/activate.csh

/sci/labs/dina/dina/localcolabfold/colabfold-conda/bin/colabfold_batch  INPUT_MSA_SEQUENCE  OUTPUT_PATH --data=/cs/labs/dina/seanco/colabfold/weights/ 
```

After you replace _INPUT_MSA_SEQUENCE_ ,_OUTPUT_PATH_ send the job with: 
```
sbatch colabfold_dina_gpu.sh
```
Note:
You can give to alphafold as input a sequence (.fasta file) or a msa (.a3m)

### Run sample MSA 

This is the sample_msa.sh file
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


python3 /sci/home/steveabecassis/colabfold_new/get_sample_msa_algn.py MSAS_FOLDER  OUTPUT_PATH
```

After you replace _MSAS_FOLDER_ , _OUTPUT_PATH_ send the job with: 
```
sbatch sample_msa.sh
```
Note:
For each msa in MSAS_FOLDER:
  1. sample sequence
  2. get large msa for this sequence
  3. filter the sequence of the large msa and keep only the sequence that were in the original msa

The original msas are msas output of ClusterMsa.


### PIPELINE

This is the pipeline_msas.sh file:
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

echo $(curl ifconfig.me)
source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh

mkdir output_pipeline_1ebo/output_get_msa 
python3 ./get_msa.py ./fasta_files/rcsb_pdb_1EBO.fasta  ./output_pipeline_1ebo/output_get_msa  -name '1ebo'
mkdir output_pipeline_1ebo/cluster_msa_output 
python3  ./ClusterMSA_moriah.py --keyword cluster -i ./output_pipeline_1ebo/output_get_msa/1ebo.a3m -o ./output_pipeline_1ebo/cluster_msa_output
mkdir output_pipeline_1ebo/sample_msa 
python3 ./get_sample_msa_algn.py ./output_pipeline_1ebo/cluster_msa_output  ./output_pipeline_1ebo/sample_msa
```
The pipeline input is a sequence.
1. Get the msa
2. Cluster the msa into severall msas
3. For each msa cluster sample a sequence and get the alignement


## Other Scripts

The repository includes additional scripts for analysis and table generation, organized as follows:

### Analysis Folder

- **`cmap_analysis.py`**: Generates the analysis of contact maps created by MSA Transformers for each fold pair.
- **`cmap_viz.py`**: Creates the contact maps that appear in the visualization notebook.
- **`esmfold_analysis.py`**: Analyzes the outputs generated by ESMFold.
- **Notebook Generator**: A script to generate Jupyter notebooks for result visualization and exploration.

### TableResults Folder

- **`gen_html_table.py`**: Generates the HTML table displayed in the GitHub repository.
- **`gen_latex_table.py`**: Generates the LaTeX table used in Overleaf.
- **`get_raw_data.py`**: Retrieves raw data at the cluster level without aggregations.
- **`get_tm_align_score.py`**: Calculates the TM score between each fold pair (e.g., to add to the results in the repository).
- **`summary_table.py`**: Generates a summary table used by other scripts to create tables in formats like LaTeX and HTML.

### Pipeline Folder

Contains `.sh` files for running scripts on the HURCS cluster. These include pipeline and task-specific scripts to facilitate parallelized or automated execution.

---

With this structure, you can easily locate the script you need based on its functionality and purpose. For more detailed usage instructions, refer to the respective script documentation or comments within the code.




