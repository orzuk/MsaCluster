# MsaCluster

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sortable Table</title>
    <style>
        table {
            border-collapse: collapse;
            width: 100%;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        th {
            cursor: pointer;
            background-color: #f2f2f2;
        }
    </style>
</head>
<body>
    <table id="myTable">
        <thead>
            <tr>
                <th onclick="sortTable(0)">Fold Pair</th>
                <th onclick="sortTable(1)">Best TM score AF 1</th>
                <th onclick="sortTable(2)">Best TM score AF 2</th> <!-- Changed header for clarity -->
            </tr>
        </thead>
        <tbody>
            <tr>
                <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/1iytA_2naoF.html" target="_blank">1iytA_2naoF</a></td>
                <td>0.79</td>
                <td>0.79</td>
            </tr>
            <tr>
                <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/5i2sA_5i2mA.html" target="_blank">5i2sA_5i2mA</a></td>
                <td>0.79</td>
                <td>0.79</td>
            </tr>
            <tr>
                <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/2qqjA_4qdsA.html" target="_blank">2qqjA_4qdsA</a></td>
                <td>0.79</td>
                <td>0.79</td>
            </tr>
            <tr>
                <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/3l5nB_2a73B.html" target="_blank">3l5nB_2a73B</a></td>
                <td>0.79</td>
                <td>0.79</td>
            </tr>
            <!-- More rows as needed -->
        </tbody>
    </table>

    <script>
        function sortTable(n) {
            var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
            table = document.getElementById("myTable");
            switching = true;
            // Set the sorting direction to ascending:
            dir = "asc"; 
            // Make a loop that will continue until no switching has been done:
            while (switching) {
                // Start by saying: no switching is done:
                switching = false;
                rows = table.rows;
                // Loop through all table rows (except the first, which contains table headers):
                for (i = 1; i < (rows.length - 1); i++) {
                    // Start by saying there should be no switching:
                    shouldSwitch = false;
                    // Get the two elements you want to compare, one from current row and one from the next:
                    x = rows[i].getElementsByTagName("TD")[n];
                    y = rows[i + 1].getElementsByTagName("TD")[n];
                    // Check if the two rows should switch place, based on the direction, asc or desc:
                    if (dir == "asc") {
                        if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                            // If so, mark as a switch and break the loop:
                            shouldSwitch = true;
                            break;
                        }
                    } else if (dir == "desc") {
                        if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
                            shouldSwitch = true;
                            break;
                        }
                    }
                }
                if (shouldSwitch) {
                    // If a switch has been marked, make the switch and mark that a switch has been done:
                    rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                    switching = true;
                    // Each time a switch is done, increase this count by 1:
                    switchcount ++;      
                } else {
                    // If no switching has been done AND the direction is "asc",
                    // set the direction to "desc" and run the while loop again.
                    if (switchcount == 0 && dir == "asc") {
                        dir = "desc";
                        switching = true;
                    }
                }
            }
        }
    </script>
</body>
</html>

    
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




