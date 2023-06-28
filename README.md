# MsaCluster


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


## HURCS cluster

### Virtual environement
Create a virtual environement.
If you are running on one of the HURCS , to [create a vitual environement](https://wiki.cs.huji.ac.il/hurcs/software/python) run:
```
virtualenv /sci/labs/pi-user-name/your-user-name/my-python-venv --python python3
```
Install the requirements
```
pip install -r requirements.txt
```



