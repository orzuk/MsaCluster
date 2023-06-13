# MsaCluster

## Virtual environement
Create a virtual environement.
If you are running on one of the HURCS , to [create a vitual environement](https://wiki.cs.huji.ac.il/hurcs/software/python) run:
```
virtualenv /sci/labs/pi-user-name/your-user-name/my-python-venv --python python3
```
Install the requirements
```
pip install -r requirements.txt
```

## Run ClusterMSA
To run ClusterMSA 
1. Modify the ClusterMsa.sh file to run the ClusterMSA on your own Msa : 
```
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=2
#SBATCH --mem=1G



python3  ./ClusterMSA_moriah.py --keyword PREFIXE_MSA_NAME -i ./msa_file.a3m -o ./output 

```
where 

- _PREFIXE_MSA_NAME_  is the prefix name you want before each msa cluster .a3m file
- _msa_file_ is your msa in .a3m format
- _output_ is your output foler that will contains all the msas (.a3m files)

2. Send the job 
```
sbatch ClusterMsa.sh
```

## Run Alphafold

To run Alphafold 
1. Modify the Alphafold.sh file to run the ClusterMSA on your own Msa : 
```
##!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=4G


python3  ./RunAF2_moriah.py --input_msa_folder ./output
```
where 

- _input_msa_folder_  is the folder where you have all your msas .a3m files
- _output_ is your output foler that will contains all the alphafold predictions .pdb files

2. Send the job 
```
sbatch Alphafold.sh
```

## Run ESM 

To run ESM for contact map predition 
1. Modify the ESM.sh file to run the ClusterMSA on your own Msa : 
```
##!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=4G


python3  ./ESM_moriah.py ./msa_file.a3m ./output
```
where 

- _msa_file_ is your msa in .a3m format
- _output_ is your output foler that will contains all the msas (.a3m files)

2. Send the job 
```
sbatch ESM.sh
```




