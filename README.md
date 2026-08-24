# TCEanalysis
The code for computational implementation of "Metagenomic mining of translational coupling elements enables programmable gene expression of polycistronic systems in E. coli".

# Project Structure  

TCEanalysis/  
├── README.md                 # Documentation and usage guide  
├── environment.yml           # Conda environment with all dependencies  
├── config.py                 # Configuration parameters (genome_paths, thresholds, etc.)  
├── utils.py                  # Utility functions (I/O, sequence processing)  
├── core_TeRe.py              # Core logic for TeRe‑TCE analysis  
├── core_UTNI.py              # Core logic for UTNI‑TCE analysis  
├── core_Leadinggene.py       # Core logic for leading gene analysis  
├── run_TeRe.py               # Entry point: run TeRe‑TCE analysis  
├── run_UTNI.py               # Entry point: run UTNI‑TCE analysis  
└── run_Leadinggene.py        # Entry point: run leading gene analysis  


#Environment Setup  
For the convenience of researchers, we have packaged all necessary dependencies into a single `environment.yml` file. Users can set up the TCEanalysis environment by running `conda env create -f environment.yml`, which installs all required packages automatically. After activation, the pipeline can be run without additional configuration.  
```
### Step 1: Clone the repository  
```bash  
git clone https://github.com/yourusername/TCEanalysis.git  
cd TCEanalysis  
```

### Step 2: Create and activate conda environment  

```
conda env create -f environment.yml  
conda activate tceanalysis  
```

### Install RNAhybrid (External Dependency)  

RNAhybrid is not available via Conda and must be compiled from source:  

```
# Download source code  
wget https://bibiserv.cebitec.uni-bielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz  

# Extract and compile  
tar -xzvf RNAhybrid-2.1.2.tar.gz  
cd RNAhybrid-2.1.2  
./configure  
sudo make  
sudo make install  

# Verify installation  
RNAhybrid --version  
```

# Dataset Preparation

All analysis scripts for TCEanalysis are provided as part of this work. The genomic datasets analyzed here are publicly available and can be downloaded using the RefSeq accession numbers listed in Supplementary Data 1.
