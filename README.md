### PPI Network Clustering Analysis

#### Overview
This GitHub repository contains all the materials used in the experimental project for the Network Sciences course at Utrecht University. It includes data, scripts, and a Jupyter notebook for analyzing clustering algorithms applied to Protein-Protein Interaction (PPI) networks.

#### Repository Structure
- **data/**: Contains network data files and clustering results in JSON format.
  - `511145.protein.physical.links.v12.0.txt.gz`: Compressed TSV file from STRING database containing the PPI network of E. Coli.
  - `511145.protein.aliases.v12.0.txt.gz`:  Compressed TSV file from STRING database to map nodes to various protein identifiers.
  - `83333.tsv`: TSV file from Complex Portal with data about all the experimentally determined protein complexes in E. Coli.
- **data/results/**: Stores JSON files with clustering results by algorithm.
  - `Girvan_Newman/`
  - `Louvain/`
  - `OH_PIN/`
- **scripts/**: Python scripts and a Jupyter notebook.
  - `OH_PIN_Final.py`: Runs the OH-PIN algorithm.
  - `GirvanNewmanParallel.py`: Runs Girvan-Newman with parallel processing.
  - `Results_Notebook.ipynb`: Reanalyzes clustering results and replicates study figures.


##### Running Scripts
You can execute the clustering algorithms we implemented in our study:
```bash
python scripts/OH_PIN_Final.py data/511145.protein.physical.links.v12.0.txt.gz
python scripts/GirvanNewmanParallel.py data/511145.protein.physical.links.v12.0.txt.gz <num_cores>
```

##### Running Jupiter Notebook
You can clone this repository and run the Results_Notebook.ipynb to reproduce the analysis we did for our project and repordice the report figures and tables, even without running the clustering algorithms yourself (we provide our the clustering outcomes). To run the scripts and the Jupyter notebook, you need the following Python packages:
- `networkx`
- `pandas`
- `matplotlib`
- `seaborn`
- `numpy`
