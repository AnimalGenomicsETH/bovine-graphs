## Instruction to install the graph pipeline 

Tested on Linux 64-bit

- Clone the repository and download the graph software

```

git clone https://github.com/danangcrysnanto/bovine-graphs.git
curl "https://zenodo.org/record/4393273/files/software.tar.gz?download=1" --output graph_software.gz 
tar -xzvf graph_software.gz 

```

- Add the graph sofware to the path (uncomment next line) or add to your $PATH manually

```
# export PATH=$PWD/software_graph:$PATH
```

- Install miniconda (uncomment next line), skip this step if conda already installed

```

# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
# faster resolving conda env using mamba (uncomment)
# conda install mamba -n base -c conda-forge
```

- Create conda environment and activate it
For manual install, required sofware are listed in `sofware_used.tsv`

```

# create conda environment
cd bovine-graphs
mamba env create -f envs/environment_test.yml
conda activate pangenome_test

# install python packages not in conda
pip install markdown-captions

```

- Run test 

```

# run test 
# remove --dry-run to run the real test, will produce a report (grtest_report.pdf) in test/reports folder
snakemake --dry-run -j 1 -s snake_graph.py 

# to produce the pipeline visualization
snakemake --rulegraph -j 1 -s snake_graph.py  | dot -Tpdf > rulegraph.pdf
```

Once successful it will generate a pdf report `(grtest_report.pdf)` in `test/reports` folder
