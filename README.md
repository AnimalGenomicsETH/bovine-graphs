[![Snakemake](https://img.shields.io/badge/snakemake-≥5.26.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Pangenome Graph Pipeline

Pipeline to integrate multiple genome assemblies into a pangenome graph representation.
This pipeline wraps the functionality [Minigraph](https://github.com/lh3/minigraph) for genome graph construction. 
We add the utility to labels nodes with the assemblies it derives, which is crucial for many pangenome analyses. 
Analyses which are common in pangenome studies are also performed, including:    

- Determination of core and flexible genome    
- Analysis of non-reference sequences    
- Extraction of the structural variations     
- Novel genes predictions and corresponding expression level     
- Variants (SNP and Indels) called from non-reference sequences     

See [Here](pipeline_scheme.pdf) for the scheme of the pipeline and [Here](reports/taurus_report.pdf) for an example of the generated report. 

Developed for analysis of bovine genomes, but should be applicable to the other species as well.      

## Citation
Please cite below when using the pipeline/scripts in your research


> Danang Crysnanto, Alexander S. Leonard, Zih-Hua Fang, Hubert Pausch. **Novel functional sequences uncovered through a bovine multi-assembly graph**. *[Biorxiv](https://www.biorxiv.org/content/10.1101/2021.01.08.425845v1.full)*

## Pipeline usage

**Input**
---

- [Minigraphs](https://github.com/lh3/minigraph) and [Gfatools](https://github.com/lh3/gfatools) need to be installed and available in `$PATH`. Please download [Here](https://doi.org/10.5281/zenodo.4393273) to get the same version used in the paper. 
Required python packages, R libraries, and bioinformatic softwares are listed [Here](envs/software_used.tsv). Alternatively, one could use `mamba / conda`
to create an environment with all softwares installed (Minigraph not included). To generate `pdf` report one need to install [weasyprint](https://weasyprint.org/start/).

```
conda env create -f envs/environment.yml
conda activate pangenome 
```

- Set all parameters required in `config/config.yaml`. All paths will be interpreted relative to the `workdir` directory. 
- Multiple assemblies in the `workdir/assembly` with the naming scheme of {population}.fa. The prefix will be used as the identifier for the assembly.
- A file `config/graph_comp.tsv` that mapping the name of the graph and the order of inclusion.
First in the order used as the graph backbone, which is usually the reference genome (e.g., UCD).       
Construction of multiple graphs can be specified in the different line e.g., 

``` 
graph1 UCD,OBV,Angus 
graph2 UCD,Angus 
```
- Cluster job specification in the pipeline designed for `LSF /bsub` system. We provide `snakemake profile` to run on LSF system, modified from [Here](https://github.com/Snakemake-Profiles/lsf) by Alex. One needs to adapt for the other computing clusters. Please follow guidelines [Here](https://github.com/snakemake-profiles). Pipeline can also be run locally without cluster executions. 


**Usage**    
---

```
# local execution
snakemake -s snake_graph.py

# LSF cluster execution
snakemake --profile "snakemake_profile/lsf" -s snake_graph.py
```

**Output**
---

- Integrated pangenome graphs in `graph` folder with the prefix set in the config e.g., `graph1.gfa`    
- Matrix that map node to the assembly it derives (i.e., node colours/labels), e.g.,    

    | Node | UCD | OBV | Angus |
    | ---- | --- | --- | ----- |
    | s1   | 1   | 0   | 0     |
    | s2   | 0   | 1   | 1     |

    *0 and 1 indicate absence and presence respectively

- Structural variations derived from graphs.      
These are large variations (fragment length > 100 bp) from bubbles in the graph that are
not part of the reference sequences. The SVs are grouped by biallelic and multiallelic. SVs crossing coding sequences visualized using `Graphviz`. Script [app.py](visualize/app.py) can be run, which will set up a local webserver (*not part of the pipeline*) to inspect SVs in a more detailed and in interactive way (*Under development*). 

- Prediction of novel /non-reference genes with corresponding expression levels from the transcriptome. 

- Variants (SNP and Indels) nested in non-reference sequences. One need to run separate pipeline for [variant calling](subworkflows/variant_calling.py). Set config for the pipeline in [Here](config/config_varcall.yaml). 

- Reports in `reports/{graph}` folder. This pdf contains summary of computational resources, statistics of core/flexible genome, structural variations, novel gene models derived from graphs. 
Will output a single pdf from each constructed graph. See the example [Here](reports/taurus_report.pdf).
