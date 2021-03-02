## Detailed instructions of the pangenome graph analysis

Step by step (manual) instruction to construct a graph from collection of genome assemblies and characterize structural variations from the pangenome graph. All these steps can be automatically invoked with the pipeline. 

### Pangenome graph construction

1.  Estimate pairwise genetic distance across assemblies 

    The genetic distance across is used to determine the order of inclusion to graph. 
    
    ```sh
    
    # sketch assembly, done separately for each assembly
    mash sketch -p {threads} -o {output} {input_assembly}
    
    # combined all sketches
    mash paste {output} {input_sketch1} {input_sketch2}
    
    # estimate distance based on combined sketches
    mash dist {input_combined_sketch} {input_combined_sketch} > {output_distance}
    
    # visualize the genetic relationship as tree (optional)
    scripts/phylo_tree_assembly.R
    ```


2.  Graph construction

    Choose one assembly as the backbone of the graph. 
    Include assemblies in the graph according to the genetic distance to the assembly backbone, 
    closer assembly added earlier than genetically more far assembly.
    
    ```sh
    
    # graph construction
    minigraph --inv no -xggs -t {threads} {input_assemb1} {input_assemb2}  > {graph.gfa}
    
    # subset graph file for easier access of nodes and edges information
    # without the need for loading the sequences in nodes
    
    # extract node information 
    awk '$1~/S/ { splt($5,chr,":"); split($6,pos,":"); split($7,arr,":");
                print $2,length($3),chr[3],pos[3],arr[3] }' {graph.gfa} > {graph_len.tsv}
    
    # extract edge information
    awk '$1 == "L"' {graph.gfa} > {graph_link.tsv}
    ```

3. Remap assemblies back to the graph
   Separately map the assembly back to the graphs and record the coverage for all nodes and edges in the graph. 
    
    ```sh
    
    minigraph -t {threads} --cov -x asm {graph.gfa} {assembly1.fa} > {graph_use_assembly1.gfa}
    
    ```

4. Combine mapping coverage information across assemblies 

    ```sh
    # custom python script
    
    scripts/comb_coverage.py -g {assemb1} {assemb2} -a {graph_name}
    
    #will output node_coverage.tsv and edge_coverage_use.tsv
    ```


5. Using coverage information to label node in the graph


    ```sh
    # custom R script
    
    scripts/colour_node.R {assemb1} {assemb2} {graph_name} 
     
    ```
    
    Output:    
    
    
    - nodecol.tsv containing node and label information 
    
    | Node |      labels      |
    |------|------------------|
    | s1   | assemb1          |
    | s2   | assemb1, assemb2 |
    
    
    - nodemat.tsv containing presence (1) and absence (0) of the label in each node
    
    | Node | assemb1 | assemb2 |
    |------|---------|---------|
    | s1   |       1 |       0 |
    | s2   |       1 |       1 |

6. Analyze the pangenome based on the label


    ```sh
    # custom R script
    
    scripts/run_core_nonref.R {graph_name}
    
    ```
    
    Four `R functions` will be invoked, which are detailed as below:
    
    |        Function        |                                Description                                |      Output file      |
    |------------------------|---------------------------------------------------------------------------|-----------------------|
    | calculate core         | calculate the proportion of core and flexible part of pangenome           | core_analysis.tsv     |
    | pangenome sampling     | evaluate the graph properties as more assemblies being added to the graph | core_flex_sim.tsv     |
    | calculate non_ref      | calculate non-ref node length from each assemblies                        | nonref_analysis.tsv   |
    | nonref_sharing_pattern | make plot of the non-ref node sharing across assemblies                   | nonref_shared_len.pdf |


### Structural variations analysis

1. Identify bubbles in the graph

    Bubbles is the regions that diverged between assemblies which has common start and final node derived from reference sequences. 
    
    ```sh
    gfatools bubble {graph.gfa} > {bubble.tsv}
    ```


2. Identify location of the structural variations from graphs

    Paths the bubble represent alleles of the structural variations. 
    This step will enumerate all possible paths based on the start and stop node of the bubbles, done separately for biallelic (2 paths/alleles) and multi-allelic  (>3 alleles) bubbles. At the end, trace the path in the origin of the assemblies using node and edge label information. 
    
    ```sh
    # custom Python scripts
    
    # biallelic SV
    scripts/get_bialsv.py -a {assemb1} {assemb2} > {output} #output: biallelic_sv.tsv
    
    # multiallelic SV
    scripts/get_multisv.py -a {assemb1} {assemb2} > {output} #multiallelic_sv.tsv
    
    #trace path in each SV according to the origin of the assemblies
    scripts/trace_path.py -g {assemb1} {assemb2} -a {graph_name} > {output}
    
    ```


3. Annotate the breakpoint location of the graph

    Annotate the breakpoint of the structural variations using start and stop node coordinate for left and right breakpoints, respectively on the    reference backbone coordinate. Require `gff` file of the backbone assembly.
    
    ```sh
    # custom Python script
    
    scripts/annot_breakpoints.py #output bubble_annot.tsv
    ```


4. Extract structural variation alleles in the bubbles

    Extract non-ref alleles (excluding path less than 100 bp, complete deletions, or paths without non-ref sequences) as a representative non-reference sequences of the pangenome. Sequences in nodes are not used directly, because multiple consecutive nodes might be part of the same allele, and they are representing a continous sequences.
    
    ```sh
    #custom Python scripts
    
    # biallelic allele extraction
    scripts/get_bialseq.py -a {assemb1} {assemb2} #output bialsv_seq.fa 
    # multiallelic allele extraction
    scripts/get_multiseq.py -a {assemb1} {assemb2} # output multisv_seq.fa
    # combined biallelic and multiallelic SV sequences as the representative of the non-ref sequences
    cat {bialsv_seq.fa} {multisv_seq.fa} > {nonref_seq.fa}
    ```

5. Visualize the structural variation in the graph 

    Visualize all paths in the bubbles from start to stop node, require `Graphviz`. 
    
    ```sh
    
    # Visualize all bubbles in the graph (max 500 bubbles) output a single pdf
    visualize/sv_viz.py # output: sv_viz.pdf
    
    # Visualize only bubbles with breakpoint crossing exon 
    visualize/sv_viz_exon.py #output: exon_viz.pdf
    
    ```
