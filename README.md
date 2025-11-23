# comparative_genomics

because our genomes were scaffolded to a c. vicina genome, we can't really use them to look at genome structure (because itll be biased by being mapped to the c. vicina genome). thus we are mostly going to look at gene family expansion/contraction etc. 

First we will use orthofinder to identify protein sequences from multiple species and figures out how their genes relate to each other.

🔹 Key Steps Inside OrthoFinder

All-vs-All Sequence Search -Compares every protein against every other protein (using DIAMOND/BLAST). -Detects which genes are “similar enough” to likely be homologs. -Clustering into Orthogroups (Gene Families) -Groups genes into “orthogroups,” which are sets of genes descended from a single ancestral gene in the last common ancestor of your species. -These orthogroups contain both orthologs (between species) and paralogs (within species).

Multiple Sequence Alignments & Gene Trees -For each orthogroup, it builds an alignment and then a phylogenetic tree (FastTree by default, ML if you ask). -This tells you how the gene family evolved (duplications, losses, species-specific expansions).

Species Tree Inference -Uses single-copy orthologs to build a species tree. -This tree is rooted and reconciled with the gene trees.

Comparative Genomics Statistics -Summarises gene family sizes per species. -Shows how many genes are single-copy, duplicated, or lost.

download musca domestica and D. melangaster sequences from genbank, use the primary_transcript.py script that comes with Ortho to get only the longest isoform. actually there is a really good ortho tutorial here: https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html
make sure that you have all of the protein files from BRAKER in fasta format plus the two outgroups in a folder like 01_hilli.aa 02_quad.aa 03_stygia.aa 04_vicina.aa 05_musca.aa 06_melangastar.aa

then run ortho on it


# I want to figure out what genes to focus on

The R script attached to this repository processes our OrthoFinder output and generates stacked barplots showing how genes are distributed into different orthogroup categories for each of the four blowfly species 

It has four major stages

1. read and cleans OrthoFinder files
   The input files are Orthogroups.GeneCount.tsv (counts of genes in each orthogroup per species) and Orthogroups_UnassignedGenes.tsv (genes that Orthofinder could not assign to any orthogroup)
- Orthofinder reads the genecount file, removes the summary rows, detects which columns correspond to species, converts everything to numeric, selects all species columns unless focal_species is specified

2. Determine orthogroup categories
   - for each orthogroup, the script decides whether it belongs to a) a 1:1:1 single copy orthologue in all species (i.e., all species has exactly 1 gene in that orthogroup, no duplications, present in all species); b) N:N:N, shared, but not single copy (i.e., present in >2 species. At least one species has >1 copy (multicopy or duplicated) but not in the 1.1.1 category; c) species specific (i.e., orthogroup is present in exactly 1 species - that species gets all the gene counts for that orthogroup); d) unassigned genes 

3. combines everything into 1 data frame, renames everything, creates normal stacked plot with each of the four categories




   
