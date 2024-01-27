# Organisation of the GitHub Repository

## Description of Code Used
There are two main scripts written and used for the project: main.py and tools.py. queries.sh is a very short bash script that we only used for making a query (through Unix commands) about our aligned sequences. grapfy.R and most_common_3_csv_maker.py were used to generate the plots helpful to us in understanding the variants. Every piece of code was written with readability in mind, so it should be fairly straightforward to understand what each of them does.

### tools.py
tools.py contains all the functions that we used to process and interpret the variants of unknown significance as well as libraries and dictionaries required to run them - it is our toolkit. It also contains some unused functions and dictionaries that were commented out, as it turned out that we did not need to use them (genetic_code, translator(), SPcalc()).

### main.py
main.py is our main piece of code, where we called the functions in tools.py to reach conclusions about the variants. Naturally, it also specifies the reference sequence to be used, the variants of unknown significance, names of the dictionaries to be built, files to write on or read from etc.

### grapfy.R
grapfy.R was used to generate final_graph.pdf.

### most_common_3_csv_maker.py
most_common_3_csv_maker.py uses information provided by tools.variant_interpreter() to build a .csv file that was later used to generate final_graph.pdf.

## Description of Files Used

### cftr_seq.fa
cftr_seq.fa contains the reference CFTR sequence used by our database (LOVD).

### 800_aligned_ref_sequences_MAFFT_MUT1_aligned.fasta
800_aligned_ref_sequences_MAFFT_MUT1_aligned.fasta contains the sequences in "800_aligned_ref_sequences_MAFFT.fa" with only the addition of Mutation 1 in the mutated_sequences.fa file.

### 800_aligned_ref_sequences_MAFFT_MUT20_aligned.fasta
800_aligned_ref_sequences_MAFFT_MUT20_aligned.fasta contains the sequences in "800_aligned_ref_sequences_MAFFT.fa" with only the addition of Mutation 1 in the mutated_sequences.fa file.

### 800_aligned_ref_sequences_MAFFT.fa
800_aligned_ref_sequences_PSIBLAST.fa contains 800 reference CFTR sequences from different species obtained from PSIBLAST and aligned via MAFFT.

### 1000_aligned_ref_sequences_MUSCLE.fa
1000_aligned_ref_sequences_MUSCLE.fa contains 1000 reference CFTR sequences from different species obtained from BLASTP and aligned via MUSCLE.

### mutated_sequences.fa
mutated_sequences.fa contains the full CFTR sequences (of course, just with the variant position different from that in the reference sequence) of the patients.

### proposal.md
proposal.md is our project proposal written earlier.

### report.md
report.md is our project report.

### results_800.txt
results_800.txt is the results of our automated analysis (through tools.vus_interpreter()) of the VUSes in the context of 800 sequences aligned with MAFFT.

### results_800_revised.txt
results_800_revised.txt is the revised version of results_800.txt.

### results_1000.txt
results_1000.txt is the results of our automated analysis (through tools.vus_interpreter()) of the VUSes in the context of 1000 sequences aligned with MUSCLE.

### results_1000_revised.txt
results_1000_revised.txt is the revised version of results_1000.txt.

### most_c_3.csv
most_c_3.csv is a comma-separated file that contains in Col. 1 a position, in Col. 2 one of the most common amino acids in that position, and in Col. 3 the percentage covered by that amino acid in that position.

### final_graph.pdf
final_graph.pdf displays a graph showing the most common 3 amino acids and the percentages they cover at each aligned position where a VUS was present.

### HGVS
The HGVS folder contains informative documents detailing the standards by which variants are named.

### Trees
The Trees folder contains the trees generated for the project and files associated with each tree.

### Figures
The Figures folder contains the figures used in the report.

### presentation.pdf
presentation.pdf is the presentation for the project.