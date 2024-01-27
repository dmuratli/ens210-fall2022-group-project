# An Inquiry About the Relationship between CFTR Variants of Unknown Significance and Cystic Fibrosis

## Authors: Deniz Muratli, Salih Kaya, A. Bera Efe

## Introduction
Cystic fibrosis is an autosomal recessive genetic disease that results from mutations in the CFTR gene, which encodes a protein that is involved in the transport of mucus, sweat, saliva, tears, and digestive enzymes through the membrane (“CFTR gene”, n.d.; “Cystic fibrosis,” 2022). Such mutation results in thick and sticky mucus that affects various organs in the human body, such as the lungs, pancreas, liver, spleen, stomach, bones, and joints (“Cystic fibrosis,” 2022).

Benign, pathogenic, and unclassified variants of CFTR are well-documented, although prone to reclassifications as a result of new genome-wide association studies (like the variants of most genes). The classification procedure is often lengthy, involving "the annotation of variants, application of frequency filters and database searches to enrich for rare variants and eliminate common variants, and prediction of functional effect." The clinical evaluation of a variant also requires "an assessment of potential effects on the function of one or more genes and the evidence supporting the attribution of the illness at presentation to the affected gene or genes." Several metrics are used for the process (Adams & Eng, 2018, p. 1355). As such, in practice, it is difficult to classify variats into a category or another on simple criteria.

Nevertheless, this work is an attempt at determining whether VUSes (variants of unknown significance) in a publicly available database are pathogenic using few but informative criteria. VUSes used were further restricted to point substitutions. The authors do not claim to have decisively determined the effect of the examined variants, as that would require significantly more experience and work - predictions provided in this paper are merely educated guesses.

## Results
Predictions of the program for 72 variants of unknown significance (more than half of those available on LOVD) can be viewed in results_800_revised.txt and results_1000_revised.txt, with the difference between the files being that they are based on differently obtained and aligned sequences. Only four of those from results_800_revised.txt will be discussed here due to time and space limitations.

An overview of most common characters (amino acids) at each of the 72 VUS positions is summarised in **Figure 1**, and a phylogenetic tree formed with aligned [healthy] CFTR sequences can be seen in **Figure 2**.

![The most common 3 amino acids and the percentages they cover at each aligned position where a VUS is present.](/Figures/figure_1.png "The most common 3 amino acids and the percentages they cover at each aligned position where a VUS is present.")
**Figure 1**: The most common 3 amino acids and the percentages they cover at each aligned position where a VUS is present.

![Phylogenetic tree of 800 reference CFTR sequences.](/Figures/figure_2.png "Phylogenetic tree of 800 reference CFTR sequences.")
**Figure 2**: Phylogenetic tree of 800 reference CFTR sequences.

### VUS 1
VUS 1 is found at position 3290 of the aligned sequence. The program reports that the most common amino acid at this position is aspartic acid, making up 70% of the amino acids at the position. Glutamic acid follows it with 26%, making a remarkable 96% of the amino acids at the position acidic in chemistry, suggesting a significance for this property.

The variant causes asparagine to replace aspartic acid, which stands out from the remaining amino acids in the position by the fact that it is a polar uncharged amino acid. Although the change in entropy this substitution causes is small in value (as will be observed in other variants too) due to the fact that this sequence only makes up a single data point among 800 sequences, it is high in comparison to those of other listed substitutions. This means that it disrupts the conservation pattern in the position. Combined with the fact that it is chemically dissimilar to the original amino acid, it may suggest pathogenicity.

However, it must be noted that there is still some probability that the substitution is benign, as the BLOSUM62 score of the substitution is not negative and the substitution has not occurred in a critical part of the protein. This observation is also reflected in the automated prediction of the algorithm. The tree in **Figure 3** contains the variant, and, as expected, it is most closely associated with the reference human sequence.

![Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 1.](/Figures/figure_3.png "Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 1.")
**Figure 3**: Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 1.

### VUS 2
VUS 2 is found at the position 2896 of the aligned sequence. Amino acid having the highest frequency at this position is glutamic acid, constituting 95% of profile here. It is followed by asparagine with 3%; and in total, these two amino acids make up 89% of this position. The two have the common feature of being acidic amino acids. Consequently, it can be concluded that being acidic is a significant part of this position’s nature.

The variant results in the replacement of glutamic acid with lysine; that is, being a positively charged amino acid, in comparison to the highly dominant glutamic acid being a negatively charged amino acid. The change in entropy is relatively high, but lower than that of VUS 1. Similarly to VUS 1, it has a positive BLOSUM62 substitution score and is not in a functionally significant part, so it seems apt to have some confidence in it being non-pathogenic. Again, the algorithm agrees with the observation.

### VUS 3
VUS 3 is found at the 3229th position of the aligned sequence. Methionine is the most common amino acid at this position with 60% - not a decisive majority. Leucine covers 32%, isoleucine covers ~7%, and valine covers 1%. Therefore, practically all (100%) of the amino acids at the position are aliphatic - most probably a significant feature.

The variant causes the amino acid at this position to become threonine, a polar uncharged amino acid. Even though it surprisingly causes a relatively small change in entropy and is not found in a critical part of the protein, this mismatch has a BLOSUM62 score of -1. Considering that it is at odds with the trademark chemical feature of the position, it is fairly certain that VUS 3 is pathogenic. This is also acknowledged by the algorithm.

### VUS 20
VUS 20 is found at position 2994 of the aligned sequence. The amino acid composition of this position is heterogeneous, with threonine covering 45%, isoleucine covering 27%, glycine covering 11%, and valine covering 8%. As a result, 45% of the amino acids here are polar uncharged, whereas 46% are aliphatic.

The reference amino acid (threonine) is substituted with isoleucine, an amino acid chemically similar to most (by bare majority) amino acids at this position. Considering that 27% of the species also have isoleucine at this position and the change in entropy is negligible, this substitution initially seems benign but is marked by the program as “Likely pathogenic”. This is despite the fact that it is not in a critical section of the protein - the algorithm is not lenient enough to tolerate a substitution that is chemically different from the most common characters and has a BLOSUM62 score of -1 (indicating that the substitution is statistically unlikely to occur in nature). 

It is not certain that the program rightly classified the variant, but the authors tend to trust the BLOSUM62 score more because it is based on statistical facts. **Figure 4** displays the variant sequence along with the reference sequences and the result is similar to that of VUS 1.

![Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 20.](/Figures/figure_4.png "Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 20.")
**Figure 4**: Phylogenetic tree of 800 sequences and the human reference CFTR sequence mutated with VUS 20.

## Discussion
The approach taken in this paper was certainly simplistic and confidence in the classifications could be much higher if data showing that the existence of a given VUS commonly corresponds to illness was available. However, it is believed that considering whether the mutation took place in a critical part of the protein and treating it accordingly was a reasonable approach, as changes in functionally significant sections will inevitably disturb the relevant function of the protein. As evident from the algorithm, the BLOSUM62 substitution score obtained was the ultimate determining factor in the automated classification of a variant - a decision the authors are confident in, since BLOSUM62 scores are based on statistics and hence are more likely to reflect the ground truth. Other helpful approaches might have included included calculating changes in identity score as substitutions occurred or taking only the most common amino acids of close relatives of *Homo sapiens* in that position into consideration.

Because of the nature of the substitutions inspected, the trees did not yield helpful clues to interpreting the variants. In the phylogenetic tree, both M1 and M20 were placed most closely to the human reference sequence. This was expected: since each mutated sequence differs from the reference sequence by only one amino acid but is otherwise the same as the reference, the closest relative was always the human CFTR sequence. Using the trees might have helped if the case was that there were several substitutions; but then it would also be easier to predict the pathogenic effect.

## Materials and Methods
The aligned sequences were used to calculate changes in Shannon's entropy due to the mutations and most common characters and associated chemical properties in substituted positions. Substitutions were further evaluated based on their BLOSUM62 score when paired with the original amino acid. Based on these pieces of information (excluding entropy, that was only meant for manual examination), the program makes an automated prediction on the likelihood of pathogenicity.

The decision-making algorithm is as follows (if statements were executed in this order):
- If the position is in a critical part of the protein (binding site, domain, motif, defined region etc.):
    - If the reference amino acid is the same as the new one, the variant is likely *non-pathogenic*.
    - If the BLOSUM62 score of the substitution is bigger than or equal to 0 and the chemical group of the reference of the amino acid is the same as that of the new one, the variant is *likely non-pathogenic*.
    - Otherwise, the variant is *likely pathogenic*.
- If the position is not in a critical part of the protein:
    - If the reference amino acid is the same as the new one, the variant is likely *non-pathogenic*.
    - If the BLOSUM62 score of the substitution is bigger than or equal to 0 and the chemical group of the reference of the amino acid is the same as that of the new one, the variant is *likely non-pathogenic*.
    - If the BLOSUM62 score of the substitution is bigger than or equal to 0 and the new amino acid is the among the most common 2 amino acids in that position, the variant is *likely non-pathogenic*.
    - If the BLOSUM62 score of the substitution is bigger than or equal to 0 and the chemical property of the new amino acid is the among the 2 most common chemical properties in that position, the variant is *likely non-pathogenic*.
    - If the BLOSUM62 score of the substitution is bigger than or equal to 0, the variant is *possibly (in lower likeliness than "likely") non-pathogenic*.
    - Otherwise, the variant is *likely pathogenic*.

Information about the segmentation of the protein (its domains, regions etc.) were obtained from UniProt and manually inputted as a matrix to the algorithm (n.d.). As evident from the flow, any deviations from the reference amino acid and its chemical group was strongly penalised when the position was within a critical section of the protein, directly classified as "likely pathogenic." On the other hand, substitutions in other sections were treated more leniently, with more drastic deviations still being classified as "likely non-pathogenic."

For the trees, approximately 800 sequences with the highest identity were obtained from PSIBLAST. These were then aligned using MAFFT alignment tool. The aligned fasta file was converted to .meg format using MEGA and a tree was created with the Maximum Parsimony method with a bootstrap number of 1000. Following that, the tree was saved in Newick format. Rectangular and polar trees were created from the raw data in Newick format, as well as from the Raw Data+M1 and Raw Data+M20 fasta files.

The plot displaying the most common 3 amino acids and the percentages they cover at each aligned position where a VUS was present was plotted via R, based on transformed VUS data.

The specific scripts and programming languages used are detailed in method.md, which is also in the GitHub repository.

## References
Adams, D. R., & Eng, C. M. (2018). Next-generation sequencing to diagnose suspected genetic disorders. *New England Journal of Medicine, 379*(14), 1353-1362. https://doi.org/10.1056/nejmra1711801

*CFTR gene*. (n.d.). MedlinePlus - Health Information from the National Library of Medicine. Retrieved December 15, 2022, from https://medlineplus.gov/genetics/gene/cftr/

Claustres, M., Thèze, C., Des Georges, M., Baux, D., Girodon, E., Bienvenu, T., Audrezet, M., Dugueperoux, I., Férec, C., Lalau, G., Pagin, A., Kitzis, A., Thoreau, V., Gaston, V., Bieth, E., Malinge, M., Reboul, M., Fergelot, P., Lemonnier, L., … Bareil, C. (2017). CFTR-France, a national relational patient database for sharing genetic and phenotypic data associated with rare CFTR variants. *Human Mutation, 38*(10), 1297-1315. https://doi.org/10.1002/humu.23276

*Cystic fibrosis*. (2022, May 9). Centers for Disease Control and Prevention. Retrieved December 15, 2022, from https://www.cdc.gov/genomics/disease/cystic_fibrosis.htm

National Center for Biotechnology Information. *BLAST: Basic local alignment search tool*. (n.d.). https://blast.ncbi.nlm.nih.gov/Blast.cgi

National Center for Biotechnology Information. *Homo sapiens CF transmembrane conductance regulator (CFTR), mRNA*. (n.d.). Retrieved January 1, 2023, from https://www.ncbi.nlm.nih.gov/nuccore/90421312

UniProt. (n.d.). *CFTR_HUMAN*. Retrieved January 7, 2023, from https://www.uniprot.org/uniprotkb/P13569/entry

*Unique variants in the CFTR gene*. (n.d.). Leiden Open Variation Database. Retrieved January 1, 2023, from https://databases.lovd.nl/shared/variants/CFTR/unique#object_id=VariantOnTranscriptUnique%2CVariantOnGenome&id=CFTR&order=VariantOnTranscript%2FDNA%2CASC&search_transcriptid=00000116&search_VariantOnTranscript/DNA=%3E&search_VariantOnGenome/ClinicalClassification=VUS&page_size=1000&page=1