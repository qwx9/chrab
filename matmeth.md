Materials and methods
=====================

16/12/2019


Data
----

Because of its stability and the large amount of available data, the Human Umbilical Vein Endothelial Cell (HUVEC) line was the sole focus for this initial study.
Also due to available data from preliminary studies, hg19 was chosen as a reference genome.

Data on histone markers was gathered from validated ChIP-seq and DNase-seq experiments in the HUVEC reference epigenome series of the ENCODE project (ref. ENCSR194DQD), namely for H3K4me3 (ENCSR578QSO), H3K27ac (ENCSR000ALB) and DHS (ENCSR000EOQ).
A recent study on HUVEC (GEO GSE94872) with GRO-seq data was used to quantify expression intensity along the genome.  The samples selected were the 4 normoxia tracks (GSM2486801 to GSM2486804).
In addition, a list of potential silencers was taken from supplemental data in another recent paper (Ovcharenko2019, Supp. Table S2).

The latest data on various characteristics of hg19 was taken from the UCSC FTP server, such as RepeatMasker output, NCBI Refseq genes, GC content, chromosome size and the ChromHMM track for HUVEC.

A preliminary study yielded the A-B profile for several cell lines, including HUVEC, providing a consensus AlwaysA/AlwaysB list of regions, as well as a list of repeat sequences (RepSeq) along hg19 (from RepBase), and their correlation to HUVEC's computed A-B profile.
A list of strong protosilencers was built from several prior publications.


Pre-processing and counts
-------------------------

In order to prepare the data for counting and modelisation, several of the input data files had to be transformed.

To build a gene counts along hg19, the RefSeq list was not taken as is, since it may contain nomenclature errors and many overlapping alternative genes.  First, genes of the same symbolic name on the same chromosome and strand were merged together.  Overlapping genes on the same chromosome and strand which shared a common prefix of at least one letter were then also merged together.  This avoided scenarios where many overlapping alternative definitions of the same gene, depending on start and end coordinate, would falsely inflate the counts.  Here, merging is defined as summarizing several intervals by a single interval with a starting position as the minimum start of its components, and an end position as the maximum end of its components.

GC content on hg19 is provided as a percentage over 5 nucleotides.  However, counts are done on hits along 100kb windows.  A mean GC content percentage was defined for every 100kb bin along each chromosome.

GRO-seq data was provided as 4 files concatenating peaks over each strand.  The GRO-seq signals for all replicates was first concatenated together, preserving the strand information.  A mean of the overlapping signals was then computed.  Then, a mean for each 100kb bin was computed as well.  Finally, outliers according to the Tukey's Fences method with k=5 upper bound were brought down to this threshold.  This served as the expression intensity parameter in linear modeling.  For counts, from the mean of overlapping signals computed before, strictly adjacent intervals (distance of 1 nucleotide) were merged together and were later used for counts.

RepSeq data for a number of families and combination of families was extracted from the hg19 RepeatMasker output based on various filters taking into account correlation thresholds for A and B regions from the results of the preliminary study.  In particular, pro-A elements were defined as having a mean correlation with A equal or greater to 0.01, and pro-B elements as those with a mean correlation with A lower than -0.01.  Elements were then further split into transposable elements (SINEs, RC, SVA, LTR, LINEs, DNA) and non-transposable elements (tRNA, snRNA, simple repeats, scRNA, satellites, rRNA and low complexity regions).  LINEs were also split to process L1s and long (>5kb) L1s separately.

The HUVEC ChromHMM track was split into individual files for each class.  Since most strong promoters are surrounded by adjacent weak promoters, weak promoters within 1000 nucleotides of a strong promoter were merged with the strong promoter and removed from the weak promoter list.  Classes 4 and 5 both indicated strong enhancers and were merged together.  The same was done for weak and poised enhancers (classes 6 and 7).

Counting was done with the BEDtools software suite along 100kb windows for each chromosome.

To avoid discrepancies between GRO-seq and ChromHMM data, these were not intersected to provide counts of "active genes".  Instead, HUVEC active promoters were defined as a merge between strong and weak promoters in ChromHMM data.

For visual representations of each element, counts were split in classes based on either an A or B tendency, then again based on gene density (none, regular (<4 genes) and high), then once more by the presence or not of active promoters as defined above.


Multivariate linear model of AB profile
---------------------------------------

To investigate whether a simple score consisting of the weighted sum of a number of counts for each bin would be sufficient to correctly predict a tendency towards A or B, a number of multivariate linear models were built from combinations between lists of elements.

A set of lists for epigenetic parameters were defined, based either on ChromHMM data, and/or silencer lists and DHS data, or histone markers, GRO-seq intensity data and DHS.  A second inherent chromosome characteristics list with various combinations of RepSeq sub-classes correlated to A or B, as well as GC percentage and gene density was added.  Then, all possible combinations of explanatory parameters based on intersections of these two lists were generated.

To avoid the effects of abrupt transitions between mostly A and mostly B regions, so called "flanking" or intermediary regions were removed from the eigenvector prior to prediction.  They were defined as the 200kb surrounding any sign transition of the HUVEC eigenvector value, i.e. any 100kb bin showing a sign change plus the next and the two bins before it.

A Linear model was built for each combination of explanatory variables, with the HUVEC eigenvector without flanking regions as a dependent variable.  Models were then compared using the adjusted square of the sample Pearson correlation coefficient to estimate the fraction of the variance in the eigenvector that is explained by the parameter vector, whilst adjusting for the number of parameters in the model.  Those models with the highest adjusted R-squared value were defined as the most efficient.


Perspectives (discussion)
-------------------------
- higher efficiency general linear models by implementing a correct variable selection algorithm based on AIC and likelihood ratio tests on nested models (rather than preselecting combinations and generating every model)
- logistic regression
- non-linear models
- neural networks and comparison to linear model efficiency
- regression trees for highly correlated parameters


!!! HUVECnoflank vector was built from HUVEC vector, not AorBvec!
!!! AorBvec not used in tabs for classification!
!!! mean null = 0 never set for groseq.mean!
