Tumour-endothelial_Contact-initiated_Phospho-signalling
=======================================================

Phospho-Proteomic analysis of tumour-endothelial bidirectional contact-initiated signalling


AIM: 

The aim of the project is to identify regulated phosphosites in breast
cancer cells (MDA-MB-231-4175-TGL, referred to as LM2 or MDA depending on the
document) and endothelial cells (HUVEC) when they interact.

EXPERIMENTAL OUTLINE:

In 15 cm petry dishes, monolayers of HUVEC were grown to confluency. After 1h30
of HUVEC starvation, twice the number of MDA-MB-231 (starved overnight
and detached with enzyme-free cell dissociation buffer) were plated onto the
HUVECs. Co-culture was conduced for 15 min at 37 degrees in a small volume
of serum free medium (3ml per 15cm dish). 
Then, non-attached LM2 were removed with the medium and counted to assess the
number of cells still attached.
SILAC labelling was used to differentiate cancer cells from endothelial cells 
and perform the quantification. Heavy label was systematically used for cells in co-
culture; medium label was used for cells in mono-culture. Non-labelled cells were 
the other cell type (not taken into account in the analysis). 
Heavy-labelled and medium-labelled cells were mixed in a 1:1 cell ratio. The 
experiment was performed 5 independent times with labelled LM2 and 4 
independent times with labelled HUVECs. 
Cell lysates were subjected to a rough membrane fragmentation, then membrane-
enriched and cytoplasmic fractions were prepared in parallel (phospho-
enrichment with TiO2, lys-C / trypsin digestion) before MS/MS. Raw data were
searched against SwissProt database using PD1.3 / Mascot. Tables were extracted
and for each identified peptide the position of the first amino-acid in the 
protein was used using the same database as the one utilised for the search. 

DATA:

- "Samples.csv" contains the sample descriptions (experiment number, labelled
cell type, fraction, sample name and raw file name).
- "SILAC1018TotV5_FAAposition.csv" is the raw table containing all the peptides
identified in the experiment and the position of the first amino-acid in the 
protein sequence.
- “HUVsignV5spF.csv", "MDAsignV5spF.csv" and "MDAsignV5spFLMNA.csv" contain the spectra corresponding to the list of regulated targets after manual inspection.  
- "HUVChecked.csv" and "MDAChecked.csv" are the equivalent tables after manual inspection of a broader range of phosphosites of interest.
- "DA-MDA-HUVEC" is the markdown document associated with the data analysis.
- “MDARegProt.csv” is the list of uniprot entry names associated to the proteins in the “Pathway” figures (manually assessed).
- “HUVRegProt.csv” is the equivalent for HUVEC regulated proteins.
- “TableRegPsites” and related document describe the construction of the
tables with median of regulation fold per experiment.

In the folder "Proteomics-20151015" are all the documents regarding proteomics analysis:
MS data from the Proteome Discoverer search are in the folders "Analysis_2" and "EPHA2-IP". 
More precised description is provided in the "TotProtAnalysis_Final" files. All the 
tables and figures from the analysis are in the folders "OutputTab" and "OutputFig".

DESCRIPTION:

The table "Samples.csv" contains the description of all samples analysed for
the project.
- experiment:
each experiment is labelled with "SILAC" and a number (from 010 to 019 - with 
no SILAC018). 
- sample:
each sample is labelled with "ps" of "psMLP" and a number. For each experiment,
there are at least two samples prepared (membrane-enriched fraction and 
cytoplasmic fraction). There can be several sample preparations / technical 
replicates, injected different days (thus several "SpectrumFile").

