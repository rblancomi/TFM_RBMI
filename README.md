# Final Master Project

**Master in Bioinformatics and Computational Biology, Universidad Autónoma de Madrid (UAM)**\
**Barcelona Biomedical Genomics Lab (BBGLab) https://bbglab.irbbarcelona.org/** \
**Institut de Recerca Biomèdica de Barcelona (IRB Barcelona)**

Necessary code to reproduce all the data from the final master project entitled as:\
*"The impact and function of degron disruptions in c-terminal truncated proteins and its role in tumorigenesis"*
<pre>
Student:        Raquel Blanco Martínez-Illescas
Directors:    Mònica Sánchez Guixé and Núria López-Bigas
Academic tutor: Luis del Peso Ovalle
Course:         2021/2022
</pre>

## Content

* [Data collection and preprocessing](#data-collection-and-preprocessing)
* [Position Weight Matrices (PWMs) degron motifs](#position-weight-matrices-pwms-degron-motifs)
* [PWM scan](#pwm-scan)
* [PWM positivity threshold, information content, specificity and discovery activity](#pwm-positivity-threshold-information-content-specificity-and-discovery-activity)
* [PWM iterative enrichment](#pwm-iterative-enrichment)
* [De novo degron identification](#de-novo-degron-identification)
* [Mutations in the last exon](#mutations-in-the-last-exon)
* [Mutation annotation in the discovered degrons](#mutation-annotation-in-the-discovered-degrons)
* [Data analysis](#data-analysis)

## Data collection and preprocessing

### 1. UbiNet Position Probability Matrices (PPMs)
>ubinet_PWMs.ipynb

The first part of this Jupyter notebook contains the code to parse the HTML of UbiNet 2.0 database and retrieve the PPMs of each degron motif.

### 2. ELM-Manual database
>create_elm_manual_database.ipynb

Jupyter notebook with the code to generate the ELM-Manual database of experimentally validated degrons.

### 3. Human proteome of Ensembl 92 canonical transcripts
>ensembl_proteome.ipynb

Jupyter notebook with the code to preprocess the downloaded human proteome of Ensembl 92 canonical transcripts.

### 4. Last exons of Ensembl 92 canonical transcripts
>ensembl_last_exons.ipynb

Jupyter notebook with the code to preprocess the downloaded exons of Ensembl 92 canonical transcripts and extract every gene's last exon.

### 5. Annotate Ensembl transcript stable ID (ENST) in CCLE and CPTAC datasets
>stabch_annotate_enst.py

Python script to annotate every mutation or WT form with the ENST of the canonical transcript.

## Position Weight Matrices (PWMs) degron motifs

### 1. PWMs from UbiNet degron motifs
>ubinet_PWMs.ipynb

The second part of this Jupyter notebook contains the code to transform the PPMs into PWMs.

### 2. PWMs from ELM-Manual degron motifs
>elm_manual_PWMs.ipynb

Jupyter notebook with the code to align degron sequences per motif and transform curated alignments into PWMs.

## PWM scan

### 1. Scan
>motifs_scan_proteome.py

Python script to scan a set of proteins (*e.g.*: proteome) with a PWM using a sliding window technique.

### 2. True degrons scan
>motifs_separate_substrates.py

Python script to divide the `motifs_scan_proteome.py` output in true degrons and the rest of proteins.

## PWM positivity threshold, information content, specificity and discovery activity
>motifs_quality_analysis.py

Python script to calculate per motif features and quality metrics.  

## PWM iterative enrichment
>motifs_iterative_enrichment_degener.py

Python script to enrich ELM-Manual PWMs using E3 ligase-substrate interactions from UbiNet database.

## De novo degron identification

### 1. Discovered degrons
>motifs_discovered_degrons.py

Python script to extract the discovered degrons from the `motifs_scan_proteome.py` output. 

### 2. Pool overlapping degrons
>pool_overlapping_degrons.ipynb

Jupyter notebook with the code to define the overlapping discovered degrons and pool them together

## Mutations in the last exon
>stabch_annotate_lastexon.py

Python script to annotate every last-exon mutation in CCLE and CPTAC datasets.

## Mutation annotation in the discovered degrons
>stabch_annotate_degrons.py

Python script to annotate CPTAC and CCLE mutations and WT forms in the discovered degrons. 

## Data analysis
>Figures.ipynb

Jupyter notebook with the code to reproduce all the plots and statistical analysis in the figures of the manuscript. 