# ImputeCC enhances integrative Hi-C-based metagenomic binning through constrained random-walk-based imputation

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [A test demo](#a-test-demo)
- [Preparation](#preparation)
- [Usage](#usage)
- [Contacts and bug reports](#contacts-and-bug-reports)
- [Copyright and License Information](#copyright-and-license-information)
- [Issues](https://github.com/dyxstat/ImputeCC/issues)

# Overview
`ImputeCC` is an integrative contig binning tool tailored for metaHi-C datasets. 
ImputeCC integrates Hi-C interactions with the inherent discriminative power of single-copy marker genes, 
initially clustering them as preliminary bins, and develops a new constrained 
random walk with restart (CRWR) algorithm to improve Hi-C connectivity among these contigs.


# System Requirements
## Hardware requirements
`ImputeCC` requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
`ImputeCC` v1.0.0 is supported and tested in *Linux* systems.

### Python Dependencies
`ImputeCC` mainly depends on the Python scientific stack:

```
numpy
pandas
biolib
scipy
scikit-learn
igraph
leidenalg
FragGeneScan
hmmer
checkm
```



# Installation Guide
We recommend using [**conda**](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) to install `ImputeCC`. 
Typical installation time is 1-5 minutes depending on your system.

### Clone the repository with git
```
git clone https://github.com/dyxstat/ImputeCC.git
```

Once complete, enter the repository folder and then create a `ImputeCC` environment using conda.


### Enter the MetaCC folder
```
cd ImputeCC
```


### Construct the conda environment in the linux or MacOS system
```
conda env create -f ImputeCC_env.yaml
```

### Enter the conda environment
```
conda activate ImputeCC_env
```


# A test demo
To test the software, please use
```
python ./ImputeCC.py test
```


# Preparation
Follow the instructions in this section to generate the input for `ImputeCC`:

### Preprocess raw metagenomic Hi-C libraries
For the shotgun library, de novo metagenome assembly is produced by an assembly software, such as MEGAHIT.
```
megahit -1 SG1.fastq.gz -2 SG2.fastq.gz -o ASSEMBLY --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95
```
Hi-C paired-end reads are aligned to assembled contigs using a DNA mapping software, such as BWA MEM. Then, samtools with parameters ‘view -F 0x904’ is applied to remove unmapped reads, supplementary alignments, and secondary alignments. BAM file needs to be sorted **by name** using 'samtools sort'.
```
bwa index final.contigs.fa
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz > MAP.sam
samtools view -F 0x904 -bS MAP.sam > MAP_UNSORTED.bam
samtools sort -n MAP_UNSORTED.bam -o MAP_SORTED.bam
```

### Generating normalized metagenomic Hi-C contact matrix by NormCC
You need to run the `NormCC` normalization module from the [MetaCC](https://github.com/dyxstat/MetaCC) software 
to generate NormCC-normalized Hi-C contact matrix. 
For instance, once you install the MetaCC software, run
```
python MetaCC.py norm -v final.contigs.fa MAP_SORTED.bam out_NormCC
```
The NormCC-normalized Hi-C contact matrix and the corresponding sorted contig information are stored in
the files ***Normalized_contact_matrix.npz*** and ***contig_info.csv*** from the output directory (i.e., out_normcc), respectively.
***Normalized_contact_matrix.npz*** is a sparse matrix of normalized Hi-C contact maps in python scipy sparse csr format, and
***contig_info.csv*** stores the information of assembled contigs with three columns (contig name, the number of restriction sites on contigs, and contig length).

Then you can opt to move the files ***Normalized_contact_matrix.npz*** and ***contig_info.csv*** to the directory of the ImputeCC folder.


# Usage
## Implement the ImputeCC pipeline
The ImputeCC binning pipeline includes two main steps, i,e, the imputation step and the clustering step.
Use the module `pipeline` to run both steps in succession:
```
python ./ImputeCC.py pipeline [Parameters] FASTA_file CONTIG_INFO HIC_MATRIX OUTPUT_directory
```
### Parameters
```
--rwr-rp: Restarting probability for CRWR (default 0.5)
--rwr-thres: Percentile threshold, determining that Hi-C contacts falling below this threshold
              will be discarded from the imputed matrix at each random walk step (default 80)
--intra: percentile threshold to assign the contigs to existing preliminary bins in the preclustering step (default 50)
--inter: percentile threshold to assign the contigs to new preliminary bins in the preclustering step (default 0)
--gene-cov: Gene coverage used to detect marker genes (default 0.9)
--max-markers: The maximum number of marker-gene-containing contigs to construct preliminary bins (default 8000) 
--cont-weight: Coefficient of completeness - cont_weight * completeness (default 2)
--min-comp: Minimum completeness of bin to consider in the final integration (default 50)
--max-cont: Maximum contamination of bin to consider in the final integration (default 10)
--report-quality: Minimum quality of bin to report in the final integration (default 10)
--min-binsize: Minimum bin size used in the output (default 100000)
--threads/-t: The number of threads (default 30)
--cover: Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v: Verbose output about more specific details of the procedure.
```

### Input File

* **FASTA_file**: a fasta file of the assembled contigs (e.g. final.contigs.fa)
* **CONTIG_INFO**: a contig information file generated by NormCC (e.g. contig_info.csv)
* **HIC_MATRIX**: normalized Hi-C contact matrix generated by NormCC (e.g. Normalized_contact_matrix.npz)

### Output File
* **FINAL_BIN**: folder containing the fasta files of draft genomic bins
* **ImputeCC.log**: the specific implementation information of ImputeCC
* **intermediate**: folder containing the intermediate files, including
  * *contigs_marker_genes.hmmout*: single copy marker genes detected from assembled contigs
  * *ImputeCC_storage.gz*: compressed format of all results from the imputation step
  * *checkm_gene_table*: lineage-specific gene tables identified by CheckM


### Example
```
python ./ImputeCC.py pipeline -v final.contigs.fa contig_info.csv Normalized_contact_matrix.npz out_directory
```


## Implement the imputation step
Use the module `impute` to only run the imputation step:
```
python ./ImputeCC.py impute --cover [Parameters] FASTA_file CONTIG_INFO HIC_MATRIX OUTPUT_directory
```
### Parameters
```
--rwr-rp: Restarting probability for CRWR (default 0.5)
--rwr-thres: Percentile threshold, determining that Hi-C contacts falling below this threshold
              will be discarded from the imputed matrix at each random walk step (default 80)
--gene-cov: Gene coverage used to detect marker genes (default 0.9)
--max-markers: The maximum number of marker-gene-containing contigs (default 8000) 
--threads/-t: The number of threads (default 30)
--cover: Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v: Verbose output about more specific details of the procedure.
```
### Input File

* **FASTA_file**: a fasta file of the assembled contigs (e.g. final.contigs.fa)
* **CONTIG_INFO**: a contig information file generated by NormCC (e.g. contig_info.csv)
* **HIC_MATRIX**: normalized Hi-C contact matrix generated by NormCC (e.g. Normalized_contact_matrix.npz)

### Output File
* **ImputeCC.log**: the specific implementation information of ImputeCC
* **intermediate**: folder containing the intermediate files, including
  * *contigs_marker_genes.hmmout*: single copy marker genes detected from assembled contigs
  * *ImputeCC_storage.gz*: compressed format of all results from the imputation step

### Example
```
python ./ImputeCC.py impute --cover -v final.contigs.fa contig_info.csv Normalized_contact_matrix.npz out_directory
```


## Implement the clustering step
Use the module `cluster` to only run the imputation step
**(the clustering step utilizes the imputed Hi-C matrix and thus must be implemented after the imputation step)**:
```
python ./ImputeCC.py cluster --cover [Parameters] FASTA_file CONTIG_INFO HIC_MATRIX OUTPUT_directory
```
### Parameters
```
--intra: percentile threshold to assign the contigs to existing preliminary bins in the preclustering step (default 50)
--inter: percentile threshold to assign the contigs to new preliminary bins in the preclustering step (default 0)
--cont-weight: Coefficient of completeness - cont_weight * completeness (default 2)
--min-comp: Minimum completeness of bin to consider in the final integration (default 50)
--max-cont: Maximum contamination of bin to consider in the final integration (default 10)
--report-quality: Minimum quality of bin to report in the final integration (default 10)
--min-binsize: Minimum bin size used in the output (default 100000)
--threads/-t: The number of threads (default 30)
--cover: Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v: Verbose output about more specific details of the procedure.
```

### Input File

* **FASTA_file**: a fasta file of the assembled contigs (e.g. final.contigs.fa)
* **CONTIG_INFO**: a contig information file generated by NormCC (e.g. contig_info.csv)
* **HIC_MATRIX**: normalized Hi-C contact matrix generated by NormCC (e.g. Normalized_contact_matrix.npz)

### Output File
* **FINAL_BIN**: folder containing the fasta files of draft genomic bins
* **ImputeCC.log**: the specific implementation information of ImputeCC
* **intermediate**: folder containing the intermediate files, including
  * *checkm_gene_table*: lineage-specific gene tables identified by CheckM


### Example
```
python ./ImputeCC.py cluster --cover -v final.contigs.fa contig_info.csv Normalized_contact_matrix.npz out_directory
```



# Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


# Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.







