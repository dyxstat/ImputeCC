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

Then you can opt to move the files ***Normalized_contact_matrix.npz*** and ***contig_info.csv*** to the directory of the ImputeCC folder.


# Usage
## Implement the ImputeCC pipeline
The ImputeCC binning pipeline includes two main steps, i,e, the imputation step and the clustering step.
You can use module `pipeline` to run both steps in succession:
```
python ./ImputeCC.py pipeline [Parameters] FASTA_file CONTIG_INFO HIC_MATRIX OUTPUT_directory
```
### Parameters
```
--rwr-rp: Restarting probability for CRWR (default 0.5)
--rwr-thres: Percentile threshold, determining that Hi-C contacts falling below this threshold
              will be discarded at each random walk step (default 80)
--intra: percentile threshold to assign the contigs to preliminary bins in the preclustering step (default 50)
--inter: percentile threshold to assign the contigs to new bins in the preclustering step (default 0)
--gene-cov: Gene coverage used to detect marker genes (default 0.9)
--max-markers: The maximum number of contigs with marker genes
--cont-weight: Coefficient of completeness - cont_weight * completeness (default 2)
--min-comp: Minimum completeness of bin to consider during the final integrative binning step (default 50)
--max-cont: Maximum contamination of bin to consider during the final integrative binning step (default 10)
--report-quality: Minimum quality of bin to report (default 10)
--min-binsize: Minimum bin size used in the output (default 100,000)
-t/--threads: The number of threads (default 30)
--cover: Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v: Verbose output about more specific details of the procedure.
```


### Input File

* **FASTA_file**: a fasta file of the assembled contigs (e.g. Test/final.contigs.fa)
* **CONTIG_INFO**: a bam file of the Hi-C alignment (e.g. Test/MAP_SORTED.bam)
* **HIC_matrix**: a bam file of the Hi-C alignment (e.g. Test/MAP_SORTED.bam)
* **OUTPUT_directory**: a bam file of the Hi-C alignment (e.g. Test/MAP_SORTED.bam)

### Output File

* **contig_info.csv**: information of assembled contigs with three columns (contig name, the number of restriction sites on contigs, and contig length).
* **Normalized_contact_matrix.npz**: a sparse matrix of normalized Hi-C contact maps in csr format and can be reloaded using Python command *'scipy.sparse.load_npz('Normalized_contact_matrix.npz')'*.
* **NormCC_normalized_contact.gz**: Compressed format of the normalized contacts and contig information by pickle. 
This file can further serve as the input of MetaCC binning module.
* **MetaCC.log**: the specific implementation information of NormCC normalization module.


### Example
```
python ./MetaCC.py norm -e HindIII -e NcoI -v final.contigs.fa MAP_SORTED.bam out_directory
```


## Implement the MetaCC binning module
**MetaCC binning module is based on the NormCC-normalized Hi-C contacts and thus must be implemented after the NormCC normalization module.**
```
python ./MetaCC.py bin --cover [Parameters] FASTA_file OUTPUT_directory
```
### Parameters
```
--min-binsize: Minimum bin size used in output (default 150,000)
--num-gene (optional): Number of marker genes detected. If there is no input, 
                       the number of marker genes will be automatically detected.
--random-seed (optional): seed for the Leiden clustering. If there is no input, a random seed will be employed.
-v (optional): Verbose output about more specific details of the procedure.
```
### Input File

* **FASTA_file**: a fasta file of the assembled contigs (e.g. Test/final.contigs.fa)
* **OUTPUT_directory**: please make sure that the output directory of the MetaCC binning module should be the same as that of the NormCC normalization module.

### Output File

* **BIN**: folder containing the fasta files of draft genomic bins
* **MetaCC.log**: the specific implementation information of MetaCC binning module


### Example
```
python ./MetaCC.py bin --cover -v final.contigs.fa out_directory
```


## Implement the post-processing step of the MetaCC binning module
Draft genomic bins are assessed using [CheckM2](https://github.com/chklovski/CheckM2)/[CheckM](https://github.com/Ecogenomics/CheckM).
Then the post-processing step of the MetaCC binning module is conducted for partially contaminated bins with completeness larger than 50% and contamination larger than 10% in order to purify the partially contaminated bins. 
```
python ./MetaCC.py postprocess --cover [Parameters] FASTA_file Contaminated_Bins_file OUTPUT_directory
```

### Parameters
```
--min-binsize: Minimum bin size used in output (default 150,000)
-v (optional): Verbose output about more specific details of the procedure.
```

### Input File
* **FASTA_file**: a fasta file of the assembled contigs (e.g. Test/final.contigs.fa).
* **Contaminated_Bins_file**: a csv file of the names of the partially contaminated bins; Bin names are arranged in columns and don't include the file formats .fa at the end of each name.
* **OUTPUT_directory**: please make sure that the output directory of the post-processing step of the MetaCC binning module should be the same as the previous steps.

Example of a Contaminated_Bins_file:
```
BIN0001
BIN0003
BIN0005
...
```

### Output File

* **BIN**: folder containing the fasta files of draft genomic bins after cleaning partially contaminated bins.
* **MetaCC.log**: the specific implementation information of the post-processing step of the MetaCC binning module.


### Example
```
python ./MetaCC.py postprocess --cover -v final.contigs.fa contaminated_bins.csv out_directory
```



# Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


# Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.







