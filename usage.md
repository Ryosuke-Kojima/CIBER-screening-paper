# How to evaluate the change of sEV release upon gene knockout
This document describes the usage of codes (count_barcodes.py and calculate_zRE.py) in our [paper] (https://www.biorxiv.org/content/10.1101/2023.09.28.559700v1) to calculate z-RE, the change of sEV release for a gene upon knockout, from raw .fastq files. 

## System requirements
We used the following packages in Windows 11 Home. Most packages are accompanied by [anaconda3]( https://www.anaconda.com/download) while other packages that is not included in anaconda such as [biopython]( https://biopython.org/) are downloaded from websites or individually pip-installed. Note that the latest or older versions have now possibly become incompatible with our codes.
- python 3.11.4
- conda 23.7.2
- biopython 1.81
- matplotlib 3.7.1
- numpy 1.24.3
- pandas 1.5.3
- scikit-learn 1.3.0

## Installation guide
The installation procedure is very simple. Go to our [Github repository]( https://github.com/Ryosuke-Kojima/CIBER-screening-paper) and download all the files (2 python scripts Åe*.pyÅf, one .gz files ÅeDTKP_CD63_sEVs_Cas9_rep1.fastq.gzÅf, 7 barcode count files Åecount_*.csvÅf and one reference file Åeref.xlsxÅf) in a same directory in your PC. 

## Demo
In our paper, we have performed 4 parallel screening with 4 subpool gRNA libraries (i.e., ACOC, DTKP, PROT, TMMO) and processed data acquired in each screening individually. We demonstrate how we processed the data acquired with DTKP library in this section. For each screening, 8 .fastq files are acquired and processed with count_barcodes.py to count barcodes. Z-RE are calculated via calculate_zRE.py using the output files of count_barcodes.py. 

### Barcode counting
You can run count_barcodes.py in your IDE just by filling line 83 with your directory where Demo_data.fastq.gz is located:
```sh
path = r'your directory'
```
After running, barcode counts will be written to an output count_*.csv file with the same name as .fastq file which is able to be used for z-RE calculation without formatting. Relevant statistics including the numbers and/or percentages of perfect matches, non perfect matches, key not found, processed reads and undetected guides will be written to statistics_*.csv. Expected running time for each .fastq file generally depends on the size of .fastq files. It typically takes approx. 15 minutes to process a .fastq file containing 20,000,000 reads. 

This Demo_data.fastq.gz is provided just for demonstration to generate truncated results and will not be used for downstream analysis. The full output files required for z-RE calculation are also given in the Demo files folder.

### Calculation of z-RE
You can run calculate_zRE.py just like count_barcodes.py by filling line 146 with your directory where count files are generated.

After running, the Release Effect (RE) at gRNA level and gene level will be written separately to .xlsx files. The Expected running time would be less than a minute. 

## Instructions for use
The only Necessary file to run your data other than .fastq.gz is reference file for barcode counting. Prepared a .xlsx file containing the barcode id in each line to be referenced with a header named ÅeidÅf(refer to ref.xlsx provided). The format of the barcode id should be Åe{target gene}_{spacer sequence}_{additional information}Åf (e.g., ÅeAADACL2_GAAAGTCAGAAACCCGA_2832.7_DTKPÅf) and save this as Åeref.xlsxÅf or change the line 85 of count_barcodes.py:
```sh
ref = 'your_ref.xlsx'
```
 Each read is assigned to a corresponding barcode if the 8 bp of scaffold sequence flanked to the spacer (GTTTAAGA; substituted to KEY in the script, please change if necedssary) is detected and 17 bp of the sequence upstream of KEY is identical to the spacer sequence.

The name of .fastq.gz files need to contain the information of library (e.g., ÅeDTKPÅf), dCas9 anchor (e.g., ÅeCD63Åf), sample origin (ÅecellÅf or ÅesEVsÅf), Cas9 expression (ÅeCas9Åf or ÅenonÅf) and replicate number (e.g., Åerep1Åf) because these information are referenced in calculate_zRE.py
