# AOPWikiXMLparsing
[![DOI](https://zenodo.org/badge/139583189.svg)](https://zenodo.org/badge/latestdoi/139583189)

This repository contains the python script to extract data from the AOP-Wiki XML file, to support the paper "Introducing WikiPathways as a data-source to support Adverse Outcome Pathways for regulatory risk assessment of chemicals and nanomaterials.".

## Step 1: Download the input files

###a: 
The AOP-Wiki provides quarterly downloads. The full AOP-Wiki XML files should be downloaded from: https://aopwiki.org/downloads. To reproduce the results of the paper, download the AOP Wiki XML of April 1st, 2018.
### b: 
The HGNC genelist was retrieved from HGNC custom download page (https://www.genenames.org/cgi-bin/download). However, this input file is also present in this repository.

## Step 2: Prepare the script

The script requires setting file paths. Make sure that all file paths are set correctly to the AOP Wiki XML, the HGNC gene name file, and the location for the output file (file name: writtenoutput.txt) as well as the SPARQL queries and the file containing all CAS Registry Numbers (CAS-RNs, file name: ListofCasrns.txt).

## Step 3: Run the script

This script was used with Python 3.5.

## Step 4: Using the output files

### a: 
The general output file (writtenoutput.txt) summarises all relevant information based on the data in the AOP-Wiki XML, readable for the human. 
### b: 
The SPARQL queries, that are stored as text files, should be entered in the WikiPathways SPARQL endpoint. In return, an output will be given as an answer to the SPARQL query.
### c: 
The ListofCasrns.txt file is used to map the CAS-RNs from the AOP-Wiki XML to ChEBI IDs, because the SPARQL endpoint does not allow queries for CAS-RNs directly.

## Step 5: Identifier mapping in R


