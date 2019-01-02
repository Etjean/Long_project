# Long Project
Nextflow pipeline and Docker integration into SHAMAN.

## Background
Quantitative metagenomics is broadly employed to identify genera or species associated with several diseases. These data are obtained by mapping the reads of each sample against operational taxonomic units (OTU) or a gene catalog. SHAMAN was one the first web application that allowed to clinician and biologist to perform an interactive analysis of quantitative metagenomics data with a dynamic-interface dedicated to the diagnostic and to the differential analysis. The interface integrates the experimental design (association of sample to one or several conditions), the statistical process for differential analysis and a real-time visualisation system.  

SHAMAN is based on R, Shiny and DESeq2. The analytical process is divided into four steps : count matrix/annotation submission, normalisation, modelisation and visualisation. The count matrix is normalised at the OTU/gene level using the DESeq2 normalisation method and then, based on the experimental design, a generalised linear model is applied to detect differences in abundance at the considered taxonomic level.  

Two years after first release, we can see a great interest from the metagenomics community with 5 publications using SHAMAN (in Nature microbiology, PNAS and Science Advances), 3 publications not involving SHAMAN authors, 74 active users per month (1430 unique visitors since first publication - 70 % are regular users) and 514 downloads of the Docker application. Several trainings were also performed to train biologist to use SHAMAN in the Pasteur Network and at ENS.  

## Project
We want to integrate a full automatized bioinformatic workflow based on Nextflow for targeted metagenomics data. This implementation will follow the current approach already implemented in bash workflow (MASQUE pipeline). The workflow should also be included in the SHAMAN Docker application for local installation on windows/mac/linux.  

## Usage
- **Download databases** 
  ```
  bash install_databases.sh
  ```

- **Run the pipeline**
  Run the script directly :
  ```
  docker run --rm \
  	-v /path/to/databases:/databases \
  	-v /path/to/data:/mydata \
  	etjean/shaman_nextflow \
  	nextflow masque.nf [OPTIONS]
  ```
  Or open an interactive container first and then run the script :
  ```
  docker run -ti \
  	-v /path/to/databases:/databases \
  	-v /path/to/data:/mydata \
  	etjean/shaman_nextflow
  
  nextflow masque.nf [OPTIONS]
  ```


## Resources
- **MASQUE pipeline** : Metagenomic Analysis with a Quantitative pipeline - <https://github.com/aghozlane/masque>  
- **SHAMAN application** : Shiny Application for Metagenomic Analysis - <http://shaman.pasteur.fr/>  
- Contact : Amine Ghozlane - amine.ghozlane@pasteur.fr  
