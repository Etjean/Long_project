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
  Download databases [here](http://dl.pasteur.fr/fop/vJlf2Krl/database.zip)

- **Run the pipeline**  
  Run the script directly :
  ```
  docker run --rm \
  	-v /path/to/databases:/databases \
  	-v /path/to/data:/mydata \
  	etjean/shaman_nextflow \
  	nextflow masque.nf --i /mydata --o /mydata/result [OPTIONS]
  ```
  Or open an interactive container first, and then run the script :
  ```
  docker run -ti \
  	-v /path/to/databases:/databases \
  	-v /path/to/data:/mydata \
  	etjean/shaman_nextflow
  
  nextflow masque.nf --i /mydata --o /mydata/result [OPTIONS]
  ```
  
- **Arguments**
  ```
  Usage:
  16S/18S:   nextflow masque.nf      --i </path/to/input/> --o </path/to/result/>
  23S/28S:   nextflow masque.nf  --l --i </path/to/input/> --o </path/to/result/>
  ITS:       nextflow masque.nf  --f --i </path/to/input/> --o </path/to/result/>
  Amplicon:  nextflow masque.nf      --a <amplicon file>   --o </path/to/result/>

  All parameters:
  --i                       Provide </path/to/input/directory/>
  --a                       Provide <amplicon file>
  --o                       Provide </path/to/result/directory/>
  --n                       Indicate <project-name>
                            (default: use the name of the input directory or meta)
  --t                       Number of <thread>
                            (default: Nextflow automatic parallelization)
  --c                       Contaminant filtering [danio,human,mouse,mosquito,phi]
                            (default: human,phi)
  --s                       Perform OTU clustering with swarm
                            (default: vsearch)
  --b                       Perform taxonomical annotation with blast
                            (default: vsearch)
  --l                       Perform taxonomical annotation
                            against LSU databases: Silva/RDP
  --f                       Perform taxonomical annotation
                            against ITS databases: Unite/Findley/Underhill/RDP
  --minreadlength           Minimum read length take in accound in the study
                            (default: 35nt)
  --minphred                Qvalue must lie between [0-40]
                            (default: minimum qvalue 20)
  --minphredperc            Minimum allowed percentage of correctly called
                            nucleotides [0-100] (default: 80)
  --nbMismatchMapping       Maximum number of mismatch when mapping end-to-end
                            against Human genome and Phi174 genome
                            (default: 1 mismatch is accepted)
  --paired                  Paired-ends reads mode
  --minoverlap              Minimum overlap when paired reads are considered
                            (default: 10)
  --maxoverlap              Maximum overlap when paired reads are considered
                            (default: 200)
  --minampliconlength       Minimum amplicon length (default: 64)
  --minotusize              Indicate minimum OTU size (default: 4)
  --prefixdrep              Perform prefix dereplication
                            (default: full length dereplication)
  --chimeraslayerfiltering  Use ChimeraSlayer database for chimera filtering
                            (default: Perform a de novo chimera filtering)
  --otudiffswarm            Number of difference accepted in an OTU with swarm
                            (default: 1)
  --evalueTaxAnnot          Evalue threshold for taxonomical annotation with blast
                            (default: evalue=1E-5)
  --maxTargetSeqs           Number of hit per OTU with blast (default: 1)
  --identityThreshold       Identity threshold for taxonomical annotation with
                            vsearch (default: 0.75)
  --conservedPosition       Percentage of conserved position in the multiple
                            alignment considered for phylogenetic tree
                            (default: 0.8)
  --accurateTree            Accurate tree calculation with IQ-TREE instead of
                            FastTree (default: FastTree)
  --help                    Print this help
  ```

## Problems
- [ ] `--t` argument is obsolete.
- [ ] `--c` argument is currently non fonctionnal.
- [ ] More arguments control is needed.


## Resources
- **MASQUE pipeline** : Metagenomic Analysis with a Quantitative pipeline - <https://github.com/aghozlane/masque>  
- **SHAMAN application** : Shiny Application for Metagenomic Analysis - <http://shaman.pasteur.fr/>  
- Contact : Amine Ghozlane - amine.ghozlane@pasteur.fr  
