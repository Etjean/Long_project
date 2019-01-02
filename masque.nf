#!/usr/bin/env nextflow
println ''

// Docker commands (not part of the pipeline, just here to remember them)
// docker run -ti -v /home/etienne/masque/masque.nf:/masque.nf -v /home/etienne/masque/databases:/databases -v /home/etienne/masque/test/data:/mydata shaman_nextflow
// docker run -v /home/etienne/masque/masque.nf:/masque.nf -v /home/etienne/masque/databases:/databases -v /home/etienne/masque/test/data:/mydata -w / --rm shaman_nextflow nextflow masque.nf
// nextflow masque.nf -with-docker shaman_nextflow



// COLORS FOR STRINGS
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";
def tored = {  str -> ANSI_RED + str + ANSI_RESET }
def toblack = {  str -> ANSI_BLACK + str + ANSI_RESET }
def togreen = {  str -> ANSI_GREEN + str + ANSI_RESET }
def toyellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def toblue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def tocyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def topurple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def towhite = {  str -> ANSI_WHITE + str + ANSI_RESET }
def say = {str -> println togreen("* $str")}





// HELP DISPLAY
params.help=false
if (params.size()==1 || params.help) {
	println("Usage:")
	println("16S/18S:   nextflow masque.nf      --i </path/to/input/> --o </path/to/result/>")
	println("23S/28S:   nextflow masque.nf  --l --i </path/to/input/> --o </path/to/result/>")
	println("ITS:       nextflow masque.nf  --f --i </path/to/input/> --o </path/to/result/>")
	println("Amplicon:  nextflow masque.nf      --a <amplicon file>   --o </path/to/result/>")
	println("")
	println("All parameters:")
	println("--i                       Provide </path/to/input/directory/>")
	println("--a                       Provide <amplicon file>")
	println("--o                       Provide </path/to/result/directory/>")
	println("--n                       Indicate <project-name>")
	println("                          (default: use the name of the input directory or meta)")
	println("--t                       Number of <thread>")
	println("                          (default: Nextflow automatic parallelization)")
	println("--c                       Contaminant filtering [danio,human,mouse,mosquito,phi]")
	println("                          (default: human,phi)")
	println("--s                       Perform OTU clustering with swarm")
	println("                          (default: vsearch)")
	println("--b                       Perform taxonomical annotation with blast")
	println("                          (default: vsearch)")
	println("--l                       Perform taxonomical annotation")
	println("                          against LSU databases: Silva/RDP")
	println("--f                       Perform taxonomical annotation")
	println("                          against ITS databases: Unite/Findley/Underhill/RDP")
	println("--minreadlength           Minimum read length take in accound in the study")
	println("                          (default: 35nt)")
	println("--minphred                Qvalue must lie between [0-40]")
	println("                          (default: minimum qvalue 20)")
	println("--minphredperc            Minimum allowed percentage of correctly called")
	println("                          nucleotides [0-100] (default: 80)")
	println("--nbMismatchMapping       Maximum number of mismatch when mapping end-to-end")
	println("                          against Human genome and Phi174 genome")
	println("                          (default: 1 mismatch is accepted)")
	println("--paired                  Paired-ends reads mode")
	println("--minoverlap              Minimum overlap when paired reads are considered")
	println("                          (default: 10)")
	println("--maxoverlap              Maximum overlap when paired reads are considered")
	println("                          (default: 200)")
	println("--minampliconlength       Minimum amplicon length (default: 64)")
	println("--minotusize              Indicate minimum OTU size (default: 4)")
	println("--prefixdrep              Perform prefix dereplication")
	println("                          (default: full length dereplication)")
	println("--chimeraslayerfiltering  Use ChimeraSlayer database for chimera filtering")
	println("                          (default: Perform a de novo chimera filtering)")
	println("--otudiffswarm            Number of difference accepted in an OTU with swarm")
	println("                          (default: 1)")
	println("--evalueTaxAnnot          Evalue threshold for taxonomical annotation with blast")
	println("                          (default: evalue=1E-5)")
	println("--maxTargetSeqs           Number of hit per OTU with blast (default: 1)")
	println("--identityThreshold       Identity threshold for taxonomical annotation with")
	println("                          vsearch (default: 0.75)")
	println("--conservedPosition       Percentage of conserved position in the multiple")
	println("                          alignment considered for phylogenetic tree")
	println("                          (default: 0.8)")
	println("--accurateTree            Accurate tree calculation with IQ-TREE instead of")
	println("                          FastTree (default: FastTree)")
	println("--help                    Print this help")
    exit 1
}



// ASSEMBLY PARAMETERS
params.i=""
params.a=""
params.o=""
readsDir=params.o+'/reads'
logDir=params.o+'/log'
errorlogDir=params.o+'/error_log'
// params.n : see below
params.t=1
params.c=["human", "phi"]
params.s=false
params.b=false
params.l=false
params.f=false
params.minreadlength=35
params.minphred=20
params.minphredperc=80
params.nbMismatchMapping=1
params.paired=false
params.minoverlap=10
params.maxoverlap=200
params.minampliconlength=64
params.minotusize=4
params.prefixdrep=false
params.chimeraslayerfiltering=false
params.otudiffswarm=1
params.evalueTaxAnnot="1E-5"
params.maxTargetSeqs=1
params.identityThreshold=0.75
params.conservedPosition=0.5
params.accurateTree=false



// ARGUMENTS CHECK
// Mandatory parameters
if (params.i=="" && params.a=="") {println tored('Please indicate an input directory or an amplicon file.'); exit 1}
if (params.i!="" && params.a!="") {println tored('Please indicate an input directory or an amplicon file, but not both.'); exit 1}
if (params.o=="") {println tored('Please indicate the output directory.'); exit 1}
// Check and create output directories
outdir = file(params.o)
if( !outdir.exists() ) {
	if( !outdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $outDir"
    } 
}
readsdir = file(readsDir)
if( !readsdir.exists() ) {
    if( !readsdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $readsDir"
    } 
} 
logdir = file(logDir)
if( !logdir.exists() ) {
    if( !logdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $logDir"
    } 
} 
errorlogdir = file(errorlogDir)
if( !errorlogdir.exists() ) {
    if( !errorlogdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $errorlogDir"
    } 
}
// Definition of params.n
if (params.i =~ /\/*([a-zA-Z-_\s]+)\/*$/) {params.n=(params.i =~ /\/*([a-zA-Z-_\s]+)\/*$/)[0][1]}
else if (params.a =~ /\/*([a-zA-Z-._]+\.[a-zA-Z~]+)$/) {params.n=(params.a =~ /\/*([a-zA-Z-._]+)\.[a-zA-Z~]+$/)[0][1]}
else {params.n='masque_run'}






// DATABASES
// Databases directory
db = '/databases'
// ChimeraSlayer reference database
// http://drive5.com/uchime/uchime_download.html
gold="$db/gold.fa"
// Alien sequences
alienseq="$db/alienTrimmerPF8contaminants.fasta"
// Filtering database
filterRef=["danio":"$db/danio_rerio.fna", "human":"$db/homo_sapiens.fna", "mosquito":"$db/anopheles_stephensi.fna", "mouse":"$db/mus_musculus.fna", "phi":"$db/NC_001422.fna"]
// Findley
// http://www.mothur.org/w/images/2/20/Findley_ITS_database.zip
findley="$db/ITSdb.findley.fasta"
// Greengenes
// ftp://greengenes.microbio.me/greengenes_release/gg_13_5/
greengenes="$db/gg_13_5.fasta"
greengenes_taxonomy="$db/gg_13_5_taxonomy.txt"
// Silva
// http://www.arb-silva.de/no_cache/download/archive/release_123/Exports/
silva="$db/SILVA_128_SSURef_Nr99_tax_silva.fasta"
silvalsu="$db/SILVA_128_LSURef_tax_silva.fasta"
underhill="$db/THFv1.3.sequence.fasta"
underhill_taxonomy="$db/THFv1.3.tsv"
unite="$db/sh_general_release_dynamic_s_20.11.2016.fasta"



// PROGRAMS
// Path to the programs
bin='/usr/local/bin'
// AlienTrimmer
alientrimmer="java -jar $bin/AlienTrimmer_0.4.0/src/AlienTrimmer.jar"
// Biom
biom="biom"
// Blastn
blastn="$bin/ncbi-blast-2.5.0+/bin/blastn"
// BMGE ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/
BMGE="java -jar $bin/BMGE-1.12/BMGE.jar"
// Bowtie2
bowtie2="$bin/bowtie2-2.2.9/bowtie2"
// Extract fasta
extract_fasta="$bin/extract_fasta/extract_fasta.py"
// Extract result
extract_result="$bin/extract_result/extract_result.py"
// Fastq2fasta
fastq2fasta="$bin/fastq2fasta/fastq2fasta.py"
// Fastqc
fastqc="$bin/FastQC/fastqc"
// Fasttree
FastTreeMP="$bin/FastTree-2.1.9/FastTree"
// FLASH
flash="$bin/FLASH-1.2.11/flash"
// mafft
mafft="$bin/mafft-linux64/mafft.bat"
// get_taxonomy
get_taxonomy="$bin/get_taxonomy/get_taxonomy.py"
// IQ-TREE
iqtree="$bin/iqtree-omp-1.5.1-Linux/bin/iqtree-omp"
// rename_otu
rename_otu="$bin/rename_otu/rename_otu.py"
// rdp classifier
rdp_classifier="java -jar $bin/rdp_classifier_2.12/dist/classifier.jar"
// swarm
swarm="$bin/swarm_bin/bin/swarm"
// swarm2vsearch
swarm2vsearch="$bin/swarm2vsearch/swarm2vsearch.py"
// vsearch
vsearch="$bin/vsearch_bin/bin/vsearch"



//DISPLAY PARAMETERS
println toblue('Project name [--n]:')
println params.n
if (params.i != '') {
	println toblue('Sample input [--i]:')
	println params.i
}
else if (params.a != '') {
	println toblue('Amplicon input [--a]:')
	println params.a
}
println toblue('Result output [--o]:')
println params.o
println toblue('Reads filtering:')
println "Minimum read length [--minreadlength]= $params.minreadlength"
println "Minimum phred quality [--minphred]= $params.minphred"
println "Minimum allowed percentage of correctly called nucleotides [--minphredperc]= $params.minphredperc"
println "Minimum number of mistach for the filtering [--NbMismatchMapping]= $params.nbMismatchMapping"
print "Filtering databases= "; params.c.each({print "$it; "}); print "\n"
if (params.paired) {
	println toblue('Merge reads parameters')
	println "Min overlap [--minoverlap] = $params.minoverlap"
	println "Max overlap [--maxoverlap] = $params.maxoverlap"
}
println toblue('OTU process:')
if (params.prefixdrep) {
	println 'Dereplication is in mode prefix [--prefixdrep]'
}
else {
	println 'Dereplication is in full length mode'
}
println "Minimum length of an amplicon [--minampliconlength]= $params.minampliconlength"
println "Minimum size of an OTU for singleton removal [--minotusize]= $params.minotusize"
if (params.chimeraslayerfiltering) {
	println "Chimera filtering use chimera slayer database for filtering [--chimeraslayerfiltering]"
}
else {
	println "Chimera filtering is in de novo mode"
}
if (params.s) {
	println "Clustering is performed with swarm [--s]"
	println "Number of difference accepted in an OTU with swarm [--otudiffswarm]= $params.otudiffswarm"
}
else {
	println "Clustering is performed with vsearch"
}
if (params.f) {
	println toblue("Fungi annotation [--f]")
}
else if (params.l) {
	println toblue("23S/28S annotation [--l]")
}
else {
	println toblue("16S/18S annotation")
}
if(params.b) {
	println "E-value with blast [--evalueTaxAnnot]= $params.evalueTaxAnnot"
	println "Maximum number of targets with blast [--maxTargetSeqs]= $params.maxTargetSeqs"
}
else {
	println "Identity threshold with vsearch [--identityThreshold]= $params.identityThreshold"
}
println "Conserved position for alignment[--conservedPosition]= $params.conservedPosition"
if (params.accurateTree) {
	println "Tree generated in accurate mode with IQ-TREE [--accurateTree]"
}
else {
	println "Tree generated in fast mode with FastTree"
}




// MAIN
say 'Start analysis'

say 'Start working on reads'

// Input channel
reads_raw_compressed = Channel.fromPath("$params.i/*.{fq.gz,fastq.gz}")
reads_raw = Channel.fromPath("$params.i/*.{fq,fastq}")



// Decompress
// Only executed if input files have the '.gz' extension.
process Decompress {
	tag "$sample"
	
	input:
		file reads_raw_compressed
		
	output:
		file "${sample}.fastq" into reads_dc
	
	script:
		sample = reads_raw_compressed.simpleName
		"""
		gzip --decompress --stdout $reads_raw_compressed > ${sample}.fastq
		"""
}



// Trimming
process Trimming {
	tag "$sample"
	publishDir logDir, mode: 'move', pattern: 'log_alientrimmer_*'
	publishDir errorlogDir, mode: 'move', pattern: 'error_log_alientrimmer_*'
	
	input:
		file reads_raw from reads_raw.mix(reads_dc)
	
	output:
		file "${sample}_alien.fastq" into reads_alien
		file "log_alientrimmer_${sample}.txt"
		file "error_log_alientrimmer_${sample}.txt"
		
	
	script:
		sample = reads_raw.simpleName
		"""
		$alientrimmer \
		-i $reads_raw \
		-o ${sample}_alien.fastq \
		-c $alienseq \
		-l $params.minreadlength \
		-p $params.minphredperc \
		-q $params.minphred \
		> log_alientrimmer_${sample}.txt \
		2> error_log_alientrimmer_${sample}.txt
		"""
}


contminants_species = Channel.from(params.c.collect( {filterRef[it]} ))

// // Filetring reads against contaminant database
// process Comtaminant {
	
// 	input:
// 		file reads_alien
// 		val contminants_species
	
// 	output:
// 		file '*.fastq'

// 	script:
// 		sample = (reads_alien.name =~ /(.+)_alien\.fastq/)[0][1]
// 		cont = 
// 		"""
// 		$bowtie2  -q \
// 		-N $params.nbMismatchMapping \
// 		-p $params.t \
// 		-x $contaminant_species \
// 		-U $reads_alien \
// 		-S /dev/null \
// 		--un ${sample}_filtered_${cont}.fastq \
// 		-t --end-to-end --very-fast  \
// 		> ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt \
// 		2>&1
// 		"""
// }

















