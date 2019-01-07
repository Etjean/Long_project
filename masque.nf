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
	// println("--c                       Contaminant filtering [danio,human,mouse,mosquito,phi]")
	// println("                          (default: human,phi)")
	println("--s                       Perform OTU clustering with swarm")
	println("                          (default: vsearch)")
	println("--b                       Perform taxonomical annotation with blast")
	println("                          (default: vsearch)")
	println("--l                       Perform taxonomical annotation against LSU")
	println("                          databases: Silva/RDP")
	println("--f                       Perform taxonomical annotation against ITS")
	println("                          databases: Unite/Findley/Underhill/RDP")
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
if (params.i =~ /\/*([^\/]+)\/*$/) {params.n=(params.i =~ /\/*([^\/]+)\/*$/)[0][1]}
else if (params.a =~ /\/*([^\/]+)\/*$/) {params.n=(params.a =~ /\/*([^\/]+)\/*$/)[0][1]}
else {params.n='masque_run'}




// DATABASES
// Databases directory
db = '/db'
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


// Input channel
reads_raw=Channel.empty()
reads_raw_compressed=Channel.empty()
reads_raw_paired=Channel.empty()
reads_raw_compressed_paired=Channel.empty()
if (params.i!="" && !params.paired && params.a=="") {
	reads_raw_compressed = Channel.fromPath("$params.i/*.{fq.gz,fastq.gz}")
	reads_raw = Channel.fromPath("$params.i/*.{fq,fastq}")
}
else if (params.i!="" && params.paired && params.a=="") {
	reads_raw_compressed_paired = Channel.fromFilePairs("$params.i/*{R1,R2}*.{fq.gz,fastq.gz}")
	reads_raw_paired = Channel.fromFilePairs("$params.i/*{R1,R2}*.{fq,fastq}")
	// reads_raw_compressed = Channel.fromPath("$params.i/*.{fq.gz,fastq.gz}")
	// reads_raw = Channel.fromPath("$params.i/*.{fq,fastq}")
}
else if (params.i=="" && !params.paired && params.a!="") {
	if (params.a =~ /(fq|fastq).gz$/) {
		reads_raw_compressed = Channel.fromPath("$params.a")
		reads_raw = Channel.empty()
	}
	else if (params.a =~ /(fq|fastq)$/) {
		reads_raw_compressed = Channel.empty()
		reads_raw = Channel.fromPath("$params.a")
	}
}

// Decompress
// Only executed if input files have the '.gz' extension.
// say 'Decompressing reads files'

process Decompress {
	tag "$sample"
	
	input:
		file reads_raw_compressed
		
	output:
		file "${sample}.fastq" into reads_dc
					
	when:
		!params.paired

	script:
		sample = reads_raw_compressed.simpleName
		"""
		gzip --decompress --stdout $reads_raw_compressed > ${sample}.fastq
		"""
}



// Trimming
// say 'Trimming reads with AlienTrimmer'

process Trimming {
	tag "$sample"
	publishDir logDir, mode: 'copy', pattern: 'log_alientrimmer_*'
	publishDir errorlogDir, mode: 'copy', pattern: 'error_log_alientrimmer_*'
	publishDir readsDir, mode:'copy', pattern: '*_alien.fastq'
	
	input:
		file reads_raw from reads_raw.mix(reads_dc)
	
	output:
		file "${sample}_alien.fastq" into reads_alien
		file "log_alientrimmer_${sample}.txt"
		file "error_log_alientrimmer_${sample}.txt"
		
	when:
		!params.paired
	
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



// Pas robuste
// reads_alien = reads_alien
// 	.groupBy { f -> (f.name=~/(.*)(R1|R2)(.*).fastq$/)[0][1]+(f.name=~/(.*)(R1|R2)(.*).fastq$/)[0][3] }
// 	.subscribe {println it}

process DecompressPaired {
	tag "$sample"
	publishDir readsDir, mode:'copy'
		
	input:
		set sample, file(rds) from reads_raw_compressed_paired
		
	// output:
	// 	set sample [file() "${sample}.fastq" into reads_dc
		
	when:
		params.paired
	
	script:
	println sample
	println rds
		"""
		gzip --decompress --stdout ${rds[0]} > ${rds[0].name.take(rds[0].name.lastIndexOf('.gz'))}
		gzip --decompress --stdout ${rds[1]} > ${rds[1].name.take(rds[1].name.lastIndexOf('.gz'))}
		"""
}



// Trimming
// say 'Trimming reads with AlienTrimmer'

// process TrimmingPaired {
// 	tag "$sample"
// 	publishDir logDir, mode: 'copy', pattern: 'log_alientrimmer_*'
// 	publishDir errorlogDir, mode: 'copy', pattern: 'error_log_alientrimmer_*'
// 	publishDir readsDir, mode:'copy', pattern: '*_alien.fastq'
	
// 	input:
// 		set key file(rds) from reads_raw.mix(reads_dc)
	
// 	output:
// 		file "${sample}_alien.fastq" into reads_alien
// 		file "log_alientrimmer_${sample}.txt"
// 		file "error_log_alientrimmer_${sample}.txt"
	
// 	when:
// 		params.paired
	
// 	script:
// 		"""
// 		$alientrimmer \
// 		-if $input1 \
// 		-ir $input2 \
// 		-of ${readsDir}/${SampleName}_alien_f.fastq \
// 		-or ${readsDir}/${SampleName}_alien_r.fastq \
// 		-os ${readsDir}/${SampleName}_alien_s.fastq \
// 		-c $alienseq \
// 		-l $minreadlength \
// 		-p $minphredperc \
// 		-q $minphred \
// 		> $log_alientrimmer_${sample}.txt  \
// 		2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
// 		"""
// }












// process Merging {
// 	tag "$sample"
// 	publishDir logDir, mode: 'copy', pattern: 'log_flash_*'
// 	// publishDir readsDir, mode:'copy', pattern: ''
	
// 	input:
// 	set sample, file(rds) from reads_alien
// 	// output:
	
	
// 	when:
// 	params.paired
	
// 	script:
// 	// """
// 	// $flash \
// 	// ${rds[0]} \
// 	// ${rds[1]} \
// 	// -M $params.maxoverlap \
// 	// -m $params.minoverlap \
// 	// -d $readsDir/ \
// 	// -o out \
// 	// -t 1  \
// 	// > log_flash_truc.txt
// 	// """
// 	"""$vsearch  \
// 	--fastq_mergepairs ${rds[0]} \
// 	--reverse ${rds[1]} \
// 	--fastqout out.fastq \
// 	--fastq_minovlen $params.minoverlap \
// 	--threads 1
// 	"""
// }








// // // Filetring reads against contaminant database
// // contminants_species = Channel.from(params.c.collect( {filterRef[it]} ))
// // process Comtaminant {
// // 	input:
// // 		file reads_alien
// // 		val contminants_species
// // 	output:
// // 		file '*.fastq'
// // 	script:
// // 		sample = (reads_alien.name =~ /(.+)_alien\.fastq/)[0][1]
// // 		cont = 
// // 		"""
// // 		$bowtie2  -q \
// // 		-N $params.nbMismatchMapping \
// // 		-p 1 \
// // 		-x $contaminant_species \
// // 		-U $reads_alien \
// // 		-S /dev/null \
// // 		--un ${sample}_filtered_${cont}.fastq \
// // 		-t --end-to-end --very-fast  \
// // 		> ${logDir}/log_mapping_${SampleName}_${contaminant[${essai}]}_${essai}.txt \
// // 		2>&1
// // 		"""
// // }





// // Convert to fasta
// // say 'Convert files from fastq to fasta'

// process Fastq2Fasta {
// 	tag "$sample"
// 	publishDir readsDir, mode:'copy', pattern:'*.fasta'
// 	publishDir errorlogDir, mode:'copy', pattern:'error_log_fastq2fasta_*'
	
// 	input:
// 	file reads_alien
	
// 	output:
// 	file "${sample}.fasta" into reads_fasta
// 	file "error_log_fastq2fasta_${sample}.txt"
	
// 	script:
// 	sample = (reads_alien.name =~ /(.+)_alien\.fastq/)[0][1]
// 	"""
// 	$fastq2fasta \
// 	-i $reads_alien \
// 	-o ${sample}.fasta \
// 	-s ${sample}  \
// 	2> error_log_fastq2fasta_${sample}.txt
// 	"""
// }



// // Concatenate files
// // say 'Combine all fasta files'

// reads_fasta
// 	.collectFile(name: params.n+'_trimmed.fasta', storeDir: readsDir)
// 	.into {reads_fasta1; reads_fasta2}



// // Dereplication
// // if (params.prefixdrep) say 'Dereplication using prefixes'
// // else say 'Dereplication using full reads lentgh'

// process Dereplication {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy', pattern:'*_drep.fasta'
		
// 	input:
// 	file reads_fasta1
	
// 	output:
// 	file "*_drep.fasta" into reads_drep
// 	// file "log_vsearch_dereplication_*"
	
// 	script:
// 	if (params.prefixdrep)
// 		"""
// 		$vsearch \
// 		--derep_prefix $reads_fasta1 \
// 		-output ${params.n}_drep.fasta \
// 		-sizeout \
// 		-minseqlength $params.minampliconlength \
// 		"""
// 	else
// 		"""
// 		$vsearch \
// 		--derep_fulllength $reads_fasta1 \
// 		-output ${params.n}_drep.fasta \
// 		-sizeout \
// 		-minseqlength $params.minampliconlength \
// 		--strand both
// 		"""
// }



// // Singleton removal
// // say 'Abundance sorting and Singleton removal'

// process SingletonRemoval {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy', pattern:'*_sorted.fasta'
// 	publishDir logDir, mode:'copy', pattern:'log_search_sort_*'
	
// 	input:
// 	file reads_drep
	
// 	output:
// 	file "*_sorted.fasta" into reads_nosing
	
// 	script:
// 	"""
// 	$vsearch \
// 	-sortbysize $reads_drep \
// 	-output ${params.n}_sorted.fasta \
// 	-minsize $params.minotusize \
//  	> log_search_sort_${params.n}.txt \
//  	2>&1
//  	"""
// }



// // Chimeras filtering
// // if (params.chimeraslayerfiltering) say 'Chimeric reads filtering with ChimeraSlayer database'
// // else say 'Chimeric reads de novo filtering'

// process ChimerasRemoval {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy', pattern:'*_nochim.fasta'
// 	publishDir readsDir, mode:'copy', pattern:'*_chim.fasta'
	
// 	input:
// 	file reads_nosing
	
// 	output:
// 	file '*_nochim.fasta' into reads_nochim
// 	file '*_chim.fasta'
	
// 	script:
// 	if (params.chimeraslayerfiltering)
// 		"""
// 		$vsearch \
// 		--uchime_ref $reads_nosing \
// 		--db $gold \
// 		--strand both \
// 		--nonchimeras ${params.n}_nochim.fasta \
// 		--chimeras ${params.n}_chim.fasta
// 		"""
// 	else
// 		"""
// 		$vsearch \
// 		--uchime_denovo $reads_nosing \
// 		--strand both \
// 		--nonchimeras ${params.n}_nochim.fasta \
// 		--chimeras ${params.n}_chim.fasta
// 		"""
// }



// // Clustering
// // if (params.s) say 'Clustering reads using Swarm'
// // else say 'Clustering reads using Vsearch'

// process Clustering {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy', pattern:'*_otu.fasta'
	
// 	input:
// 	file reads_nochim
	
// 	output:
// 	file '*_otu.fasta' into otu
	
// 	script:
// 	if (params.s)
// 		"""
// 		$swarm \
// 		-t 1 \
// 		-f \
// 		-z \
// 		-w ${params.n}_otu_tmp.fasta \
// 		-o ${params.n}_swarm_clustering.txt \
// 		-s ${params.n}_swarm_stats.txt \
// 		-u ${params.n}_swarm_uclust.txt \
// 		$reads_nochim
		
// 		python $swarm2vsearch \
// 		-i ${params.n}_otu_tmp.fasta \
// 		-c ${params.n}_swarm_clustering.txt \
// 		-o ${params.n}_otu.fasta \
// 		-oc ${params.n}_otu_swarm_clustering.txt \
// 		-u ${params.n}_swarm_uclust.txt \
// 		-ou ${params.n}_otu_swarm_uclust.txt
// 		"""
// 	else
// 		"""
// 		$vsearch \
// 		--cluster_size $reads_nochim \
// 		--id 0.97 \
// 		--centroids ${params.n}_otu.fasta \
// 		--sizein \
// 		--strand both
// 		"""
// }
// otu.into {otu1; otu2; otu3; otu4; otu5; otu6}



// // Mapping back reads to OTUs
// // say 'Mapping back reads to OTUs'

// process Mapping {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy'
	
// 	input:
// 	file otu1
// 	file reads_fasta2
	
// 	output:
// 	file "${params.n}_otu_table.tsv"
// 	file "${params.n}_count.biom" into count
	
// 	script:
// 	"""
// 	$vsearch \
// 	-usearch_global $reads_fasta2\
// 	-db $otu1 \
// 	--strand both \
// 	--id 0.97 \
// 	--otutabout ${params.n}_otu_table.tsv \
// 	--biomout ${params.n}_count.biom
// 	"""
// }
// count.into {count1; count2; count3}



// // Taxonomy Annotation
// // if (params.b) say 'Taxonomy annotation using Blast'
// // else say 'Taxonomy annotation using Vsearch'
// // if (!params.l && !params.f) say 'Databases : RDP, Greengenes, Silva'
// // else if (params.l && !params.f) say 'Databases : RDP, Silva'
// // else if (params.f) say 'Databases : RDP, Findley, Unite, Underhill'

// process TaxonomyAnnotationRDP {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy'
	
// 	input:
// 	file otu2
	
// 	output:
// 	file "${params.n}_vs_rdp_annotation.tsv" into annotation_rdp
	
// 	script:
// 	// RDP
// 	"""
// 	$rdp_classifier classify \
// 	-q $otu2 \
// 	-o ${params.n}_vs_rdp_annotation.tsv
// 	"""
// }

// process TaxonomyAnnotationGreengenes {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy'
	
// 	input:
// 	file otu3
// 	file count1
	
// 	output:
// 	file "${params.n}_vs_greengenes_annotation_*.tsv" into annotation_greengenes
// 	file "${params.n}_greengenes_*.biom"
	
// 	when:
// 	!params.l && !params.f
	
// 	script:
// 	// Grenngenes for SSU with Vsearch
// 	if (!params.b) {
// 		"""
// 	    $vsearch \
// 	    --usearch_global $otu3 \
// 	    --db $greengenes \
// 	    --id $params.identityThreshold \
// 	    --blast6out ${params.n}_vs_greengenes_id_${params.identityThreshold}.tsv \
// 	    --strand both
	    
//         python $get_taxonomy \
//         -i ${params.n}_vs_greengenes_id_${params.identityThreshold}.tsv \
//         -d $greengenes \
//         -o ${params.n}_vs_greengenes_annotation_id_${params.identityThreshold}.tsv \
//         -dtype greengenes \
//         -t $greengenes_taxonomy \
//         -ob ${params.n}_vs_greengenes_annotation_id_${params.identityThreshold}.biomtsv \
//         -u $otu3
        
//         $biom add-metadata \
//         -i $count1 \
//         -o ${params.n}_greengenes_id_${params.identityThreshold}.biom \
//         --observation-metadata-fp ${params.n}_vs_greengenes_annotation_id_${params.identityThreshold}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json
// 		"""
// 	}
	
// 	// Greengenes for SSU with Blast
// 	else if (params.b) {
// 		"""
//         $blastn \
//         -query $otu3 \
//         -db $greengenes \
//         -evalue $params.evalueTaxAnnot \
//         -num_threads 1 \
//         -out ${params.n}_vs_greengenes_eval_${params.evalueTaxAnnot}.tsv \
//         -max_target_seqs $params.maxTargetSeqs \
//         -task megablast \
//         -outfmt "6 qseqid sseqid  pident qcovs evalue" \
//         -use_index true
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_greengenes_eval_${params.evalueTaxAnnot}.tsv \
//         -d $greengenes \
//         -u $otu3 \
//         -o ${params.n}_vs_greengenes_annotation_eval_${params.evalueTaxAnnot}.tsv \
//         -dtype greengenes \
//         -t $greengenes_taxonomy \
//         -ob ${params.n}_vs_greengenes_annotation_eval_${params.evalueTaxAnnot}.biomtsv
        
//         $biom add-metadata \
//         -i $count1 \
//         -o ${params.n}_greengenes_eval_${params.evalueTaxAnnot}.biom \
//         --observation-metadata-fp ${params.n}_vs_greengenes_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json
// 		"""
// 	}
// }

// process TaxonomyAnnotationSilva {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy'
	
// 	input:
// 	file otu4
// 	file count2
	
// 	output:
// 	file "${params.n}_vs_silva_annotation_*.tsv" into annotation_silva
// 	file  "${params.n}_silva_*.biom"
	
// 	when:
// 	!params.f
	
// 	script:
// 	// Silva with Vsearch
// 	if (!params.b) {
// 		// SSU
// 		if (!params.l) {
// 			"""
// 	        $vsearch \
// 	        --usearch_global $otu4 \
// 	        --db $silva \
// 	        --id $params.identityThreshold \
// 	        --blast6out ${params.n}_vs_silva_id_${params.identityThreshold}.tsv \
// 	        --strand both
	        
// 	        python $get_taxonomy \
// 	        -i ${params.n}_vs_silva_id_${params.identityThreshold}.tsv \
// 	        -u $otu4 \
// 	        -d $silva \
// 	        -o ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.tsv \
// 	        -ob ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.biomtsv
	        
// 		    $biom add-metadata \
// 		    -i $count2 \
// 		    -o ${params.n}_silva_id_${params.identityThreshold}.biom \
// 		    --observation-metadata-fp ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.biomtsv \
// 		    --observation-header id,taxonomy \
// 		    --sc-separated taxonomy \
// 		    --output-as-json
// 			"""
// 		}
// 		// LSU
// 		else if (params.l) {
// 			"""
// 			$vsearch \
// 			--usearch_global $otu4 \
// 			--db $silvalsu \
// 			--id $params.identityThreshold \
// 			--blast6out ${params.n}_vs_silva_id_${params.identityThreshold}.tsv \
// 			--strand both
			
// 	        python $get_taxonomy \
// 	        -i ${params.n}_vs_silva_id_${params.identityThreshold}.tsv \
// 	        -u $otu4 \
// 	        -d $silvalsu \
// 	        -o ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.tsv \
// 	        -ob ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.biomtsv
	        
// 	        $biom add-metadata \
// 		    -i $count2 \
// 		    -o ${params.n}_silva_id_${params.identityThreshold}.biom \
// 		    --observation-metadata-fp ${params.n}_vs_silva_annotation_id_${params.identityThreshold}.biomtsv \
// 		    --observation-header id,taxonomy \
// 		    --sc-separated taxonomy \
// 		    --output-as-json
// 			"""
// 		}
// 	}
	
// 	// Silva with Blast
// 	else if (params.b) {
// 		// SSU
// 		if (!params.l) {
// 			"""
// 	        $blastn \
// 	        -query $otu4 \
// 	        -db $silva \
// 	        -evalue $params.evalueTaxAnnot \
// 	        -num_threads 1 \
// 	        -out ${params.n}_vs_silva_eval_${params.evalueTaxAnnot}.tsv \
// 	        -max_target_seqs $params.maxTargetSeqs \
// 	        -task megablast \
// 	        -outfmt "6 qseqid sseqid  pident qcovs evalue" \
// 	        -use_index true 
	        
// 	        python $get_taxonomy \
// 	        -i ${params.n}_vs_silva_eval_${params.evalueTaxAnnot}.tsv \
// 	        -d $silva \
// 	        -u $otu4 \
// 	        -o ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.tsv \
// 	        -ob ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.biomtsv
	        
// 		    $biom add-metadata \
// 		    -i $count2 \
// 		    -o ${params.n}_silva_eval_${params.evalueTaxAnnot}.biom \
// 		    --observation-metadata-fp ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
// 		    --observation-header id,taxonomy \
// 		    --sc-separated taxonomy \
// 		    --output-as-json
// 			"""
// 		}
// 		// LSU
// 		else if (params.l) {
// 			"""
// 	        $blastn \
// 	        -query $otu4 \
// 	        -db $silvalsu \
// 	        -evalue $params.evalueTaxAnnot \
// 	        -num_threads 1 \
// 	        -out ${params.n}_vs_silva_eval_${params.evalueTaxAnnot}.tsv \
// 	        -max_target_seqs $params.maxTargetSeqs \
// 	        -task megablast \
// 	        -outfmt "6 qseqid sseqid  pident qcovs evalue" \
// 	        -use_index true 
	        
// 	        python $get_taxonomy \
// 	        -i ${params.n}_vs_silva_eval_${params.evalueTaxAnnot}.tsv \
// 	        -d $silvalsu \
// 	        -u $otu4 \
// 	        -o ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.tsv \
// 	        -ob ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.biomtsv
	        
// 	        $biom add-metadata \
// 		    -i $count2 \
// 		    -o ${params.n}_silva_eval_${params.evalueTaxAnnot}.biom \
// 		    --observation-metadata-fp ${params.n}_vs_silva_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
// 		    --observation-header id,taxonomy \
// 		    --sc-separated taxonomy \
// 		    --output-as-json
// 			"""	
// 		}
// 	}
// }

// process TaxonomyAnnotationITS {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy'
	
// 	input:
// 	file otu5
// 	file count3
	
// 	output:
// 	file "${params.n}_vs_*_annotation_*.tsv" into annotation_ITS
// 	file "${params.n}_*.biom"
	
// 	when:
// 	params.f
	
// 	script:
// 	// Findley/Unite/Underhill with Vsearch
// 	if (!params.b) {
// 		"""
//         $vsearch --usearch_global $otu5 \
//         --db $findley \
//         --id $params.identityThreshold \
//         --blast6out ${params.n}_vs_findley_id_${params.identityThreshold}.tsv \
//         --strand both
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_findley_id_${params.identityThreshold}.tsv \
//         -d $findley  \
//         -u $otu5 \
//         -o ${params.n}_vs_findley_annotation_id_${params.identityThreshold}.tsv \
//         -ob ${params.n}_vs_findley_annotation_id_${params.identityThreshold}.biomtsv  \
//         -dtype findley 
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_findley_id_${params.identityThreshold}.biom \
//         --observation-metadata-fp ${params.n}_vs_findley_annotation_id_${params.identityThreshold}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json 

//         $vsearch --usearch_global $otu5 \
//         --db $unite \
//         --id $params.identityThreshold \
//         --blast6out ${params.n}_vs_unite_id_${params.identityThreshold}.tsv \
//         --strand both
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_unite_id_${params.identityThreshold}.tsv \
//         -d $unite \
//         -u  $otu5  \
//         -o ${params.n}_vs_unite_annotation_id_${params.identityThreshold}.tsv \
//         -ob ${params.n}_vs_unite_annotation_id_${params.identityThreshold}.biomtsv \
//         -dtype unite
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_unite_id_${params.identityThreshold}.biom \
//         --observation-metadata-fp ${params.n}_vs_unite_annotation_id_${params.identityThreshold}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json

//         $vsearch --usearch_global $otu5 \
//         --db $underhill \
//         --id $params.identityThreshold \
//         --blast6out ${params.n}_vs_underhill_id_${params.identityThreshold}.tsv \
//         --strand both
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_underhill_id_${params.identityThreshold}.tsv \
//         -d $underhill \
//         -u  $otu5 \
//         -t $underhill_taxonomy \
//         -o ${params.n}_vs_underhill_annotation_id_${params.identityThreshold}.tsv \
//         -ob ${params.n}_vs_underhill_annotation_id_${params.identityThreshold}.biomtsv \
//         -dtype underhill
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_underhill_id_${params.identityThreshold}.biom \
//         --observation-metadata-fp ${params.n}_vs_underhill_annotation_id_${params.identityThreshold}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json
// 		"""
// 	}
	
// 	// Findley/Unite/Underhill with Blast
// 	else if (params.b) {
// 		"""
//         $blastn \
//         -query $otu5 \
//         -db $findley \
//         -evalue $params.evalueTaxAnnot \
//         -num_threads 1 \
//         -out ${params.n}_vs_findley_eval_${params.evalueTaxAnnot}.tsv \
//         -max_target_seqs $params.maxTargetSeqs \
//         -task megablast \
//         -outfmt "6 qseqid sseqid  pident qcovs evalue" \
//         -use_index true
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_findley_eval_${params.evalueTaxAnnot}.tsv \
//         -d $findley \
//         -u  $otu5  \
//         -o ${params.n}_vs_findley_annotation_eval_${params.evalueTaxAnnot}.tsv \
//         -ob ${params.n}_vs_findley_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         -dtype findley
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_findley_eval_${params.evalueTaxAnnot}.biom \
//         --observation-metadata-fp ${params.n}_vs_findley_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json

// 		$blastn \
// 		-query $otu5 \
// 		-db $unite \
// 		-evalue $params.evalueTaxAnnot \
// 		-num_threads 1 \
// 		-out ${params.n}_vs_unite_eval_${params.evalueTaxAnnot}.tsv \
// 		-max_target_seqs $params.maxTargetSeqs \
// 		-task megablast \
// 		-outfmt "6 qseqid sseqid  pident qcovs evalue" \
// 		-use_index true
		
//         python $get_taxonomy \
//         -i ${params.n}_vs_unite_eval_${params.evalueTaxAnnot}.tsv \
//         -d $unite \
//         -u  $otu5  \
//         -o ${params.n}_vs_unite_annotation_eval_${params.evalueTaxAnnot}.tsv \
//         -ob ${params.n}_vs_unite_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         -dtype unite
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_unite_eval_${params.evalueTaxAnnot}.biom \
//         --observation-metadata-fp ${params.n}_vs_unite_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json

//         $blastn \
//         -query $otu5 \
//         -db $underhill \
//         -evalue $params.evalueTaxAnnot \
//         -num_threads 1 \
//         -out ${params.n}_vs_underhill_eval_${params.evalueTaxAnnot}.tsv \
//         -max_target_seqs $params.maxTargetSeqs \
//         -task megablast \
//         -outfmt "6 qseqid sseqid  pident qcovs evalue" \
//         -use_index true
        
//         python $get_taxonomy \
//         -i ${params.n}_vs_underhill_eval_${params.evalueTaxAnnot}.tsv \
//         -d $underhill \
//         -u  $otu5 \
//         -t $underhill_taxonomy  \
//         -o ${params.n}_vs_underhill_annotation_eval_${params.evalueTaxAnnot}.tsv \
//         -ob ${params.n}_vs_underhill_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         -dtype underhill
        
//         $biom add-metadata \
//         -i $count3 \
//         -o ${params.n}_underhill_eval_${params.evalueTaxAnnot}.biom \
//         --observation-metadata-fp ${params.n}_vs_underhill_annotation_eval_${params.evalueTaxAnnot}.biomtsv \
//         --observation-header id,taxonomy \
//         --sc-separated taxonomy \
//         --output-as-json
// 		"""
// 	}
// }


// // Mixing all annotations into the same channel
// annotation = annotation_rdp.mix(annotation_greengenes,annotation_silva,annotation_ITS)


// // Phylogeny with FastTree
// process Phylogeny {
// 	tag "$params.n"
// 	publishDir readsDir, mode:'copy', pattern: "${params.n}_otu_*"
// 	publishDir logDir, mode:'copy', pattern: "log_*.txt"

// 	input:
// 	file otu6
// 	file annot from annotation

// 	output:
// 	file "${params.n}_otu_*"
// 	file "log_*.txt"
	
// 	script:
// 	annot_db= (annot.name =~ /_vs_(rdp|greengenes|silva|findley|unite|underhill)_annotation.*\.tsv/)[0][1]
// 	"""
// 	python $extract_fasta \
// 	-d $otu6 \
// 	-i $annot \
// 	-o ${params.n}_otu_${annot_db}.fasta
	
// 	$mafft \
// 	--adjustdirectionaccurately \
// 	--thread 1 \
// 	--genafpair \
// 	--maxiterate 1000 \
// 	--ep 0  \
// 	${params.n}_otu_${annot_db}.fasta \
// 	> ${params.n}_otu_${annot_db}.ali \
// 	2> log_mafft_${params.n}_${annot_db}.txt
	
// 	$BMGE \
// 	-i ${params.n}_otu_${annot_db}.ali \
// 	-t DNA \
// 	-m ID \
// 	-h 1 \
// 	-g $params.conservedPosition \
// 	-w 1 \
// 	-b 1 \
// 	-of ${params.n}_otu_${annot_db}_bmge.ali
	
// 	$FastTreeMP \
// 	-nt ${params.n}_otu_${annot_db}_bmge.ali \
// 	> ${params.n}_otu_${annot_db}_bmge.ali.treefile \
// 	2> log_fasttree_${annot_db}.txt	
// 	"""
// }


