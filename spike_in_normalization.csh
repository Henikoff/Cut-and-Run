#!/bin/csh 

#Cut-and-Run spike-in normalization example C-shell script
#Seven arguments
#                   1           2            3        4                  5             6      7
#spike_norm.csh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens_file  min_len max_len

#Assumes Illumina paired-end sequencing
#Packages used: bowtie2 or other alignment program, picard and bedtools
#Requires bed files of alignments to the experimental genome as well as to the spike-in genome
#If there is no spike-in, put "none" in argument 2.
#The only information used from the spike-in alignment is the number of fragments aligned
#(the number of lines in spike_genome.bed).
#min_len and max_len refer to the fragment lengths in genome.bed

#Format for genome_chr_lens_file is (this information is in the sam file headers):
#chr1	249250621
#chr2	243199373
#...

#Before running this script:
#1. Align Illumina fastq read files to both genomes producing two sam files.
#   Any alignment program that produces a sam/bam file can be used.
#   However the alignments are done, be careful that not too many reads are aligned
#   to both genomes.
#2. Optionally remove PCR and optical duplicates (use picard MarkDuplicates) from the sam file.
#3. Extract properly aligned fragments from both genomes.
#	a. Convert sam to bam format (picard or samtools)
#	b. Extract aligned fragments from bam file (bedtools bamtobed)
#		For subsequent processing, the length of each fragment is
#		required in the bed file, this will probably need to be
#		added. awk can be used to add fragment length, e.g.:
#		cat bamtobed.bed | awk -v OFS='\t' '{len = $3 - $2; print $0, len }'
#		Assuming that $2 is the actual start position - 1

#----------------------------------------------------------------------------------------------
#4. spike-in normalization

if ($#argv < 7) then
	echo "USAGE spike_norm.csh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens min_len max_len"
	echo "Spike-in calibration using bedtools genomecov"
	echo "scale is an arbitrary large number used as multiplier (e.g. 10000)"
	echo "If there is no spike-in, use spike_genome.bed = none"
	echo "min_len and max_len refer to the lengths of fragments in genome.bed to normalize"
	echo "To normalize all fragments use min_len = 1 and max_len = 1000"
	echo "Output will be placed in the current directory"
	exit(-1)
endif

set genome_bed = $1
set spike_bed = $2
set scale = $3
set report = $4
set genome_len = $5
set min_len = $6
set max_len = $7

echo $genome_bed $spike_bed $scale $report $min_len $max_len
if (!(-e $genome_bed) || (-z $genome_bed)) then
	echo "$genome_bed not found or is empty"
	exit(-1)
endif
if (!(-e $genome_len) || (-z $genome_len)) then
	echo "$genome_len not found or is empty"
	exit(-1)
endif

set temp = $genome_bed:t
set name = $temp:r
set output = $name.${min_len}-${max_len}.$report
echo "Output is in file $output in the current directory"
if (-e $spike_bed && !(-z $spike_bed)) then
        set spike_count = `wc -l $spike_bed | awk '{print $1}'`
        set scale_factor = `echo "$scale / $spike_count" | bc -l`
else
        set spike_count = 0
        set scale_factor = $scale
endif
echo scale_factor=$scale_factor

#Select fragments within the length range, assumes fragment length is in column 4 of the bed file
cat $genome_bed | awk -v min=$min_len -v max=$max_len '{if ($4 >= min && $4 <= max) print}' > $$.temp.bed

#Use genomecov to compute normalized genome coverage
#see http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
#The first position is start-1 and the last is end for bed files, -bg and -bga
#-bga prints zero intervals -bg doesn't, -d prints each bp starting from 1 (not 0)

#This is for IGV which doesn't recognize .bg or .bga files
if ($report == "bg" || $report == "bga") then
	set output = $name.${min_len}-${max_len}.bedgraph
	echo track type=bedGraph name=$name > $output
endif
bedtools genomecov -$report -scale $scale_factor -i $$.temp.bed -g $genome_len >> $output

unalias rm
rm $$.*
exit
