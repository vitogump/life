picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp -jar $picard_tools_path/SortSam.jar TMP_DIR=/home/liurui/tmp INPUT=${.bam} OUTPUT=${output=/home/liurui/data/bamfiles/OASLref|suffix=sorted.bam} SORT_ORDER=coordinate
