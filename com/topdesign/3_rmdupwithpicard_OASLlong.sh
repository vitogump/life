picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp -jar $picard_tools_path/MarkDuplicates.jar INPUT=${.sorted.bam} OUTPUT=${output=/home/liurui/data/bamfiles/OASLref|suffix=sorted.dedup.bam} METRICS_FILE=${output=/home/liurui/data/bamfiles/OASLref|suffix=sorted.dedup.metrics} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 REMOVE_DUPLICATES=true
