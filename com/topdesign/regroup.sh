picard_tools_path=/home/liurui/software/picard-tools-1.119
java -Xmx40g -jar $picard_tools_path/AddOrReplaceReadGroups.jar I=/home/liurui/data/bamfiles/shaoxingqingkeegg/shaoxingpool27.dedup.realigned.mergesorted.uniqmap.bam O=/home/liurui/data/bamfiles/shaoxingqingkeegg/shaoxingpool27.dedup.realigned.mergesorted.uniqmap.regroupintoone.bam  ID=shaoxing LB=lib1 PL=illumina PU=pool SM=shaoxing27
java -Xmx40g -jar $picard_tools_path/BuildBamIndex.jar I=/home/liurui/data/bamfiles/shaoxingqingkeegg/shaoxingpool27.dedup.realigned.mergesorted.uniqmap.regroupintoone.bam
