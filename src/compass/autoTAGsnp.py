'''
Created on 2017年11月14日

@author: liurui
'''
"""
nohup java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar ~/software/Haploview.jar -nogui -pedfile filtered_OutDSW33216_chr1.ped -info filtered_OutDSW33216_chr1.info -blockoutput GAB -log filtered_OutDSW33216_chr1.log -out filtered_OutDSW33216_r2_6_chr1 -memory 1240000 -tagrsqcutoff 0.6 -aggressiveTagging > chr1.log 2>&1 &
plink -file filtered_OutDSW33216_chr6 --chr 6 --out filtered_OutDSW33216_chr6_1 --from-bp 0 --to-bp 360000 --recode
"""

from optparse import OptionParser


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-c", "--chrlist", dest="chrlist",help="chromosome")
parser.add_option("-v", "--variantfilewithref", dest="variantfilewithref",action="append",nargs=2, help="vcflikefile corresponding_ref")
parser.add_option("-r", "--objectref", dest="objectref", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-s", "--chrsignal", dest="chrsignal",help="chromosome")
parser.add_option("-f", "--flanklen", dest="flanklen",default='70',help="ref fa file mode2")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
if __name__ == '__main__':
    pass