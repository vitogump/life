'''
Created on 2017年11月14日

@author: liurui
'''
"""
plink -file filtered_OutDSW33216_chr6 --chr 6 --out filtered_OutDSW33216_chr6_1 --from-bp 0 --to-bp 360000 --recode
sed 's/-9/0/' filtered_OutDSW33216_chr6_1.ped>temp
mv temp filtered_OutDSW33216_chr6_1.ped
awk 'BEGIN{OFS="\t"}{print $1,$4}' filtered_OutDSW33216_chr6_1.map>filtered_OutDSW33216_chr6_1.info
nohup java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar ~/software/Haploview.jar -nogui -pedfile filtered_OutDSW33216_chr1.ped -info filtered_OutDSW33216_chr1.info -blockoutput GAB -log filtered_OutDSW33216_chr1.log -out filtered_OutDSW33216_r2_6_chr1 -memory 1240000 -tagrsqcutoff 0.6 -aggressiveTagging > chr1.log 2>&1 &

"""

from multiprocessing.dummy import Pool
from optparse import OptionParser
import os
import re

from NGS.BasicUtil import VCFutil


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
# parser.add_option("-p", "--prefix", dest="prefix",help="prefix1.map prefix2.map ....,prefix1.ped,prefix2.ped,")
parser.add_option("-v", "--vcffile", dest="vcffile", help="vcffilename")
parser.add_option("-s", "--sizeOfchip", dest="sizeOfchip", help="400000")
parser.add_option("-t", "--threads", dest="threads", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-o", "--outputfilename", dest="outputfilename",help="chromosome")
parser.add_option("-r", "--rmindvd", dest="rmindvd", help="DSW33216")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
NUMBER=10
def selectTAGsnp(filename):
    os.system("java     -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+options.outputfilename+filename+" -memory 80000 -maxDistance 500 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
if __name__ == '__main__':
    #count total size of 
    genomesnpspansize=0
    vcfobj=VCFutil.VCF_Data(options.vcffile)
#     print(vcfobj.VcfIndexMap)
    vcfFile=open(options.vcffile,'r')
    for curchr in vcfobj.chromOrder:
        vcfFile.seek(vcfobj.VcfIndexMap[curchr][0])
        firstline=vcfFile.readline().strip()
#         print(firstline)
        startlinelist=re.split(r'\s+',firstline)
        startpos = int(startlinelist[1].strip())
        
        vcfFile.seek(vcfobj.VcfIndexMap[curchr][1]-10000)
        lastline=vcfFile.readline().strip()
        if re.search(r"^"+curchr+"\t",lastline)==None:
            lastline=vcfFile.readline().strip()
        endlinelist=re.split(r'\s+',lastline)
        endpos = int(endlinelist[1].strip())
        
#         vcfrecOfcurchr=vcfobj.getVcfListByChrom(curchr,  MQfilter=0)
        genomesnpspansize+=(endpos-startpos)
        winsize=int(genomesnpspansize/int(options.sizeOfchip))
        sizetoSelectTAG=genomesnpspansize/(int(options.threads)*NUMBER)
    print(genomesnpspansize,winsize)
    #cut into snp record pieces whose amount equal threads * NUMBER to select TAGs . size that each piece span should bigger than winsize
    mappedlist=[]
    for curchr in vcfobj.chromOrder:
        vcfFile.seek(vcfobj.VcfIndexMap[curchr][0])
        firstline=vcfFile.readline().strip()
#         print(firstline)
        startlinelist=re.split(r'\s+',firstline)
        startpos = int(startlinelist[1].strip())
        
        vcfFile.seek(vcfobj.VcfIndexMap[curchr][1]-10000)
        lastline=vcfFile.readline().strip()
        if re.search(r"^"+curchr+"\t",lastline)==None:
            lastline=vcfFile.readline().strip()
        endlinelist=re.split(r'\s+',lastline)
        endpos = int(endlinelist[1].strip())
        i=1
        while startpos<endpos:
            if not os.path.exists(re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".ped"):
                os.system("vcftools --vcf "+options.vcffile+" --recode --recode-INFO-all --remove-indv DSW33216 --chr "+curchr +" --from-bp "+ str(startpos) +" --to-bp "+ str(startpos+sizetoSelectTAG) + " --out "+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
                os.system("vcftools --vcf "+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".recode.vcf"+" --out "+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+"  --plink")
            os.system("""awk 'BEGIN{OFS="\t"}{print $2,$4,$1}' """+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".map > "+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".info")
                
            mappedlist.append(re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
            i+=1;startpos+=sizetoSelectTAG
    # run haplovew to selectTAG    
    for j in range(0,len(mappedlist),NUMBER):
        pool=Pool(int(len(options.threads)))        
        pool.map(selectTAGsnp,mappedlist[j:j+NUMBER])
        pool.close()
        pool.join()
    #extract two tags by winsize, if no enough TAG in a win then seleced two snps whose AF approxmate to 0.5