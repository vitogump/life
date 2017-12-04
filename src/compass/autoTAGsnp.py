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

import copy
from functools import reduce
from optparse import OptionParser
import os, math, re
import sys
# import threading
from time import ctime

from NGS.BasicUtil import VCFutil, Util, Caculators


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
# parser.add_option("-p", "--prefix", dest="prefix",help="prefix1.map prefix2.map ....,prefix1.ped,prefix2.ped,")
parser.add_option("-c", "--chrlist", dest="chrlist",help="if it is not suppled , then use vcffile's information")
parser.add_option("-v", "--vcffile", dest="vcffile", help="vcffilename")
parser.add_option("-s", "--sizeOfchip", dest="sizeOfchip", help="400000")
parser.add_option("-N", "--numbertosplic", dest="numbertosplic", help="DSW33216")
parser.add_option("-t", "--threads", dest="threads", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-o", "--outputfilename", dest="outputfilename",help="chromosome")
parser.add_option("-r", "--rmindvd", dest="rmindvd", help="DSW33216")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
print(ctime())
outvcfmappedPRE=re.search(r'[^/]*$',options.outputfilename).group(0)
tempfile=open(outvcfmappedPRE+"numberinfo.txt",'w')
NUMBER=math.ceil(int(options.numbertosplic)/int(options.threads))
def splicVcfbyChr(curchr):
    mappedlistOfoneChrmOrdered=[]
    vcfFile.seek(vcfobj.VcfIndexMap[curchr][0])
    firstline=vcfFile.readline().strip()
    startlinelist=re.split(r'\s+',firstline)
    startpos = int(startlinelist[1].strip())
    
    vcfFile.seek(vcfobj.VcfIndexMap[curchr][1]-10000)#10000 will different for different vcf ,different loci
    lastline=vcfFile.readline().strip()
    if re.search(r"^"+curchr+"\t",lastline)==None:
        lastline=vcfFile.readline().strip()
    endlinelist=re.split(r'\s+',lastline)
    endpos = int(endlinelist[1].strip())
    i=1
    outvcfmappedPRE=re.search(r'[^/]*$',options.outputfilename).group(0)
    while startpos<endpos:
        if True or not os.path.exists(re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".ped"):
            os.system("vcftools --vcf "+options.vcffile+" --recode --recode-INFO-all --remove-indv DSW33216 --chr "+curchr +" --from-bp "+ str(startpos) +" --to-bp "+ str(startpos+sizetoSelectTAG) + " --out "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
            os.system("vcftools --vcf "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".recode.vcf"+" --out "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+"  --plink")
        os.system("""awk 'BEGIN{OFS="\t"}{print $2,$4,$1}' """+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".map > "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".info")
            
        mappedlistOfoneChrmOrdered.append(re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
        i+=1;startpos+=sizetoSelectTAG
    return mappedlistOfoneChrmOrdered
def selectTAGsnp(filename):
    print("java -XX:-UseGCOverheadLimit  -Xmx124g  -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+options.outputfilename+filename+" -memory 120000 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
    try:
        os.system("java -XX:-UseGCOverheadLimit  -Xmx124g  -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+options.outputfilename+filename+" -memory 120000 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
    except:
        print("java -XX:-UseGCOverheadLimit  -Xmx156g  -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+options.outputfilename+filename+" -memory 150000 -maxDistance 200 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
        os.system("java -XX:-UseGCOverheadLimit  -Xmx156g  -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+options.outputfilename+filename+" -memory 150000 -maxDistance 200 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
        return
if __name__ == '__main__':
    #count total size of 
    genomesnpspansize=0
    vcfobj=VCFutil.VCF_Data(options.vcffile)
#     print(vcfobj.VcfIndexMap)
    vcfFile=open(options.vcffile,'r')
    chrmap={}
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
        chrmap[curchr]=(startpos,endpos)
#         vcfrecOfcurchr=vcfobj.getVcfListByChrom(curchr,  MQfilter=0)
        genomesnpspansize+=(endpos-startpos)
    winsize=int(genomesnpspansize/int(options.sizeOfchip))
    sizetoSelectTAG=genomesnpspansize/(int(options.numbertosplic))
    print("genomesizewithSNP:",genomesnpspansize,". for the sizeOfchip every",winsize,"should have a tag SNP. cut genome into ",sizetoSelectTAG,"to select TAG (ie run haploview)\n",NUMBER,file=tempfile)
    tempfile.close()
    #cut into snp record pieces whose amount equal threads * NUMBER to select TAGs . size that each piece span should bigger than winsize
#     mappedlistordered=[]
#     for curchr in vcfobj.chromOrder:
        #
#     mappedlistordered=[]
#     for chr_idx in range(0,len(vcfobj.chromOrder),NUMBER):
        
    pool=Pool(int(options.threads)) 
    mappedlistordered=pool.map(splicVcfbyChr,vcfobj.chromOrder)
    pool.close()
    pool.join()
#         mappedlistordered+=t_mappedlistordered
#         print(len(mappedlistordered),mappedlistordered,sep="\n")
#         sys.stdout.flush()
#     mappedlistordered=reduce(lambda x,y:x+y,map(splicVcfbyChr,vcfobj.chromOrder))
    print(NUMBER,len(mappedlistordered),mappedlistordered,sep="\n")
    sys.stdout.flush()
    #temp . this file has been TAGed as a test
    # run haplovew to selectTAG
#     for i in range(0,)
#     for j in range(0,len(mappedlistordered),NUMBER):
#         t_list=[]
#         for i in range(j,j+NUMBER):
#             if i<len(mappedlistordered):
#                 t_list.append(threading.Thread(target=selectTAGsnp,args=(mappedlistordered[i],)))
#                 t_list[-1].setDaemon(True)
#                 t_list[-1].start()
#         else:
#             print(j,i)
#         for t in t_list:
#             t.join()
    
    for j in range(0,len(mappedlistordered),int(options.threads)):
        pool=Pool(int(options.threads))        
        pool.map(selectTAGsnp,mappedlistordered[j:j+int(options.threads)])
        pool.close()
        print("thread from",j,"to",j+int(options.threads))
        pool.join()
    print("finish haploview TAGing")    
    #extract two tags by winsize, if no enough TAG in a win then seleced two snps whose AF approxmate to 0.5
    TAGSNP={}
    for tagfile in mappedlistordered:
        tagfilename=options.outputfilename+tagfile+".TAGS"
        if not os.path.exists(tagfilename):
            continue
        with open(tagfilename,'r') as tf:
            for line in tf:
                if re.search(r'^Test\s+Alleles\s+Captured',line.strip())!=None:
                    break
            for line in tf:
                snpID=re.split(r":",re.split(r"\s+",line.strip())[0])
                if snpID[0] in TAGSNP:
                    TAGSNP[snpID[0]].append(int(snpID[1]))
                else:
                    snppos=re.search(r'^\d+',snpID[1]).group(0)
                    TAGSNP[snpID[0]]=[int(snppos)]
    tagsMap={}           
    with open(options.outputfilename+".ALLTAGS",'w') as atf:
        for chrom in sorted(TAGSNP.keys()):
            TAGSNP[chrom].sort()
            tagsMap[chrom]=[]
            for pos in TAGSNP[chrom]:
                tagsMap[chrom].append((pos,"."))
                print(chrom,pos,file=atf)
    win = Util.Window()
    selectedTAGsnps={}
    for c_chrom in chrmap.keys():
        selectedTAGsnps[c_chrom]=[]            
        findtagcaculator=Caculators.CaculatorToFindTAGs(c_chrom,vcfobj,winsize,winsize,chrmap[c_chrom][0])
        win.slidWindowOverlap(tagsMap[c_chrom],chrmap[c_chrom][1],winsize,winsize,findtagcaculator,chrmap[c_chrom][0])
        selectedTAGsnps[c_chrom]=copy.deepcopy(win.winValueL)
    with open(options.outputfilename+".selectedTAGS",'w') as atf:
        for c_chrom in selectedTAGsnps.keys():
            for i in range(len(selectedTAGsnps[c_chrom])):
                print(c_chrom + "\t" + str(i) + "\t" + str(selectedTAGsnps[c_chrom][i][2]) + "\t" + str(selectedTAGsnps[c_chrom][i][3][0]) + "\t" +str(selectedTAGsnps[c_chrom][i][3][1]), file=atf)
    vcfFile.close()