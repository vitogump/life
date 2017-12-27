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
# import threading

import copy
from functools import reduce
from multiprocessing.dummy import Pool
from optparse import OptionParser
import os, math, re
import shlex  
import signal
import sys, subprocess
from time import ctime, time

from NGS.BasicUtil import VCFutil, Util, Caculators


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
# parser.add_option("-p", "--prefix", dest="prefix",help="prefix1.map prefix2.map ....,prefix1.ped,prefix2.ped,")
parser.add_option("-c", "--chrlist", dest="chrlist",help="if it is not suppled , then use vcffile's information")
parser.add_option("-v", "--vcffile", dest="vcffile", help="vcffilename")
parser.add_option("-s", "--sizeOfchip", dest="sizeOfchip", help="400000")
parser.add_option("-N", "--numbertosplic", dest="numbertosplic", help="800 ? 300? ...")
parser.add_option("-t", "--threads", dest="threads", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-o", "--outputfilename", dest="outputfilename",help="chromosome")
parser.add_option("-r", "--rmindvd", dest="rmindvd",default=None, help="includeTagsFilestr")
parser.add_option("-i", "--includeTag", dest="includeTag", help="DSW33216")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

"need vcftools java Haploview.jar"
vcftools="vcftools"
Haploview="~/software/Haploview.jar"
(options,args)=parser.parse_args()
print(ctime());sys.stdout.flush()
if options.includeTag!=None:
    includeTagsFilestr=" -includeTagsFile "+options.includeTag
else:
    includeTagsFilestr=""
outvcfmappedPRE=re.search(r'[^/]*$',options.outputfilename).group(0)
outvcfmappedpath=os.path.dirname(options.outputfilename)
tempfile=open(outvcfmappedPRE+"numberinfo.txt",'w')
NUMBER=math.ceil(int(options.numbertosplic)/int(options.threads))
regionMap={}
def splicVcfbyChr(curchr):
    mappedlistOfoneChrmOrdered=[]
    vcfFile=open(options.vcffile,'r')
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
    while startpos<endpos:
        if  not os.path.exists(outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".ped"):
            os.system(vcftools+" --vcf "+options.vcffile+" --recode --recode-INFO-all --remove-indv DSW33216 --chr "+curchr +" --from-bp "+ str(startpos) +" --to-bp "+ str(startpos+sizetoSelectTAG) + " --out "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
            os.system(vcftools+" --vcf "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".recode.vcf"+" --out "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+"  --plink")
            os.system("""awk 'BEGIN{OFS="\t"}{print $2,$4,$1}' """+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".map > "+outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".info")
        else:
            print("use exist file",outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)+".ped")
        regionMap[outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i)]=(startpos,startpos+sizetoSelectTAG) 
        print(outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
        print(startpos,startpos+sizetoSelectTAG)
        mappedlistOfoneChrmOrdered.append(outvcfmappedPRE+re.search(r"^[^.]*",re.search(r"[^/]*$",options.vcffile).group(0)).group(0)+"chr"+curchr+"_"+str(i))
        i+=1;startpos+=sizetoSelectTAG
    vcfFile.close()
    return mappedlistOfoneChrmOrdered
def selectTAGsnp(filename):

    #this block is for 800 fen only
    mainchrno=re.search(r'(\d+)_(\d+)$',filename).group(1)

    
    x8end=int(re.search(r'_(\d+)$',filename).group(1)) * sizetoSelectTAG
    x8start=(int(re.search(r'_(\d+)$',filename).group(1))-1)*sizetoSelectTAG
    x3splitno=math.ceil( x8end/1336284.436)
    X3mod=x8end%1336284.436

    ##########and the if block below##########################

#     if re.search(r'[68]$',mainchrno)==None:
#         print("skip other chrom:",filename)
#         return  
    """It is recommended that Haploview be run on a machine with at least 128M of memory. The Haploview
jarfile should now automatically allocate extra memory when starting up, so the -Xmx flag is no longer
required when running the program from the command line.
    """      
#     print("java -XX:-UseGCOverheadLimit   -jar "+Haploview+" -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+outvcfmappedpath+"/"+filename+" -memory 124000 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
#     try:
    if (os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS") or (os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(int( x8end/1336284.436))) and X3mod>=winsize))  and (x8start>=(x3splitno-1)*1336284.436 or os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno-1)+".TAGS")):#this if block is specific program for 800 fen, which were used to complement the running of 300 shares split 
        print("skip"+ outvcfmappedpath+"/"+filename+".TAGS as corresponding 300fen exist:"+"first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS")
        sys.stdout.flush()
        return########blockend 
    elif not os.path.exists(outvcfmappedpath+"/"+filename+".TAGS") and (not os.path.exists(outvcfmappedpath+"/"+filename+"_1.TAGS")) :
        tagcomdstr="java -XX:-UseGCOverheadLimit   -jar "+Haploview+" -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+outvcfmappedpath+"/"+filename+" -memory 50000  -minMAF 0.05 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.8 -pairwiseTagging "+includeTagsFilestr
#                    "java -XX:-UseGCOverheadLimit   -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+outvcfmappedpath+"/"+filename+" -memory 124000 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging" 
#         tagcomd=shlex.split(tagcomdstr)
        print(tagcomdstr)
        p=subprocess.Popen(tagcomdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        while p.poll() is None:
            line=p.stdout.readline()
            if "Exception" in str(line) or "Error" in str(line):
                print(line,p.pid,filename)
                try:
                    p.terminate()
                    p.kill()
                    k=os.kill(p.pid,signal.SIGTERM)
                    print(k)
                except:
                    print("still runing?")
                
                print(line)

                try:
                    pass
                except:
                    pass
            else:
                print(line)
            print(p.stdout.readline())
        try:
            if p.returncode == 0:
                return
#             elif a==0 and b==0:
#                 return
            else:
                os.system(vcftools+ " --vcf "+filename+".recode.vcf --recode --recode-INFO-all --chr "+ mainchrno+" --from-bp "+str(regionMap[filename][0]) +" --to-bp "+str((regionMap[filename][1]+regionMap[filename][0])/2)+" --out "+filename+"_1")
                os.system(vcftools+ " --vcf "+filename+"_1"+".recode.vcf"+" --out "+filename+"_1 --plink")
                os.system("""awk 'BEGIN{OFS="\t"}{print $2,$4,$1}' """+filename+"_1.map >" +filename+"_1.info")
                
                os.system(vcftools+ " --vcf "+filename+".recode.vcf --recode --recode-INFO-all --chr "+ mainchrno+" --from-bp "+str((regionMap[filename][1]+regionMap[filename][0])/2) +" --to-bp "+str(regionMap[filename][1])+" --out "+filename+"_2")
                os.system(vcftools+ " --vcf "+filename+"_2"+".recode.vcf"+" --out "+filename+"_2 --plink")
                os.system("""awk 'BEGIN{OFS="\t"}{print $2,$4,$1}' """+filename+"_2.map >" +filename+"_2.info")
                print("start splited threads",filename)
                a=os.system("java -XX:-UseGCOverheadLimit  -jar "+Haploview+" -nogui -pedfile "+filename+"_1.ped -info "+filename+"_1.info  -blockoutput GAB -log "+filename+"_1.log -out "+outvcfmappedpath+"/"+filename+"_1   -memory 150000 -maxDistance 150 -taglodcutoff 3 -tagrsqcutoff 0.6 -pairwiseTagging "+includeTagsFilestr)
                b=os.system("java -XX:-UseGCOverheadLimit  -jar "+Haploview+" -nogui -pedfile "+filename+"_2.ped -info "+filename+"_2.info  -blockoutput GAB -log "+filename+"_2.log -out "+outvcfmappedpath+"/"+filename+"_2   -memory 150000 -maxDistance 150 -taglodcutoff 3 -tagrsqcutoff 0.6 -pairwiseTagging "+includeTagsFilestr)
                print("should run again? mv the split block here? os.system(vcftools+...)")
        except:
            print(filename+"_1.ped or"+filename+"_1.ped with exception" )
#         a=os.system("java -XX:-UseGCOverheadLimit   -jar ~/software/Haploview.jar -nogui -pedfile "+filename+".ped -info "+filename+".info -blockoutput GAB -log "+filename+".log -out "+outvcfmappedpath+"/"+filename+" -memory 124000 -maxDistance 300 -taglodcutoff 3 -tagrsqcutoff 0.6 -aggressiveTagging")
        
    else:
        print(outvcfmappedpath+"/"+filename+".TAGS exist")
        a=0
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
    print("genomesizewithSNP:",genomesnpspansize,". for the sizeOfchip every",winsize,"should have a tag SNP. cut genome into ",sizetoSelectTAG,"to select TAG (ie run haploview)\n need batches",NUMBER,file=tempfile)
    tempfile.close()

        
    pool=Pool(int(options.threads)) 
    mappedlistordered=reduce(lambda x,y:x+y,pool.map(splicVcfbyChr,vcfobj.chromOrder))
    pool.close()
    pool.join()


    print(NUMBER,len(mappedlistordered),mappedlistordered,sep="\n")
    sys.stdout.flush()
    pool=Pool(int(options.threads))        
    pool.map(selectTAGsnp,mappedlistordered)
    pool.close()

    pool.join()
    print("finish haploview TAGing")    
    #extract two tags by winsize, if no enough TAG in a win then seleced two snps whose AF approxmate to 0.5
    TAGSNP={}
    for tagfile in mappedlistordered:
        tflist=[]
        tagfilename=outvcfmappedpath+"/"+tagfile+".TAGS"
        if not os.path.exists(tagfilename):
            mainchrno=re.search(r'(\d+)_(\d+)$',tagfile).group(1)
            x8end=int(re.search(r'_(\d+)$',tagfile).group(1)) * sizetoSelectTAG
            x8start=(int(re.search(r'_(\d+)$',tagfile).group(1))-1)*sizetoSelectTAG
            x3splitno=math.ceil( x8end/1336284.436)
            X3mod=x8end%1336284.436
#             (os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS") or (os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(int( x8end/1336284.436))) and X3mod>=winsize))  and (x8start>=(x3splitno-1)*1336284.436 or os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno-1)+".TAGS"))
            
            if os.path.exists(outvcfmappedpath+"/"+tagfile+"_1.TAGS") or os.path.exists(outvcfmappedpath+"/"+tagfile+"_2.TAGS"):
                if os.path.exists(outvcfmappedpath+"/"+tagfile+"_1.TAGS"):
                    tflist.append(open(outvcfmappedpath+"/"+tagfile+"_1.TAGS",'r'))
                elif os.path.exists(outvcfmappedpath+"/"+tagfile+"_2.TAGS"):
                    tflist.append(open(outvcfmappedpath+"/"+tagfile+"_2.TAGS",'r'))
            elif os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS"):
                tflist.append(open(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS",'r'))
                print("use "+outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno)+".TAGS")
            elif os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(int( x8end/1336284.436))) and X3mod>=winsize:
                tflist.append(open(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(int( x8end/1336284.436)),'r'))
                print("use "+outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(int( x8end/1336284.436)))

            
            else:
                print(tagfilename,"does not exist and no complement file")
                continue
            #
            if not (os.path.exists(outvcfmappedpath+"/"+tagfile+"_1.TAGS") or os.path.exists(outvcfmappedpath+"/"+tagfile+"_2.TAGS")) and (x8start<(x3splitno-1)*1336284.436 and os.path.exists(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno-1)+".TAGS")):
                print("use "+outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno-1)+".TAGS")
                tflist.append(open(outvcfmappedpath+"/first_200ksites_300fenfilteredchr"+str(mainchrno)+"_"+str(x3splitno-1)+".TAGS",'r'))
            else:
                tfe=None
        else:
            tflist.append(open(tagfilename,'r'))
        for tf in tflist:
            for line in tf:
                if re.search(r'^Test\s+Alleles\s+Captured',line.strip())!=None:
                    break
            for line in tf:
                snpID=re.split(r":",re.split(r"\s+",line.strip())[0])
                if snpID[0] in TAGSNP:
                    TAGSNP[snpID[0]].append(int(re.split(r",",snpID[1])[0]))
                else:
                    snppos=int(re.split(r",",snpID[1])[0])
                    TAGSNP[snpID[0]]=[snppos]
            tf.close()

        
    tagsMap={}           
    with open(options.outputfilename+".ALLTAGS",'w') as atf:
        for chrom in sorted(TAGSNP.keys()):
            TAGSNP[chrom].sort()
            tagsMap[chrom]=[]
            for pos in TAGSNP[chrom]:
                tagsMap[chrom].append((pos,"."))
                print(chrom,pos,sep="\t",file=atf)
    win = Util.Window()
    selectedTAGsnps={}
    for c_chrom in sorted(chrmap.keys()):
        selectedTAGsnps[c_chrom]=[]
        findtagcaculator=Caculators.CaculatorToFindTAGs(mod="randomvcf")
        vcflist=vcfobj.getVcfListByChrom(c_chrom,MQfilter=0)#randomselecter.vcfobj.
        win.slidWindowOverlap(vcflist,chrmap[c_chrom][1],winsize,winsize,findtagcaculator,chrmap[c_chrom][0])
        findtagcaculator.mod="selectTAG"
        win.slidWindowOverlap(tagsMap[c_chrom],chrmap[c_chrom][1],winsize,winsize,findtagcaculator,chrmap[c_chrom][0])
        selectedTAGsnps[c_chrom]=copy.deepcopy(win.winValueL)
    with open(options.outputfilename+".selectedTAGS",'w') as atf:
        print("chrNo\twinNo\tNoOfTags\tTAG1\tTAG2",file=atf)
        for c_chrom in sorted(selectedTAGsnps.keys()):
            for i in range(len(selectedTAGsnps[c_chrom])):
                print(c_chrom + "\t" + str(i) + "\t" + str(selectedTAGsnps[c_chrom][i][2]) + "\t" + str(selectedTAGsnps[c_chrom][i][3][0]) + "\t" +str(selectedTAGsnps[c_chrom][i][3][1]), file=atf)
    vcfFile.close()