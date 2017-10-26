# -*- coding: UTF-8 -*-
from optparse import OptionParser
import re,string, os,pickle,random
'''
Created on 2017年9月29日

@author: liurui
'''
from NGS.BasicUtil import VCFutil


parser = OptionParser()
parser.add_option("-t","--Morganperbp",dest="Morganperbp",default="5.4696217209617786e-8")
parser.add_option("-c", "--csvfile", dest="csvfile",help="used only when -d 0")
parser.add_option("-o", "--objV", dest="objV",help="used only when -d 0")# action="callback",type="string",callback=useoptionvalue_previous1,
parser.add_option("-r","--refV", dest="refV", help="if value is 0 convert -t type -i inputfile into vcf files. or value =1  covert -v assigned vcf files into -t type files")
parser.add_option("-x", "--xp-clrpath", dest="xppath",help="used only when -d 0")
parser.add_option("-w", "--withheader", dest="withheader",action="store_true",help="if vcf file with # header ")
(options, args) = parser.parse_args()


genotypesep=" "
gWin=" 0.05 "
snpWin=" 200 "
gridSize=" 2000 "
p=" -p0 "#specify this if the genotype data is unphased.

corrLevel=" 0.95 "
sampleIDpopmapfile=open(options.csvfile,'r')
sampleIDpopmapfile.readline();sampleID_to_popmap={}
for line in sampleIDpopmapfile:
    linelist=re.split(r',',line.strip())
    sampleID_to_popmap[linelist[3].strip()]=linelist[6].strip()
    
    
chrlist=[1,2,3,4,5,6,7,8,9,10,11,12]

objpop=VCFutil.VCF_Data(options.objV)
refpop=VCFutil.VCF_Data(options.refV)

print("produce geno snp ind for obj")
VCFutil.VCF_Data.Vcf2geno_snp_ind(options.objV, sampleID_to_popmap, re.search(r"[^/]*$",options.objV).group(0)[:-4], "GATK", options.Morganperbp, "allchr", objpop.VcfIndexMap, genotypesep,options.withheader)
print("produce geno snp ind for ref")
VCFutil.VCF_Data.Vcf2geno_snp_ind(options.refV, sampleID_to_popmap, re.search(r"[^/]*$",options.refV).group(0)[:-4], "GATK", options.Morganperbp, "allchr", refpop.VcfIndexMap, genotypesep,options.withheader)
def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])
if __name__ == '__main__':
    randstr=random_str()
    tmp=open(randstr+"tempXPCLRmergefile","w")
    tmp.close()
    
    for chrom in chrlist:
        #os.popen()
        print("""awk '$2=="""+'"'+str(chrom)+'"{OFS="\t";print " "$0}'+"' "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+"allchr.snp > "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp")

        os.system("""awk '$2=="""+'"'+str(chrom)+'"{OFS="\t";print " "$0}'+"' "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+"allchr.snp > "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp")
        print("wc -l "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp")
        a=os.popen("less "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp|wc -l")
        totallines=int(a.readline().strip())
        a.close()
        lasttotalines=1
        print("sed -n '"+str(lasttotalines)+","+str(lasttotalines-1+totallines)+"p' "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+"allchr.geno >"+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".geno")
        print("sed -n '"+str(lasttotalines)+","+str(lasttotalines-1+totallines)+"p' "+re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".geno")
        os.system("sed -n '"+str(lasttotalines)+","+str(lasttotalines-1+totallines)+"p' "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+"allchr.geno >"+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".geno")
        os.system("sed -n '"+str(lasttotalines)+","+str(lasttotalines-1+totallines)+"p' "+re.search(r"[^/]*$",options.refV).group(0)[:-4]+"allchr.geno >"+re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".geno")
        
        lasttotalines=totallines+1
        print(options.xppath+" -xpclr "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".geno "+ re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".geno "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+"vs"+re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+" -wl "+gWin+snpWin+gridSize+str(chrom)+p+corrLevel)
        os.system(options.xppath+" -xpclr "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".geno "+ re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".geno "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+"vs"+re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+" -wl "+gWin+snpWin+gridSize+str(chrom)+p+corrLevel)
        os.system("cat "+randstr+"tempXPCLRmergefile "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+"vs"+re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".xpclr.txt > "+randstr+"tempXPCLRmergefileo")
        os.system("mv "+randstr+"tempXPCLRmergefileo "+randstr+"tempXPCLRmergefile")
        os.system("rm "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".geno "+ re.search(r"[^/]*$",options.refV).group(0)[:-4]+str(chrom)+".geno "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+str(chrom)+".snp ")
    os.system("mv "+randstr+"tempXPCLRmergefile "+re.search(r"[^/]*$",options.objV).group(0)[:-4]+"vs"+re.search(r"[^/]*$",options.refV).group(0)[:-4]+"allchr")