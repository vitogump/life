'''
Created on 2018年5月9日

@author: Dr.liu
'''
import os,sys,re,pickle,copy
from NGS.BasicUtil import  VCFutil
from optparse import OptionParser
from multiprocessing.dummy import Pool
from functools import reduce

parser = OptionParser()
parser.add_option("-c", "--chrommap", dest="chrommap",default=None,# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="first col corresponding vcf's, second corresponding new. all vcf's chrom is not necessary. transchr only occur in outfile")
parser.add_option("-a","--affectedlist",dest="affectedlist",help="default infile1_infile2")
parser.add_option("-u","--unaffectedlist",dest="unaffectedlist")
parser.add_option("-g","--geno",dest="geno",help="0~1")#
parser.add_option("-e", "--excludesites", dest="excludesites",help="first two cols are chr pos")
parser.add_option("-m","--maf",dest="maf")
parser.add_option("-v","--vcffile",dest="vcffile",action="append", default=[],help="default infile1_infile2")
parser.add_option("-o", "--output", dest="output", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

def mergeMapPed(totalMapPed,mapped):
    totalSnpPoslist=totalMapPed[0]+mapped[0]
    for indname in totalMapPed[1].keys():
        totalMapPed[1][indname]+=mapped[1][indname]
    return copy.deepcopy(totalSnpPoslist),copy.deepcopy(totalMapPed[1])
if __name__ == '__main__':
    spath = os.path.dirname(os.path.realpath(sys.argv[0]))
    print(spath)
    affectedlist=[];unaffectedlist=[];affectunaffectmark={}
    f=open(options.affectedlist,'r')
    for line in f:
        affectedlist.append(line.strip())
        affectunaffectmark[line.strip()]="2"
    f.close()
    f=open(options.unaffectedlist,'r')
    for line in f:
        affectunaffectmark[line.strip()]="1"
    f.close()
    
    mapfile = open(options.output+".map", "w")
    pedfile = open(options.output+ ".ped", "w")
    vcfobj=VCFutil.VCF_Data(options.vcffile[0])
    chromchangemap={};excludesitesMapBchr={}
    
    f=open(options.chrommap,'r')
    for line in f:
        linelist=re.split(r'\s+',line.strip())
        chromchangemap[linelist[0].strip()]=linelist[1].strip()
    f.close()
    f=open(options.excludesites,'r')
    for line in f:
        linelist=re.split(r'\s+',line.strip())
        if linelist[0] in excludesitesMapBchr:
            excludesitesMapBchr[linelist[0]].append(int(linelist[1]))
        else:
            excludesitesMapBchr[linelist[0]]=[int(linelist[1])]
    f.close()
    parameterstuples_list=[]
    for chromNo in vcfobj.chromOrder:
        excludesitesAchr=excludesitesMapBchr[chromNo] if chromNo in excludesitesMapBchr else []
#         parameterstuples_list.append((options.vcffile[0],"GATK",vcfobj.VcfIndexMap,True,chromNo,options.affectedlist,options.unaffectedlist,chromchangemap,float(options.geno),float(options.maf),excludesitesMapBchr[chromNo]))
        parameterstuples_list.append({"vcfFileName":options.vcffile[0],"software":"GATK","VcfIndexMap":vcfobj.VcfIndexMap,"withheader":True,"chrom":chromNo,"affectedlist":options.affectedlist,"unaffectedlist":options.unaffectedlist,"chromchangemap":chromchangemap,"geno":float(options.geno),"maf":float(options.maf),"excludesits":excludesitesAchr})
    print("vcfobj.chromOrder",len(vcfobj.chromOrder))
    pool=Pool(36)
    
#     print("parameterstuples_list",parameterstuples_list,len(parameterstuples_list[0]),len(parameterstuples_list[1]))
    positionlist,pedmap=reduce(mergeMapPed,pool.map(VCFutil.VCF_Data.Vcf2Ped_WapperForpoolthreads,parameterstuples_list))
#     pool.map(VCFutil.VCF_Data.Vcf2Ped_WapperForpoolthreads,parameterstuples_list)
    pool.close()
    pool.join()

    #collection and print final
#     for chromNo in vcfobj.chromOrder:
#         positionlist=pickle.load(open(chromNo+".positionlist","rb"))
#         pedmap=pickle.load(open(chromNo+".pedmap","rb"))
    for elem in positionlist:
        print(elem[0], elem[1], elem[2], elem[3], sep='\t', file=mapfile)
    i = 1
    for name in sorted(pedmap.keys()):
        print(i,name,"0","0","1",affectunaffectmark[name],"\t".join(pedmap[name]),sep="\t",file=pedfile)
        i+=1
    mapfile.close()
    pedfile.close()