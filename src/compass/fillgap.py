'''
Created on 2018年6月13日

@author: Dr.liu
'''

import re,sys,os,copy,pickle
from NGS.BasicUtil import Util, VCFutil,Caculators

if len(sys.argv)<6:
    print("python fillgap.py [ref.fa] [vcf] [gap.bed] [scored.snp]  [formatedOutnamepre] [winsize]")
    exit(-1)
chrlenm={
"1":42145699,
'2':49200776,
'3':50652576,
"4":40408058,
"5":47253416,
"6":36015257,
"7":35964515,
"8":40690061,
"9":58970518}
winsize=int(sys.argv[6])
reff=open(sys.argv[1],'r')
try:
    refidx=pickle.load(open(sys.argv[1]+".myfasteridx",'rb'))
except IOError:
    Util.generateFasterRefIndex(sys.argv[1], sys.argv[1]+".myfasteridx")
    refidx=pickle.load(open(sys.argv[1]+".myfasteridx",'rb'))
vcftools="vcftools"
gapf=open(sys.argv[3],'r')
scoredsnp=open(sys.argv[4],'r')
scoredsnp.readline()
sitesingap=open(sys.argv[5],'w')
if __name__ == '__main__':
    win = Util.Window()
    i=0
    interferf=open(sys.argv[5]+".InterferingTEMP",'w')
    for gapregion in gapf:
        i+=1
        filledsites=[]
        gaplist=re.split(r"\s+",gapregion.strip())
        if not os.path.exists(sys.argv[5]+"temp"+str(i)+".recode.vcf"):
            os.system(vcftools+" --vcf "+sys.argv[2]+" --recode --recode-INFO-all --remove-indv DSW33216 --chr "+gaplist[0] +" --from-bp "+ str(gaplist[1]) +" --to-bp "+ str(gaplist[2]) + " --out "+sys.argv[5]+"temp"+str(i))
        vcfobj=VCFutil.VCF_Data(sys.argv[5]+"temp"+str(i)+".recode.vcf")
        vcflist=vcfobj.getVcfListByChrom(gaplist[0],MQfilter=0)
        
        findtagcaculator=Caculators.CaculatorToFindTAGs(mod="randomvcf",Interferingf=interferf)
        findtagcaculator.curchrom=gaplist[0]
        win.slidWindowOverlap(vcflist,int(gaplist[2]),winsize,winsize,findtagcaculator,int(gaplist[1]))
        filledsites=copy.deepcopy(win.winValueL)
        for s,e, n,poss in filledsites:
            if poss[0]!="NA":
                RefSeqMap = Util.getRefSeqBypos_faster(reff,refidx,gaplist[0],poss[0][0]-35,poss[0][0]+35,chrlenm[gaplist[0]])
                if RefSeqMap[gaplist[0]][36].upper() !=poss[0][-2]:
                    print("".join(RefSeqMap[gaplist[0]][1:]),file=open("inconsistent.txt",'a'))
                print("Setaria italica","Chr"+gaplist[0]+"_"+ str(poss[0][0]),"".join(RefSeqMap[gaplist[0]][1:36]).upper()+"["+RefSeqMap[gaplist[0]][36].upper()+"/"+poss[0][-1]+"]"+"".join(RefSeqMap[gaplist[0]][37:]).upper(),"1","Standard",gaplist[0],str(poss[0][0]),poss[0][-2],poss[0][-1],sep="\t",file=sitesingap)
        sitesingap.flush()
    interferf.close()
    FFF=open(sys.argv[5]+".InterferingTEMP",'r')
    OFO=open(sys.argv[5]+".Interfering",'w')
    for line in FFF:
        l=re.split(r"\s+",line.strip())
        if float(l[2])>0.01:
            RefSeqMap = Util.getRefSeqBypos_faster(reff,refidx,l[0],int(l[1])-35,int(l[1])+35,chrlenm[l[0]])
            print("Setaria italica","Chr"+l[0]+"_"+ str(l[1]),"".join(RefSeqMap[l[0]][1:36]).upper()+"["+l[5].upper()+"/"+l[6].upper()+"]"+"".join(RefSeqMap[l[0]][37:]).upper(),sep="\t",file=OFO)
    #gap fill over
    for oldmarker in scoredsnp:
        print(oldmarker.strip(),file=sitesingap)
#         oldmarkerlist=re.split(r"\s+",oldmarker.strip())
#         print("Setaria italica",oldmarkerlist[0],oldmarkerlist[1],"2","Standard",oldmarkerlist[4],oldmarkerlist[5],sep="\t",file=sitesingap)
    
    f=open("diao603GWASng.sorted.snp",'r')
    for line in f:
        linelist=re.split(r'\s+',line.strip())
        RefSeqMap = Util.getRefSeqBypos_faster(reff,refidx,linelist[0],int(linelist[1])-35,int(linelist[1])+35,chrlenm[linelist[0]])
        print("Setaria italica","Chr"+linelist[0]+"_"+linelist[1],"".join(RefSeqMap[linelist[0]][1:36]).upper()+"["+linelist[-3].upper()+"/"+linelist[-4]+"]"+"".join(RefSeqMap[linelist[0]][37:]).upper(),"1","Standard",linelist[0],linelist[1],linelist[-3],linelist[-4],sep="\t",file=sitesingap)
    sitesingap.close();f.close();scoredsnp.close();OFO.close()