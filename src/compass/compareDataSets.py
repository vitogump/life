# -*- coding: UTF-8 -*-
'''
Created on 2017年10月23日

@author: liurui
'''
from optparse import OptionParser
import pickle
import re

from NGS.BasicUtil import Util, VCFutil


parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa",action="append",
                  help="reference.fa")
parser.add_option("-M", "--myformatfilesuffix", dest="myformatfilesuffix", help="program must be run under this folder, myformatfilesuffix")
parser.add_option("-m", "--myformatNamelistfile", dest="myformatNamelistfile", help="variants")
parser.add_option("-V", "--vcffilename", dest="vcffilename", help="variants")

parser.add_option("-o", "--pathtooutputpre", dest="pathtooutputpre", help="default infile1_infile2")

parser.add_option("-n", "--namemapfile", dest="namemapfile")

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()

refFastaName1=options.reffa[0]
refFastaName2=options.reffa[1]
reffastaidxName1 = refFastaName1 + ".myfasteridx"
reffastaidxName2 = refFastaName2 + ".myfasteridx"
try:
    refidxByChr2 = pickle.load(open(reffastaidxName2, 'rb'))
    refidxByChr1 = pickle.load(open(reffastaidxName1, 'rb'))
except IOError:
    Util.generateFasterRefIndex(refFastaName1, reffastaidxName1)
    Util.generateFasterRefIndex(refFastaName2, reffastaidxName2)
    refidxByChr1 = pickle.load(open(reffastaidxName1, 'rb'))
    refidxByChr2 = pickle.load(open(reffastaidxName2, 'rb'))
commsample_idxlistinM=[]
commsample_idxlistinV=[]
degenerateM={"R":"AG",
             "Y":"CT",
             "M":"AC",
             "K":"GT",
             "S":"GC",
             "W":"AT",
             
             "A":"AA","T":"TT","C":"CC","G":"GG"
    }
outfile=open(options.pathtooutputpre+".compared","w")
if __name__ == '__main__':
    
    #find the same indvd between two data sets
    """
    first two col must be longer than the following col
    in this case , i put vcf name to first two col
    """
    MfileNameMap={};VfileNameMap={}
    VtoMmap={}#{vcfname:myformatname}
    nameMapFile=open(options.namemapfile,'r')
    nameMapFile.readline()
    for line in nameMapFile:
        linelist=re.split(r"\s+",line.strip())
        print(linelist)
        VfileNameMap[linelist[1]]=linelist[0];
        if len(linelist)>=4:
            MfileNameMap[linelist[3]]=linelist[2] 
        else:
            pass
    nameMapFile.close()
    mf=open(options.pathtooutputpre+"temp",'w')
    print("nameInVcf\tnameInMyFormat",file=mf)
    for k,v in MfileNameMap.items():
        if k in VfileNameMap:
            VtoMmap[VfileNameMap[k]]=v
            print(VfileNameMap[k],v,sep="\t",file=mf)
    mf.close()
    
    myformatNamelist=[]
    f=open(options.myformatNamelistfile,"r")
    for my_Sample_name in f:
        myformatNamelist.append(my_Sample_name.strip())
    f.close()
    

    #find the same indvd END
    refFastahandle1=open(refFastaName1,'r')
    refFastahandle2=open(refFastaName2,'r')
    vcfdataset=VCFutil.VCF_Data(options.vcffilename)
    for k,v in VtoMmap.items():
        commsample_idxlistinM.append(myformatNamelist.index(v)) 
        commsample_idxlistinV.append(vcfdataset.VcfIndexMap["title"].index(k))
    
    #
#     bbb=Util.getRefSeqBypos_faster(refFastahandle, refidxByChr, "1", 1, 1)
#     print(bbb)
#     exit()
    for chrom in vcfdataset.chromOrder:
        vcfRecOfAChrom=vcfdataset.getVcfListByChrom(chrom, MQfilter=None)
        MFfRecOfAChrom=[]
        #read -M file and change format
        try:
            curMyFormatfile=open('Chr0'+chrom+options.myformatfilesuffix.strip(),'r')
        except:
            continue
        for line in curMyFormatfile:
            linelist=re.split(r"\s+",line.strip())
            allindalleles=""
            for geno in linelist[2:]:
                if geno!="-":
                    allindalleles+=degenerateM[geno]
            alleles=list(set(list(allindalleles)))
            allelefreq={}
            for allele in alleles:
                allelefreq[allele]=allindalleles.count(allele)/(len(allindalleles))
            if len(alleles)==0 : 
                INFO="AC=.;AF=。;AN=0"
                MFfRecOfAChrom.append((int(linelist[1]),"-","-",INFO,"GT:AD:DP:GQ:PL",["./."]))
            elif len(alleles)==1:
                INFO="AC=1;AF=1;AN="+str(len(allindalleles))
                MFfRecOfAChrom.append((int(linelist[1]),alleles[0],"-",INFO,"GT:AD:DP:GQ:PL",["./."]))
            else:
                INFO="AC="+str(allindalleles.count(alleles[1]))+";AF="+str(allelefreq[alleles[1]])[0:5]+";AN="+str(len(allindalleles))
                MFfRecOfAChrom.append((int(linelist[1]),alleles[0],",".join(alleles[1:]),INFO,"GT:AD:DP:GQ:PL",["./."]))
            
        curMyFormatfile.close()
        #read -M file and change format end

        #outjoin those recs
        outjoinedrecs=Util.alinmultPopSnpPos_diffrefalt([{chrom:MFfRecOfAChrom},{chrom:vcfRecOfAChrom}], "o")
        REF1={chrom:[0,"x"]}
        REF2={chrom:[0,"x"]}
        #test 
        testf=open("Mtesto.txt",'w')
        for recM in outjoinedrecs[chrom]:
            print(recM,file=testf)
        testf.close()
        #test end 
        print("position\tdataMgeno\tdataMfreq\tdataVgeno\tdataMfreq\trefV1\trefV2",file=outfile)       
        for rec in outjoinedrecs[chrom]:

            REF1=Util.getRefSeqBypos_faster(refFastahandle1, refidxByChr1, chrom, rec[0], rec[0])#,seektuple=(refFastahandle1.tell(),REF1[chrom][0])
            REF2=Util.getRefSeqBypos_faster(refFastahandle2, refidxByChr2, chrom, rec[0], rec[0])#,seektuple=(refFastahandle2.tell(),REF2[chrom][0])
#             RefSeqMap=Util.getRefSeqBypos_faster(refFastahandle, refindex, chrom, rec[0]-1, rec[0], currentChromNOlen, seektuple)
            if rec[1]==None:


                print('Chr0'+chrom,str(rec[0]),"-","-","/".join(rec[2][:2]),rec[2][2][:24],REF1[chrom][1],REF2[chrom][1],sep="\t",file=outfile)#rec[2][4:]
            elif rec[2]==None:
                print('Chr0'+chrom,str(rec[0]),"/".join(rec[1][:2]),rec[1][2],"-","-",REF1[chrom][1],REF2[chrom][1],sep="\t",file=outfile)#["./."]*len(commsample_idxlistinV),
            else:
                
                print('Chr0'+chrom,str(rec[0]),"/".join(rec[1][:2]),rec[1][2],"/".join(rec[2][:2]),rec[2][2][:24],REF1[chrom][1],REF2[chrom][1],sep="\t",file=outfile)#rec[2][4:],
    refFastahandle1.close();refFastahandle2.close()
    outfile.close()