'''
Created on 2020年9月24日

@author: RuiLiu
'''
from optparse import OptionParser
import re
import pandas as pd
parser = OptionParser()

parser.add_option("-K", "--KEGGAnnotation", dest="KOAnnot", help="gotable title with :GeneID    KO    GeneName    Definition    GO domain    GO Term Name    GO Term Definition,order and upper/lower case is arbitrarily")
parser.add_option("-G","--geneabundance",dest="geneabund",default=None)

parser.add_option("-o","--outpre",dest="outpre")

(options, args) = parser.parse_args()
"""
awk '{print $1}' KEGG_Annotation.txt|sort|uniq -d|less -S 发现没有重复！！！！！
"""
outputfile=open(options.outpre,'w')
geneKeyKOvalue={}
"""
geneID:KO
"""
KOcount={};sample1idx=4
"""
KO:[G1SUM,G2SUM,,,,,S1SUM,S2SUM,,,]
"""
if __name__ == '__main__':
    
    abundata=pd.read_table(options.geneabund,header=0)#open(options.geneabund,'r')
    titlelist=list(abundata.columns)#re.split(r"\t",abundfile.readline().strip())
    KOAnnofile=open(options.KOAnnot,'r')
    KOAnnofile.readline()
    ValueNumb=len(titlelist)-sample1idx
    print(ValueNumb)
    # store KOAnnofile into geneKeyKOvalue structure
    for line in KOAnnofile:
        linelist=re.split(r"\t",line.strip())
        if linelist[0] in geneKeyKOvalue:
            print("error,something wrong")
        else:
            geneKeyKOvalue[linelist[0]]=linelist[1];KOcount[linelist[1]]=[0]*ValueNumb
    KOAnnofile.close()
    #均一化
    abundatafrac=pd.DataFrame(columns=[abundata.columns[0]]+abundata.columns.tolist()[sample1idx:])#[geneID,G1,G2,,,,S1,S2,,,,] order as the same
    abundatafrac[abundata.columns[0]]=abundata[abundata.columns[0]]
    for id in titlelist[4:]:
        abundatafrac[id]=abundata[id]/abundata[id].sum()
    
    #print out
    print("gut","\t".join(titlelist[sample1idx:]),sep="\t",file=outputfile)
    i=0
    #Sum in KO && place KO|geneID to print
    for idx,gabu in abundatafrac.iterrows():
        gabul=gabu.tolist()
        if gabul[0] in geneKeyKOvalue:# abundence recod gene has a KO 
            for samp_idx in range(len(gabul[1:])):
                KOcount[geneKeyKOvalue[gabul[0]]][samp_idx]+=int(gabul[1+samp_idx])
            print(geneKeyKOvalue[gabul[0]]+"|"+gabul[0],*gabul[1:],sep="\t",file=outputfile)
    else:
        for KO in KOcount.keys():
            print(KO,*KOcount[KO],sep="\t",file=outputfile)
            

    outputfile.close()