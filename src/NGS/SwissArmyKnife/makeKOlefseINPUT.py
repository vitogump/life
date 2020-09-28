'''
Created on 2020年9月24日

@author: RuiLiu
'''
from optparse import OptionParser
import re
from scipy import stats
import pandas as pd
parser = OptionParser()

parser.add_option("-K", "--KEGGAnnotation", dest="KOAnnot", help="gotable title with :GeneID    KO    GeneName    Definition    GO domain    GO Term Name    GO Term Definition,order and upper/lower case is arbitrarily")
parser.add_option("-P","--KEGGPathway",dest="Kpathwayfile")
parser.add_option("-G","--geneabundance",dest="geneabund",default=None)

parser.add_option("-o","--outpre",dest="outpre")

(options, args) = parser.parse_args()
"""
awk '{print $1}' KEGG_Annotation.txt|sort|uniq -d|less -S 发现没有重复！！！！！
"""
outputfile=open(options.outpre,'w');koout=open(options.outpre+".kopathway",'w')
geneKeyKOvalue={};geneKeykova_rev={}
"""
geneID:KO
reverse pathwayfile => geneID:[ko1,ko2,,,]
"""
KOcount={};kocount={}
sample1idx=4
"""
KO1:[G1SUM,G2SUM,,,,,S1SUM,S2SUM,,,]
"""
KOKeyGeneidxvalue={}
"""
KO1:[g1idx1,g1idx2,,],KO2:[],,,
"""
        
        
if __name__ == '__main__':
    
    abundata=pd.read_table(options.geneabund,header=0)#open(options.geneabund,'r')
    titlelist=list(abundata.columns)#re.split(r"\t",abundfile.readline().strip())
    KOAnnofile=open(options.KOAnnot,'r')
    KOAnnofile.readline()
    ValueNumb=len(titlelist)-sample1idx
    print(ValueNumb)
    # read pathway file
    pathwfile=open(options.Kpathwayfile,'r')
    pathwfile.readline()
    for line in pathwfile:
        linelist=re.split(r'\t',line.strip().strip("\""))
        if re.search(r'ko\d+\s[\w\s]+',linelist[0]):
            genelist=re.split(r";",linelist[2])
            if len(linelist)>3:
                KOlist=re.split(r";",linelist[3])
            else:
                print(linelist)
            for geneid in genelist:
                if geneid in geneKeykova_rev:
                    geneKeykova_rev[geneid].append(linelist[0])
                else:
                    geneKeykova_rev[geneid]=[linelist[0]]
            kocount[linelist[0]]=[0]*ValueNumb    
    
    # store KOAnnofile into geneKeyKOvalue structure
    for line in KOAnnofile:
        linelist=re.split(r"\t",line.strip())
        if linelist[0] in geneKeyKOvalue:
            print("error,something wrong")
        else:
            geneKeyKOvalue[linelist[0]]=linelist[1]
            KOKeyGeneidxvalue[linelist[1]]=[];KOcount[linelist[1]]=[0]*ValueNumb
    KOAnnofile.close()
    #均一化
    abundatafrac=pd.DataFrame(columns=[abundata.columns[0]]+abundata.columns.tolist()[sample1idx:])#[geneID,G1,G2,,,,S1,S2,,,,] order as the same
    abundatafracp=pd.DataFrame(columns=["KOid|geneid"]+abundata.columns.tolist()[sample1idx:]+["pvalue"])# gene passed kruskal test
    abundatafrac[abundata.columns[0]]=abundata[abundata.columns[0]]
    abundatafracp["KOid|geneid"]=abundata[abundata.columns[0]]
    for id in titlelist[4:]:
        abundatafrac[id]=abundata[id]/abundata[id].sum()
    
    #print out:geneid|KOid    G1    G2    G3    G4    G5    S1    S2    S3    S4    S5
    print("KOid|geneid","\t".join(titlelist[sample1idx:]),"pvalue",sep="\t",file=outputfile)
    i=0
    #Sum in KO && place KO|geneID to print
    for idx,gaburow in abundatafrac.iterrows():
        gabul=gaburow.tolist()
        #Kruskal test for each gene
        try:
            statistc,pvalue=stats.kruskal(gabul[1:6],gabul[6:])
        except:
            pvalue=1
        if gabul[0] in geneKeyKOvalue:# abundence recod gene has a KO 
            KO=geneKeyKOvalue[gabul[0]]
            KOKeyGeneidxvalue[KO].append(idx)
            for samp_idx in range(len(gabul[1:])):#every sample for one gene rec
                KOcount[KO][samp_idx]+=gabul[1+samp_idx]
            print(KO+"|"+gabul[0],*gabul[1:],pvalue,sep="\t",file=outputfile)
#             print(pvalue<=0.05,len(abundatafracp))
            if pvalue<=0.05:
                abundatafracp.loc[len(abundatafracp)]=[KO+"|"+gabul[0],*gabul[1:],pvalue]
#                 print([KO+"|"+gabul[0],*gabul[1:],pvalue])
#                 abundatafracp.append(pd.DataFrame([KO+"|"+gabul[0],*gabul[1:],pvalue]),ignore_index=True)
        if gabul[0] in geneKeykova_rev:
            for ko in geneKeykova_rev[gabul[0]]:
                for samp_idx in range(len(gabul[1:])):
                    kocount[ko][samp_idx]+=gabul[1+samp_idx]
    else:
        for KO in KOcount.keys():
            statistic,pvalue=stats.kruskal(KOcount[KO][:5],KOcount[KO][5:])
            print(KO,*KOcount[KO],pvalue,sep="\t",file=outputfile)
        for ko in kocount.keys():
            statistic,pvalue=stats.kruskal(kocount[ko][:5],kocount[ko][5:])
            print(ko,*kocount[ko],pvalue,sep="\t",file=koout)
            
    #colect gene that the KO passed Kruskal test 
#     for KO in KOcount.keys():
#         statistic,pvalue=stats.kruskal(KOcount[KO][:5],KOcount[KO][5:])
#         if pvalue<=0.05:
#             for geneidx in KOKeyGeneidxvalue[KO]:
#                 print(abundatafrac[geneidx],)
#             print()
    #wilcoxon test 
    print("no subclass for sample, no wilcoxon test need")
    #lda test
    outputfile.close();koout.close()
    print(abundatafracp)
    #read to pandas dataframe and select gene that passed Kruskal test <0.05
    