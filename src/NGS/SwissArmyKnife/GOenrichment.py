from optparse import OptionParser
import os
import re, sys, time
from scipy import stats
from scipy.stats import hypergeom
from NGS.BasicUtil import *


'''
Created on 2015-4-10

@author: liurui
'''
parser = OptionParser()
parser.add_option("-g", "--gotablefile", dest="gotablefile", help="gotable title with :Ensembl Gene ID    Ensembl Transcript ID    GO Term Accession    GO Term Evidence Code    GO domain    GO Term Name    GO Term Definition,order and upper/lower case is arbitrarily")
parser.add_option("-G","--genelist",dest="genelist",default=None)
parser.add_option("-T","--trscptlist",dest="trscptlist",default=None)
parser.add_option("-o","--outpre",dest="outpre")

(options, args) = parser.parse_args()
if __name__ == '__main__':
    gotablefile=open(options.gotablefile,'r')
    title = gotablefile.readline()
    titlelist= [e.strip().lower() for e in re.split(r"\t",title)]
    geneididx=titlelist.index("ensembl gene id")
    tpididx=titlelist.index("ensembl transcript id")
    gotermaccessionidx=titlelist.index("go term accession")
    gotermNameidx=titlelist.index("go term name")
    godomainidx=titlelist.index("go domain")
    goternDefinition=titlelist.index("go term definition")
    genenameidx=titlelist.index("associated gene name")    
    
    if options.genelist!=None:
        ensemblIDlistfile=open(options.genelist,"r")
        IDidx=geneididx
    elif options.trscptlist!=None:
        ensemblIDlistfile=open(options.trscptlist,'r')
        IDidx=tpididx
    
    sampledIDlist=ensemblIDlistfile.readlines()
    print(sampledIDlist,sep="\n")
    
    ensemblIDlistfile.close()
    

    gotable={}
    """
    gotable={tp_id1:(geneID,geneName,bp,cc,mf),tp_id2:(geneID,geneName,bp,cc,mf),,,,,,}
    """
    oneGO2manyID={}
    """
    oneGO2manyID={go_Accession1:[tp1,tp2,...],go_Accession2:[],....}
    """
    goTermMap={}
    """
    goTermMap={go_Accession1:[go term name,go domain],go_Accession2:[],,,,}
    """
    
    bp="";cc="";mf="";geneName="";geneID=""
    genelist=[]
    for termline in  gotablefile:
        termlist=re.split(r"\t",termline)
        if termlist[gotermaccessionidx].strip() in oneGO2manyID:
            goTermMap[termlist[gotermaccessionidx].strip()]+=[termlist[gotermNameidx],termlist[godomainidx]]
            oneGO2manyID[termlist[gotermaccessionidx].strip()].append(termlist[IDidx].strip())
        else:
            goTermMap[termlist[gotermaccessionidx].strip()]=[termlist[gotermNameidx],termlist[godomainidx]]
            oneGO2manyID[termlist[gotermaccessionidx].strip()]=[termlist[IDidx].strip()]
        if termlist[geneididx].strip() not in genelist:
            genelist.append(termlist[geneididx].strip())
        if termlist[tpididx].strip() in gotable:
            if termlist[godomainidx].lower().strip()=="biological_process":
                bp+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
            elif termlist[godomainidx].lower().strip()=="cellular_component":
                cc+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"                 
            elif termlist[godomainidx].lower().strip()=="molecular_function":
                mf+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
           
        else:#new gene start

            geneID=termlist[geneididx]
            geneName=termlist[genenameidx]
            print(geneName,file=open("test.txt",'a'))
            bp="";cc="";mf=""
            if termlist[godomainidx].lower().strip()=="biological_process":
                bp+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
            elif termlist[godomainidx].lower().strip()=="cellular_component":
                cc+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"            
            elif termlist[godomainidx].lower().strip()=="molecular_function":
                mf+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
            if geneName.split():
                gotable[termlist[IDidx].strip()]=(geneID,geneName,bp,cc,mf)
            else:
                gotable[termlist[IDidx].strip()]=(geneID,"unknow",bp,cc,mf)            
    GOAnnationForGene_out_fileName=options.outpre.strip()+".GO_annotion"
    GOenrichment_fileName=options.outpre.strip()+".GO_enrichment"
    annf=open(GOAnnationForGene_out_fileName,'w')
    enrichfile=open(GOenrichment_fileName,'w')
    all_IDlist=list(gotable.keys());m_n=len(genelist);del genelist
    for id in sampledIDlist:
        id=id.strip()
        if id not in gotable:
            print(id,"don't have go annotion")
            continue
        print(id,gotable[id][0].strip(),gotable[id][1].strip(),gotable[id][2].strip(),gotable[id][3].strip(),gotable[id][4].strip(),sep="\t",file=annf)
    outlist=[]
    """
    outlist=[(go_Accession1,go_term_name,go_domain,p-value,FDR,sampled_inTerm,termsize),(),,,,]
    """
    testGOONETOMANY=open("GOONETOmany.txt",'w')
    for goID in oneGO2manyID.keys():
        print(goID,*oneGO2manyID[goID],sep="\t",file=testGOONETOMANY)
    testGOONETOMANY.close()
    k=len(sampledIDlist)
    for goassecesion in sorted(oneGO2manyID.keys()):
        x=0
        containingtrscript=[]
        genetermlist=[]
        for id in sampledIDlist:
            id=id.strip()
            if id in oneGO2manyID[goassecesion]:
                containingtrscript.append(id)
                genetermlist.append(gotable[id][1])
                x+=1
        m=len(oneGO2manyID[goassecesion])
        n=m_n - m
        pvalue=stats.hypergeom.sf(x-1,m_n,m,k)
        if len(goTermMap[goassecesion])<2:
            continue
        outlist.append((goassecesion,goTermMap[goassecesion][0],goTermMap[goassecesion][1],pvalue,"FDR",x,len(oneGO2manyID[goassecesion]),containingtrscript,genetermlist))
    outlist.sort(key=lambda listRec:listRec[3])
    for e in outlist:
        print(*e,sep="\t",file=enrichfile)
    enrichfile.close()
    os.system("""awk 'BEGIN{FS="\t"}$3~/biological_process/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_biological_process")
    os.system("""awk 'BEGIN{FS="\t"}$3~/cellular_component/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_cellular_component")
    os.system("""awk 'BEGIN{FS="\t"}$3~/molecular_function/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_molecular_function")
    annf.close()
    gotablefile.close()
    print("finish")