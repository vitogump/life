'''
Created on 2015-12-27

@author: liurui
'''
from optparse import OptionParser
import os
import re
from scipy import stats
# import rpy2.robjects as robjects

parser = OptionParser()
parser.add_option("-g", "--gotablefile", dest="gotablefile", help="minvalue maxvalue breaks", metavar="FILE")
parser.add_option("-t", "--type", dest="type",help="g:go k:kegg")
parser.add_option("-o","--outpre",dest="outpre",help="outfilepreName with path")
parser.add_option("-E","--EntrezGene_ID",dest="EntrezGene_ID",help="outfilepreName with path")
parser.add_option("-U","--UniProtlistfile",dest="UniProtlistfile",help="file depthfile")#
parser.add_option("-T","--tplistfile",dest="tplist",help="file depthfile")#
parser.add_option("-G","--genelistfile",dest="genelistfile",help="file depthfile")#                                                                                                                                                      
(options, args) = parser.parse_args()
def keggenrichmentForHuman(outpre,unipIDlistfile):
    curpath=os.getcwd()
    r = robjects.r
    r("setwd('" + curpath + "')")
    r('.libPaths("/opt/Rpackages/")')
    r("library('org.Hs.eg.db')")
    r("library('GSEABase')")
    r("library('GOstats')")
    r("library('KEGG.db')")
    r("library('pathview')")
    r("library('Category')")
    r('genes<-read.table("'+options.EntrezGene_ID+'",header=T)')
    r('entrezID <- as.character(genes$EntrezGene_ID)')
    r('keggAnn <- get("org.Hs.egPATH")')
    r('universe <- Lkeys(keggAnn)')
    r("""params <- new("KEGGHyperGParams", 
                    geneIds=entrezID, 
                    universeGeneIds=universe, 
                    annotation="org.Hs.eg.db", 
                    categoryName="KEGG", 
                    pvalueCutoff=0.05,
                    testDirection="over")""")
    r("over <- hyperGTest(params)")
    r("kegg <- summary(over)")
    r("glist <- geneIdsByCategory(over)")
    r("""glist <- sapply(glist, function(.ids) {
     .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
     .sym[is.na(.sym)] <- .ids[is.na(.sym)]
     paste(.sym, collapse=";")
     })""")
    r("glist <- geneIdsByCategory(over)")
    r("""glist <- sapply(glist, function(.ids) {
     .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
     .sym[is.na(.sym)] <- .ids[is.na(.sym)]
     paste(.sym, collapse=";")
     })""")
    r("kegg$Symbols <- glist[as.character(kegg$KEGGID)]")
    r('write.table(kegg,file="'+outpre.strip()+'pathway_enrichment.txt")')
def GOenrichment(gotablefile,outpre,genelist=None,trscptlist=None,UniProtlist=None):
    gotablefile=open(gotablefile,'r')
    title = gotablefile.readline()
    titlelist= [e.strip().lower() for e in re.split(r"\t",title)]
    UniProtidx=titlelist.index("uniprot/trembl accession")
    geneididx=titlelist.index("ensembl gene id")
    tpididx=titlelist.index("ensembl transcript id")
    gotermaccessionidx=titlelist.index("go term accession")
    gotermNameidx=titlelist.index("go term name")
    godomainidx=titlelist.index("go domain")
    goternDefinition=titlelist.index("go term definition")
    genenameidx=titlelist.index("associated gene name")    
    
    if genelist!=None:
        sampledIDlist=genelist
#         ensemblIDlistfile=open(genelist,"r")
        IDidx=geneididx
    elif trscptlist!=None:
        sampledIDlist=trscptlist
#         ensemblIDlistfile=open(trscptlist,'r')
        IDidx=tpididx
    elif UniProtlist!=None:
        sampledIDlist=UniProtlist
        IDidx=UniProtidx
    
#     sampledIDlist=ensemblIDlistfile.readlines()
    print(sampledIDlist,sep="\n")
    
#     ensemblIDlistfile.close()
    

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
#             print(geneName,file=open("test.txt",'a'))
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
    GOAnnationForGene_out_fileName=outpre.strip()+".GO_annotion"
    GOenrichment_fileName=outpre.strip()+".GO_enrichment"
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
if __name__ == '__main__':
    uniontpidlist=[]
    if options.type =="g":
        outpre=options.outpre
        if options.UniProtlistfile is not None:
            UniProtlistfile=options.UniProtlistfile
            f=open(UniProtlistfile,"r")
            for line in f:
                uniontpidlist.append(line.strip())
            f.close()
            GOenrichment(options.gotablefile,options.outpre,None,None,list(set(uniontpidlist)))
        elif options.tplist is not None:
            tplistfile=options.tplist
            f=open(tplistfile,"r")
            for line in f:
                uniontpidlist.append(line.strip())
            f.close()
            GOenrichment(options.gotablefile,options.outpre,None,list(set(uniontpidlist)),None)
        else:
            genelistfile=options.genelistfile
            f=open(genelistfile,"r")
            for line in f:
                uniontpidlist.append(line.strip())
            f.close()
            GOenrichment(options.gotablefile,options.outpre,list(set(uniontpidlist)),None,None)