import copy
import os,numpy
import re,pickle
import time,math,config

from Bio import SeqIO
from scipy import stats


import src.NGS.BasicUtil.DBManager as dbm


SLEEP_FOR_NEXT_TRY=3
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
    print("GOenrichmentinput:",sampledIDlist,sep="\n")
    
#     ensemblIDlistfile.close()
    

    gotable={}
    """
    gotable={tp_id1:(geneID,geneName,bp,cc,mf,description),tp_id2:(geneID,geneName,bp,cc,mf,description),,,,,,}
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
                gotable[termlist[IDidx].strip()]=(geneID,geneName,bp,cc,mf,termlist[14])
            else:
                gotable[termlist[IDidx].strip()]=(geneID,"unknow",bp,cc,mf,termlist[14])     
    else:
        print()       
    GOAnnationForGene_out_fileName=outpre.strip()+".GO_annotion"
    GOenrichment_fileName=outpre.strip()+".GO_enrichment"
    annf=open(GOAnnationForGene_out_fileName,'w')
    print("ensembl trscptID","ensembl geneID","gene symbol","go number","description",file=annf)
    enrichfile=open(GOenrichment_fileName,'w')
    all_IDlist=list(gotable.keys());m_n=len(genelist);del genelist
    for id in sampledIDlist:
        id=id.strip()
        if id not in gotable:
            print(id,"don't have go annotion",file=annf)
            continue
        print(id,gotable[id][0].strip(),gotable[id][1].strip(),gotable[id][2].strip(),gotable[id][3].strip(),gotable[id][4].strip(),gotable[id][5].strip(),sep="\t",file=annf)
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
        if len(goTermMap[goassecesion])<2:
            continue        
        for id in sampledIDlist:
            id=id.strip()
            if id in oneGO2manyID[goassecesion]:
                containingtrscript.append(id)
                genetermlist.append(gotable[id][1])
                x+=1
        m=len(oneGO2manyID[goassecesion])
        n=m_n - m
        pvalue=stats.hypergeom.sf(x-1,m_n,m,k)

        outlist.append((goassecesion,goTermMap[goassecesion][0],goTermMap[goassecesion][1],pvalue,"FDR",x,len(oneGO2manyID[goassecesion]),containingtrscript,genetermlist))
    outlist.sort(key=lambda listRec:listRec[3]);NoOftestingBP=0
    for e in outlist:
        if e[3]<1 and e[2]=="biological_process":
            NoOftestingBP+=1
        print(*e,sep="\t",file=enrichfile)
    enrichfile.close()
    os.system("""awk 'BEGIN{FS="\t"}$3~/biological_process/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_biological_process")
#     os.system("""awk 'BEGIN{FS="\t"}$3~/cellular_component/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_cellular_component")
#     os.system("""awk 'BEGIN{FS="\t"}$3~/molecular_function/{print $0}' """+GOenrichment_fileName+">"+GOenrichment_fileName+"_molecular_function")
    #multiple test control of FDR for p<0.05
    GObpwithP=[]
    enrichGObpf=open(GOenrichment_fileName+"_biological_process",'r')
    GObpwithP=enrichGObpf.readlines()
    i=1;print(NoOftestingBP)
    for testTerm in GObpwithP:
        tTliest=re.split(r'\t',testTerm.strip())
        FDRv=(i/NoOftestingBP)*0.05
        if float(tTliest[3]) <= FDRv:
            GObpwithP[i-1]=GObpwithP[i-1].replace('\tFDR\t',"\t"+str(FDRv)+"\tTRUE\t")
        else:
            GObpwithP[i-1]=GObpwithP[i-1].replace('\tFDR\t',"\t"+str(FDRv)+"\tFALSE\t")
        i+=1
    enrichGObpf.close()
    OUTagainGObp=open(GOenrichment_fileName+"_biological_process",'w')
    for testTerm in GObpwithP:
        print(testTerm.strip(),file=OUTagainGObp)
    OUTagainGObp.close()
    annf.close()
    gotablefile.close()
def collectSNP_locatInRegion(MultipleVcfMap,chrom,startpos,endpos):
    """
    return [[pos1,ref,alt,(),(),()],[POS2],[],,,]
    """
    low=0;high=len(MultipleVcfMap[chrom])-1
    if MultipleVcfMap[chrom]==[] or startpos>=MultipleVcfMap[chrom][-1][0] or endpos<=MultipleVcfMap[chrom][0][0]:
        print("collectSNP_locatInRegion",chrom,startpos,endpos,"empty")
        return []
    while low<=high:
        mid=(low + high)>>1
        if MultipleVcfMap[chrom][mid][0]<endpos:
            low=mid+1
        elif MultipleVcfMap[chrom][mid][0]>endpos:

            high=mid-1
        else:
            end_idx=mid
            break
    else:
        if low>=len(MultipleVcfMap[chrom]):
            
            low=len(MultipleVcfMap[chrom])-1
        if MultipleVcfMap[chrom][low][0]>endpos and low>0:
            low-=1
        print(MultipleVcfMap[chrom][low][0])
        end_idx=low
    low=0;high=len(MultipleVcfMap[chrom])-1
    while low<=high:
        mid=(low + high)>>1
        if MultipleVcfMap[chrom][mid][0]<startpos:
            low=mid+1
        elif MultipleVcfMap[chrom][mid][0]>startpos:
            high=mid-1
        else:
            start_idx=mid
            break
    else:
        start_idx=high
        if MultipleVcfMap[chrom][high][0]<startpos and len(MultipleVcfMap[chrom])>high+1:
            start_idx+=1
    
    if MultipleVcfMap[chrom][start_idx][0]>=startpos and MultipleVcfMap[chrom][start_idx][0]<=endpos and len(MultipleVcfMap[chrom])>(end_idx+1) and MultipleVcfMap[chrom][end_idx+1][0]>=startpos and MultipleVcfMap[chrom][end_idx+1][0]<=endpos:
        return copy.deepcopy(MultipleVcfMap[chrom][start_idx:end_idx+1])
    elif MultipleVcfMap[chrom][end_idx][0]>=startpos and MultipleVcfMap[chrom][end_idx][0]<=endpos:
        return copy.deepcopy(MultipleVcfMap[chrom][start_idx:end_idx])
    else:
        return []
def mapWinvaluefileToChrOfReletiveSpecie(anchorfile,winfileinName,winwidth,slidesize,standardsexseperately=True,mapfile=None):
    newanchorfilehandler=anchorfile
    if mapfile:
        scaffoldmap={}
        mapfile=open(mapfile,'r')
        for line in mapfile:
            linelist=re.split(r"\s+",line.strip())
            scaffoldmap[linelist[0].strip().lower()]=linelist[1].strip()
        oldanchorfilehandler=open(anchorfile,'r')
        newanchorfilehandler=open(anchorfile+"changed",'w')
        for line in oldanchorfilehandler:
            linelist=re.split(r"\s+",line.strip())
            if linelist[3] in scaffoldmap:
                linelist[3]=scaffoldmap[linelist[3].strip().lower()]
            print(*linelist,sep="\t",file=newanchorfilehandler)
        oldanchorfilehandler.close()
        newanchorfilehandler.close()
            
        
    anchorDATASTRUCTURE={}
    """
    {chr1:[(53353,53806,scaffold451,558997,558537,-),(57200,62371,scaffold451,553669,548504,-),(),,],chr2:[],,,,}
    """
    reverseAnchorDATASTRUCTURE={}
    """
    {scaffold451:{chr1:[0,1,2,,,,]},C17734302:{chr1:[idx]}}  idx is idx in the list of anchorDATASTRUCTURE[chr1] 
    """
    if mapfile:
        
        newanchorfilehandler=open(anchorfile+"changed",'r')
    else:
        newanchorfilehandler=open(anchorfile,'r')
    for line in newanchorfilehandler:
        linelist=re.split(r"\s+",line.strip())
        if linelist[0].strip() in anchorDATASTRUCTURE:
            anchorDATASTRUCTURE[linelist[0].strip()].append((int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip()))
        else:
            anchorDATASTRUCTURE[linelist[0].strip()]=[(int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip())]
        #fill reverseAnchorDATASTRUCTURE
        if linelist[3].strip() in reverseAnchorDATASTRUCTURE:
            if linelist[0].strip() in reverseAnchorDATASTRUCTURE[linelist[3].strip()]:
                reverseAnchorDATASTRUCTURE[linelist[3].strip()][linelist[0].strip()].append(len(anchorDATASTRUCTURE[linelist[0].strip()])-1)
            else:
                reverseAnchorDATASTRUCTURE[linelist[3].strip()][linelist[0].strip()]=[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]
        else:
            reverseAnchorDATASTRUCTURE[linelist[3].strip()]={linelist[0].strip():[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]}
    newanchorfilehandler.close()
#     if __name__ == '__main__':
    winfile=open(winfileinName,'r')
    title=winfile.readline()
    winMap={}#{scaffold:[(startpos,endpos,noofsnp,winvalue,zvalue),(),(),,,]}
    for line in winfile:
        linelist=re.split(r"\s+",line.strip())
        if linelist[0].strip()  in winMap:
            winMap[linelist[0].strip()].append((int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6]))
        else:
            winMap[linelist[0].strip()]=[(int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6])]
    winfile.close()
    #winfile has been loaded into memonery
    #print a new winfile mark auto or sex chromosome
    winMapMarked=copy.deepcopy(winMap)

    outwinfile=open(winfileinName+"arrangemented",'w')
    print(title.strip()+"\tmark",file=outwinfile)
    for chrom in sorted(anchorDATASTRUCTURE.keys()):
        idx=0#it seems not useful
        if anchorDATASTRUCTURE[chrom][idx][5]=="-":
            regionstart=anchorDATASTRUCTURE[chrom][idx][4]
            regionend=anchorDATASTRUCTURE[chrom][idx][3]
        elif anchorDATASTRUCTURE[chrom][idx][5]=="+":
            regionstart=anchorDATASTRUCTURE[chrom][idx][3]
            regionend=anchorDATASTRUCTURE[chrom][idx][4]
        lastscaffold=anchorDATASTRUCTURE[chrom][idx][2]
        if lastscaffold not in winMap:
            for startpos,endpos,scaffold,sstartpos,sendpos,foward_reverse in anchorDATASTRUCTURE[chrom]:
                idx+=1
                if scaffold in winMap:
                    lastscaffold=scaffold
                    break

        for startpos,endpos,scaffold,sstartpos,sendpos,foward_reverse in anchorDATASTRUCTURE[chrom][idx:]:
            if lastscaffold==scaffold:
#                 print("this code block counting the continue stretch . and Should consider wither there are some gap in it")
                if foward_reverse=="-":
                    regionstart=min(regionstart,sendpos)
                    regionend=max(regionend,sstartpos)
                else:
                    regionstart=min(regionstart,sstartpos)
                    regionend=max(regionend,sendpos)
                
            else:
                print("after determining the start or the end of the stretch, arranging the scaffolds which in the stretch")
                if regionstart < int(winwidth):
                    winstartNo=0
                else:
                    winstartNo=math.ceil((regionstart-int(winwidth))/int(slidesize))
                if regionend<int(winwidth):
                    winendNo=0
                else:
                    winendNo=math.ceil((regionend-int(winwidth))/int(slidesize))
                if anchorDATASTRUCTURE[chrom][idx-1][5]=="-":#last scaffold so here is idx-1
                    for i in range(winstartNo,winendNo+1)[::-1]:
#                         if lastscaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[lastscaffold] and  "z" not in reverseAnchorDATASTRUCTURE[lastscaffold] and  "W" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "w" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "X" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "x" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "Y" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "y" not in reverseAnchorDATASTRUCTURE[lastscaffold]):
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+[chrom])
#                         else:
#                             winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                        print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
                else:
                    for i in range(winstartNo,winendNo+1):
#                         if lastscaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[lastscaffold] and  "z" not in reverseAnchorDATASTRUCTURE[lastscaffold] and  "W" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "w" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "X" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "x" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "Y" not in reverseAnchorDATASTRUCTURE[lastscaffold] and "y" not in reverseAnchorDATASTRUCTURE[lastscaffold]):
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+[chrom])
#                         else:
#                             winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                        print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
                print("new region")
                if foward_reverse=="-":
                    regionstart=sendpos
                    regionend=sstartpos
                elif foward_reverse=="+":
                    regionstart=sstartpos
                    regionend=sendpos
                if scaffold not in winMap:
                    continue
                lastscaffold=scaffold
            idx+=1
        else:
            if regionstart < int(winwidth):
                winstartNo=0
            else:
                winstartNo=math.ceil((regionstart-int(winwidth))/int(slidesize))
            if regionend<int(winwidth):
                winendNo=0
            else:
                winendNo=math.ceil((regionend-int(winwidth))/int(slidesize))
            if anchorDATASTRUCTURE[chrom][idx-1][5]=="-":
                for i in range(winstartNo,winendNo+1)[::-1]:
#                     if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                    winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+[chrom])
#                     else:
#                         winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                    print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
            else:
                for i in range(winstartNo,winendNo+1):
#                     if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                    winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+[chrom])
#                     else:
#                         winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                    print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
            
#             for i in range(winstartNo,winendNo+1):
#                 print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],scaffold,sep="\t",file=outwinfile)
    outwinfile.close()
    winCrossGenomeMap={"autosome":[],"sexchromosome":[]}#{"autosome":[],"Z":[],"W":[],"X":[],"Y":[]}
    outwinfile=open(winfileinName+"marked",'w')
    print(title.strip()+"\tmark",file=outwinfile)
    for scaffold in sorted(winMapMarked.keys()):
        for winNo in range(len(winMapMarked[scaffold])):

            if len(winMapMarked[scaffold][winNo])==5 :
                if (scaffold in reverseAnchorDATASTRUCTURE) and ((len(reverseAnchorDATASTRUCTURE[scaffold])==1 and  re.search(r"[zwxy]" , "".join(reverseAnchorDATASTRUCTURE[scaffold].keys()).lower())!=None) or (len(reverseAnchorDATASTRUCTURE[scaffold])>1 and re.search(r"([zwxyZWXY]+)" , "".join(reverseAnchorDATASTRUCTURE[scaffold].keys()))!=None and len(reverseAnchorDATASTRUCTURE[scaffold][re.search(r"([zwxyZWXY]+)" , "".join(reverseAnchorDATASTRUCTURE[scaffold].keys())).group(1)[-1]])>3)):
                    print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"sexchromosome",sep="\t",file=outwinfile)#,"unknow"
                    signal="sexchromosome"
                else:
                    print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"autosome",sep="\t",file=outwinfile)#,"unknow"
                    signal="autosome"
            elif winMapMarked[scaffold][winNo][5].upper()=="Z" or winMapMarked[scaffold][winNo][5].lower()=="w" or winMapMarked[scaffold][winNo][5].upper()=="X" or winMapMarked[scaffold][winNo][5].lower()=="y":
                print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"sexchromosome",sep="\t",file=outwinfile)#winMapMarked[scaffold][winNo][5],
                signal="sexchromosome"
            else:
                print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"autosome",sep="\t",file=outwinfile)#winMapMarked[scaffold][winNo][5],
                signal="autosome"
            if  re.search(r"^[1234567890\.e-]+$",winMapMarked[scaffold][winNo][3])!=None:
                winCrossGenomeMap[signal].append(float(winMapMarked[scaffold][winNo][3]))
    autoexception=numpy.mean(winCrossGenomeMap["autosome"])
    autostd1=numpy.std(winCrossGenomeMap["autosome"],ddof=1)
    sexexception=numpy.mean(winCrossGenomeMap["sexchromosome"])
    sexstd1=numpy.std(winCrossGenomeMap["sexchromosome"],ddof=1)
    print("autoexception,autostd",autoexception,autostd1,"sexchromosome:",sexexception,sexstd1)
    outwinfile.close()
    if standardsexseperately:
        markfilname=winfileinName+"marked.sexchromseperatestandard"
        markedfile=open(winfileinName+"marked","r")
        markedseperatelyfile=open(markfilname,'w')
        title=markedfile.readline()
        print(title,end="",file=markedseperatelyfile)
        for line in markedfile:
            linelist=re.split(r"\s+",line.strip())
            if re.search(r"^[1234567890\.e-]+$",linelist[5])!=None:
                if linelist[7] =="sexchromosome":
                    zscore=(float(linelist[5])-sexexception)/sexstd1
                elif linelist[7] =="autosome":
                    zscore=(float(linelist[5])-autoexception)/autostd1
                else:
                    print("what's wrong");exit(-1)
                print(linelist[0],linelist[1],linelist[2],linelist[3],linelist[4],linelist[5],zscore,linelist[7],sep="\t",file=markedseperatelyfile)
            else:
                print(line,end="",file=markedseperatelyfile)
        markedseperatelyfile.close();markedfile.close()
    else:
        markfilname=winfileinName+"marked.sexchromseperatestandard"
    
    return markfilname,winfileinName+"arrangemented"
class WinInGenome():           
    def __init__(self, dbname, winFileName8Field,Nocol=7, tableName=None):
        super().__init__()
        self.dbname = dbname
        self.chromOrder, self.windbtools, self.wintablewithoutNA, self.wintabletextvalueallwin = self.loadWinDataIntoDB(dbname, winFileName8Field,Nocol, tableName)
        self.winContainTrscptMap = {}
    def loadWinDataIntoDB(self, dbname, winFileName8Field,Nocol="7", tableNamewithoutNA=None):
        chromOrder = []
        
        tempdbtools = dbm.DBTools(config.ip, config.username, config.password, config.dbname)
        if tableNamewithoutNA == None:
            tableNamewithoutNA = config.random_str()
#             return chromOrder, tempdbtools, tableNamewithoutNA, tableNametextValueForappendGeneName 
        tableNametextValueForappendGeneName = tableNamewithoutNA + "textField"
        
        TABLES = {}
        TABLES[tableNamewithoutNA] = (
            "CREATE TABLE " + tableNamewithoutNA + " ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` varchar(128) NOT NULL,"
            " `bp_start` varchar(128) NOT NULL,"
            " `bp_end` varchar(128) NOT NULL,"
            " `snpcount` int(11) NOT NULL,"
            " `winvalue` double NOT NULL,"  #########why?
            " `zvalue` double NOT NULL,"  ##########
            " `mark` varchar(30) NOT NULL DEFAULT 'unknown', "
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")"
            )
        TABLES[tableNametextValueForappendGeneName] = (
            "CREATE TABLE " + tableNametextValueForappendGeneName + " ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` varchar(128) NOT NULL,"
            " `bp_start` varchar(128) NOT NULL,"
            " `bp_end` varchar(128) NOT NULL,"
            " `snpcount` int(11) NOT NULL,"
            " `winvalue` text NOT NULL,"  ##############why?
            " `zvalue` text NOT NULL,"  ###############
            " `mark` varchar(30) NOT NULL DEFAULT 'unknown', "
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")"
            )        
        print(TABLES)
        tempdbtools.create_table(TABLES)
        a = os.popen("awk '{print $1}' " + winFileName8Field + "|uniq")
        for chromNo in a:
            chromOrder.append(chromNo.strip())
        a.close()
        a = os.system("awk '$"+str(Nocol)+"!~/NA/ && NR!=1{print $0}' " + winFileName8Field + ">" + winFileName8Field + "_tmpfile")
        if a != 0:
            print("awk '$"+str(Nocol)+"!~/NA/ && NR!=1{print $0}' " + winFileName8Field + ">" + winFileName8Field + "_tmpfile" + ": failed")
            exit(-1)
        print("awk '$"+str(Nocol)+"!~/NA/ && NR!=1{print $0}' " + winFileName8Field + ">" + winFileName8Field + "_tmpfile" + ": ok")
        loaddatasql = "load data local infile '" + winFileName8Field + "_tmpfile' into table " + tableNamewithoutNA + " fields terminated by '\\t'"
        
        shellstatment = "mysql -uroot -p"+config.password+" -D" + dbname.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        if a != 0:
            print("Util : loadWinDataIntoDB func os.system return not 0")
            exit(-1)
        print(shellstatment + ":ok")
        os.system("rm " + winFileName8Field + "_tmpfile")
        
        loaddatasql = "load data local infile '" + winFileName8Field + "' into table " + tableNametextValueForappendGeneName + " fields terminated by '\\t'"
        
        shellstatment = "mysql -h"+config.ip+" -uroot -p"+config.password+" -D" + dbname.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        
        if a != 0:
            print(shellstatment + ":failed")
            exit(-1)
        print(shellstatment + ":ok")    
        tempdbtools.operateDB("delete", "delete from " + tableNametextValueForappendGeneName + " where chrID='chrNo' and winNo='winNo' and winvalue='winvalue' ")    

        return chromOrder, tempdbtools, tableNamewithoutNA, tableNametextValueForappendGeneName 
    def appendGeneName(self, TranscriptGenetable, genomedbtools, winwidth, slideSize, outfileName,upextend=0, downextend=0,findNearestGene=(5,"m")):
        outfile = open(outfileName, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnps\twinvalue\tzvalue\tmark\tgeneName\ttrscptID", file=outfile)

        allwins = self.windbtools.operateDB("select", "select * from " + self.wintabletextvalueallwin )
        self.windbtools.operateDB("callproc", "mysql_sp_add_column", data=(self.dbname, self.wintabletextvalueallwin, "geneName", "varchar(128)", "default null"))
        self.windbtools.operateDB("callproc", "mysql_sp_add_column", data=(self.dbname, self.wintabletextvalueallwin, "trscptID", "varchar(128)", "default null"))  
        for win in allwins:
            region = (win[0], int(win[1]) * slideSize, int(win[1]) * slideSize + winwidth, win[1], win[5])
            geneNames = "";trscptIDs = ""
            recs=self.collectTrscptInWin(genomedbtools, TranscriptGenetable, region, upextend, downextend)
            for rec in recs:
                trscptIDs += rec[0].strip() + ";"
                if rec[2].strip() != "":
                    geneNames += (rec[2].strip() + ";")
            self.windbtools.operateDB("update", "update " + self.wintabletextvalueallwin + " set geneName = '" + geneNames[0:-1] + "', trscptID= '" + trscptIDs[0:-1] + "' where chrID= '" + win[0] + "' and winNo=" + win[1])
        #process outliers win
        
        total_outliers=findNearestGene[0]
        if findNearestGene[1]=="m":
            outlierwins=self.windbtools.operateDB("select","select * from "+ self.wintablewithoutNA+" order by zvalue desc limit 0,"+str(total_outliers))
        elif  findNearestGene[1]=="l":
            outlierwins=self.windbtools.operateDB("select","select * from "+ self.wintablewithoutNA+" order by zvalue asc limit 0,"+str(total_outliers))
        print(total_outliers,outlierwins)
        for win in outlierwins:
            region = (win[0], int(win[1]) * slideSize, int(win[1]) * slideSize + winwidth, win[1], win[5])
            geneNames = "";trscptIDs = ""
            recs=self.collectTrscptInWin(genomedbtools, TranscriptGenetable, region,upextend, downextend,True)
            for rec in recs:
                trscptIDs+=rec[0].strip() + ";"
                if rec[2].strip()!="":
                    geneNames+=(rec[2].strip() + ";")
            if recs==[]:
                geneNames+="top"+str(total_outliers)
                trscptIDs+="NA"
            self.windbtools.operateDB("update","update " + self.wintabletextvalueallwin + " set geneName = '" + geneNames[0:-1] + "', trscptID= '" + trscptIDs[0:-1] + "' where chrID= '" + win[0] + "' and winNo=" + win[1])
        allwins = self.windbtools.operateDB("select", "select * from " + self.wintabletextvalueallwin)
        for win in allwins:
            if win[-2] == "":
                if win[-1] == "":
                    print(*(win[:-2] + ("NA", "NA")), sep="\t", file=outfile)
                else:
                    print(*(win[:-2] + ("NA", win[-1])), sep="\t", file=outfile)
            else:
                print(*win, sep="\t", file=outfile)
        outfile.close()
    @staticmethod
    def collectTrscptInWin( genomedbtools, trscptableName, region, upextend=0, downextend=0,extendtodistal=0,treatallasPROTEINGENE=True):
        """select trscpt overlaped with the region
        reture a list of trscpts [tp_generecord1+overlapcode,tp_generecord2+overlapcode,,,]
        """
        trscptlist = []
        transcripttable = trscptableName
        chrID = region[0]
        Region_start = region[1]
        Region_end = region[2]
        """
        region=(chrom,Region_start,Region_end,Nwin,extremeValue,maxsnp,mixsnp)
        
        """
        selectType1OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos >= " + str(Region_start) + " and trscpt_end_pos <= " + str(Region_end)
        selectType2OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_start) + " and trscpt_end_pos > " + str(Region_end)
        selectType3OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_start) + " and trscpt_end_pos > " + str(Region_start) + " and trscpt_end_pos < " + str(Region_end)
        selectType4OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos > " + str(Region_start) + " and trscpt_start_pos < " + str(Region_end) + " and trscpt_end_pos > " + str(Region_end)
        selectType5OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_end_pos > " + str(Region_start - upextend) + " and trscpt_end_pos < " + str(Region_start)
        selectType6OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_end + downextend) + " and trscpt_start_pos > " + str(Region_end)
        #1
        findPROTEINGENE=False
        result = genomedbtools.operateDB("select", selectType1OverlapGenesql)
        for row in result:
            row += tuple([1])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
        #2
        result = genomedbtools.operateDB("select", selectType2OverlapGenesql)
        for row in result:
            row += tuple([2])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
        #3
        result = genomedbtools.operateDB("select", selectType3OverlapGenesql)
        for row in result:
            row += tuple([3])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
        #4
        result = genomedbtools.operateDB("select", selectType4OverlapGenesql)
        for row in result:
            row += tuple([4])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
        #5
        result = genomedbtools.operateDB("select", selectType5OverlapGenesql)
        for row in result:
            row += tuple([5])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
        #6
        result = genomedbtools.operateDB("select", selectType6OverlapGenesql)
        for row in result:
            row += tuple([6])
            trscptlist.append(row)
            if row[0].find("ENS")==0 or treatallasPROTEINGENE:
                findPROTEINGENE=True
            
        if not findPROTEINGENE and extendtodistal>max(upextend,downextend):
            result=genomedbtools.operateDB("select","select * from "+ transcripttable + " where transcript_ID regexp  'ENS' and  chrID='" + chrID +  "' and trscpt_end_pos < "+str(Region_start) + " order by trscpt_end_pos ")
            if len(result)!=0:#result is a list
                row=list(result[-1])
                if Region_start-int(row[6])<extendtodistal:
                    row[2]=("<"+str(Region_start-int(row[6])))+row[7]+row[2]
                    trscptlist.append(tuple(row)+tuple([7]))
            result=genomedbtools.operateDB("select","select * from "+ transcripttable + " where transcript_ID regexp  'ENS' and chrID='" + chrID +  "' and trscpt_start_pos > "+ str(Region_end)+" order by trscpt_start_pos desc")
            if len(result)!=0:#result is a list
                row=list(result[-1])
                if int(row[5])-Region_end<extendtodistal:
                    row[2]=(">"+str(int(row[5])-Region_end))+row[7]+row[2]
                    trscptlist.append(tuple(row)+tuple([8]))
        return trscptlist
def findTrscpt(winfile,outbedfilename,upextend,downextend,winwidth,slideSize,winType,morethan_lessthan,threshold_title_list=None,percentage=None,mergeNA=False,extendtodistal=0,anchorfile=None,found=False,mapfile=None):

    if percentage!=None and threshold_title_list!=None:
        print("-t conflict with -p")
        exit(-1)
    threshold_title_list
    if anchorfile:
#         winfile=standardseparately(anchorfile,winfile)
        winfilemark,winfilearrangement=mapWinvaluefileToChrOfReletiveSpecie(anchorfile, winfile, winwidth, slideSize, True,mapfile)
    else:
#         winfile=standardseparately(anchorfile,winfile)
        os.system("awk ' {if(NR=1){print $0"+'"\tmark"'+"}else{print $0"+'"\tunknown"'+"}}' "+winfile+">"+winfile+"marked.sexchromseperatestandard")
    winFileName8Field = winfile+"marked.sexchromseperatestandard"
    f=open(winFileName8Field,"r")
    title=re.split(r"\s+",f.readline().strip())
    f.close()
    Nocol=title.index(winType)+1
    re.search(r"[^/]*$",winFileName8Field).group(0)
    if re.search(r'^.*/',outbedfilename)!=None:
        path=re.search(r'^.*/',outbedfilename).group(0)
    else:
        a = os.popen("pwd")
        path=a.readline().strip()+"/"
        a.close()
    if found:
        outfileNameWINwithGENE=path+re.search(r"[^/]*$",winFileName8Field).group(0)+".wincopywithgene"
        return outfileNameWINwithGENE   
    outfile=open(outbedfilename+".bed.selectedgene",'w')
    print("chrNo\tRegion_start\tRegion_end\tNoofWin\textram"+winType+"\tminNoSNP\tmaxNoSNP\ttranscpt\toverlapcode\tgeneID",file=outfile)
    outfileNameWINwithGENE=path+re.search(r"[^/]*$",winFileName8Field).group(0)+".wincopywithgene"
    print(config.ip, config.username, config.password, config.genomeinfodbname)
    genomedbtools = dbm.DBTools(config.ip, config.username, config.password, config.genomeinfodbname)
    
    winGenome = WinInGenome(config.ghostdbname, winFileName8Field,Nocol)
    
    time.sleep(SLEEP_FOR_NEXT_TRY)
    selectWinNos="threshold method"
    totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintablewithoutNA)[0][0]  
#     selectWinNos = int(float(percentage) * totalWin)  
    if anchorfile:
        wherestatmentmt=" where (mark='autosome' and "+winType+">=" + threshold_title_list[0]+") or (mark='sexchromosome' and "+winType+">=" +threshold_title_list[-1]+")"
#         wherestatmentmp=" where 1 order by "+winType+" desc limit 0," + str(selectWinNos)
        wherestatmentlt=" where (mark='autosome' and "+winType+"<=" + threshold_title_list[0]+") or (mark='sexchromosome' and "+winType+"<=" +threshold_title_list[-1]+")"
#         wherestatmentlp=" where 1 order by "+winType+" asc limit 0," + str(selectWinNos)
    else:
        wherestatmentmt= " where 1 and "+winType+">=" + threshold_title_list[0]
#         wherestatmentmp=" where 1 order by "+winType+" desc limit 0," + str(selectWinNos)
        wherestatmentlt=" where "+winType+"!= 'NA' and "+winType+"<=" + threshold_title_list[0]
#         wherestatmentlp=" where 1 order by "+winType+" asc limit 0," + str(selectWinNos)
    winGenome.appendGeneName(config.TranscriptGenetable, genomedbtools, winwidth, slideSize, outfileNameWINwithGENE,upextend,downextend,(10,morethan_lessthan))
#    should be rewrite in a clear statment
    if percentage!=None:
        
        
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" desc limit 0," + str(selectWinNos))
            print("select * from "+winGenome.wintablewithoutNA + " where 1 order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" asc limit 0," + str(selectWinNos))
            print("select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" asc limit 0," + str(selectWinNos))
    elif threshold_title_list!=None:
        if morethan_lessthan=="m" or morethan_lessthan=="M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + wherestatmentmt)
            
        elif morethan_lessthan=="l" or morethan_lessthan=="L":
#             print("select", "select * from " + winGenome.wintablewithoutNA + " where "+winType+"!= 'NA' and "+winType+"<=" + threshold)
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + wherestatmentlt)
        selectWinNos=len(selectedWins)
    selectedWins.sort(key=lambda listRec:float(listRec[5]))
    if selectWinNos==0:
        outfile.close()
        print("selectWinNos==0")
        exit(0)
    print(outbedfilename+".bed.selectgene",selectWinNos,"~=",len(selectedWins),selectedWins[0],selectedWins[-1])
    selectedWinMap={}
    for win in selectedWins:
        if win[0] in selectedWinMap:
            selectedWinMap[win[0]].append(win)
        else:
            selectedWinMap[win[0]]=[win]

    selectedRegion={}

    for chrom in selectedWinMap:
        selectedWinMap[chrom].sort(key=lambda listRec: int(listRec[1]))
        selectedRegion[chrom]=[]
        mergedRegion=[selectedWinMap[chrom][0]]
        i=1
        while i < len(selectedWinMap[chrom]):
#             print(chrom,selectedWinMap[chrom][i])
#             try:
            if int(selectedWinMap[chrom][i-1][1])+1==int(selectedWinMap[chrom][i][1]) or int(selectedWinMap[chrom][i-1][1])*slideSize+winwidth>=int(selectedWinMap[chrom][i][1])*slideSize:#continues win
                mergedRegion.append(selectedWinMap[chrom][i])
            else:#not continues
                #process last region
                Region_start=int(mergedRegion[0][1])*slideSize
                Region_end=int(mergedRegion[-1][1])*slideSize+winwidth
                Nwin=len(mergedRegion)
                extremeValues=[]
                noofsnps=[]
                for e in mergedRegion:
                    if winType=="winvalue":
                        extremeValues.append(float(e[5]))
                    elif winType=="zvalue": 
                        extremeValues.append(float(e[6]))
                    noofsnps.append(int(e[4]))    
                        
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue=min(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=max(extremeValues)
                maxNoSNP=max(noofsnps)
                mixNoSNP=min(noofsnps)  
                selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue,mixNoSNP,maxNoSNP))
                #process this win
                mergedRegion=[selectedWinMap[chrom][i]]
            i+=1
#             except IndexError:
#                 print(i,len(selectedWinMap[chrom]),selectedWinMap[chrom])
#                 exit(-1)
        else:
            Region_start=int(mergedRegion[0][1])*slideSize
            Region_end=int(mergedRegion[-1][1])*slideSize+winwidth
            Nwin=len(mergedRegion)
            extremeValues=[]
            noofsnps=[]
            for e in mergedRegion:
                if winType=="winvalue":
                    extremeValues.append(float(e[5]))
                elif winType=="zvalue": 
                    extremeValues.append(float(e[6]))
                noofsnps.append(int(e[4]))
            if morethan_lessthan == "m" or morethan_lessthan == "M":
                extremeValue=min(extremeValues)
            elif morethan_lessthan == "l" or morethan_lessthan == "L":
                extremeValue=max(extremeValues)  
            maxNoSNP=max(noofsnps)
            mixNoSNP=min(noofsnps)                      
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue,mixNoSNP,maxNoSNP))
    if mergeNA!=False and int(mergeNA)>0:
        for chrom in selectedRegion:
            selectedRegion[chrom].sort(key=lambda listRec: int(listRec[1]))
            i=1
            idxlist_to_pop=[]
            while i <len(selectedRegion[chrom]):
                winNo_end=str(int(selectedRegion[chrom][i][1]/slideSize))
                winNo_start=str(int((selectedRegion[chrom][i-1][2]-winwidth)/slideSize))
                print("select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and  winNo<"+winNo_end)
                wincount_to_determine=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
                wincount_to_add=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintabletextvalueallwin + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
                if len(wincount_to_determine)==0 and len(wincount_to_add)<= int(mergeNA):
                    if morethan_lessthan == "m" or morethan_lessthan == "M":
                        extremeValue=min(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    elif morethan_lessthan == "l" or morethan_lessthan == "L":
                        extremeValue=max(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    maxNoSNP=max(selectedRegion[chrom][i][3],selectedRegion[chrom][i-1][3])
                    mixNoSNP=min(selectedRegion[chrom][i][3],selectedRegion[chrom][i-1][3])
                    selectedRegion[chrom][i]=(chrom,selectedRegion[chrom][i-1][1],selectedRegion[chrom][i][2],selectedRegion[chrom][i-1][3]+selectedRegion[chrom][i][3]+len(wincount_to_add),extremeValue,mixNoSNP,maxNoSNP)
                    idxlist_to_pop.append(i-1)
                i+=1
            else:
                idxlist_to_pop.reverse()
                for idx_to_pop in idxlist_to_pop:
                    selectedRegion[chrom].pop(idx_to_pop)
    else:
        for chrom in selectedRegion:
            selectedRegion[chrom].sort(key=lambda listRec: int(listRec[1]))
#    get final table
    print("getting final table")
    final_table={}
    for chrom in selectedRegion:
        for region in selectedRegion[chrom]:
            print(chrom,region)
            if extendtodistal>0:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,config.TranscriptGenetable,region,upextend,downextend,extendtodistal)
            else:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,config.TranscriptGenetable,region,upextend,downextend)
#process top outlier values
    print("fill bedselectedtable")
    for chrom in winGenome.chromOrder:
        if chrom not in selectedRegion:
            continue
        for region in selectedRegion[chrom]:
            if chrom.strip()==region[0].strip():
                tcpts=""
                tpcode=""
                gnames=""
                for tcpt in final_table[region]:
                    tcpts+=(tcpt[0]+",")
                    tpcode+=(str(tcpt[-1])+",")
                    if tcpt[2].strip()!="":
                        gnames+=(tcpt[2]+",")
                print("\t".join(map(str,region)),tcpts[:-1],tpcode[:-1],gnames[:-1],sep="\t",file=outfile)                  

    winGenome.windbtools.drop_table(winGenome.wintabletextvalueallwin)
    winGenome.windbtools.drop_table(winGenome.wintablewithoutNA)
    outfile.close()
    return outfileNameWINwithGENE

class genes():
    def __init__(self, gtfList, pos, RefSeqList, minintervalbetweengenes_basesperfaline=60):
        super().__init__()
        self.lastgenesRearpos = 0
        self.minintervalbetweengenes_basesperfaline = minintervalbetweengenes_basesperfaline
        self.geneOverlapList = self.getNearestGeneOverlapList(gtfList, pos)
        self.tscptSeqAllCds = {}
        self.cds_frame = {}  # {transcript_id:{cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}}
        
        
        for gene in self.geneOverlapList:
            
            genename = gene[0]
            self.tscptSeqAllCds[genename] = []
            self.cds_frame[genename] = {}  # {cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}

            cdsidx = 3
            for feature, elemStart, elemEnd, frame in gene[4:]:
                cdsidx += 1
                if feature == 'CDS':
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]  # ???如果不够呢
                elif  feature == "stop_codon":  # feature == 'start_codon' or
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]

#         print(genename,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,RefSeqList[0],len(RefSeqList))
#         print(self.cds_frame)
    def getNearestGeneOverlapList(self, gtfList, pos):
        """
        input:for a chrom,contain all transcript of this chrom
        gtfList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        return: the first gene that after the pos and the genes contain in or overlap with or contact with this gene indirect
        geneOverlapList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        order by "start"
        """
        if gtfList == None:
            return []
        for i in range(len(gtfList)):
            if gtfList[i][2] >= pos and (pos == 1 or self.lastgenesRearpos + self.minintervalbetweengenes_basesperfaline <= gtfList[i][2]):  # the distance between two genes must longer than minintervalbetweengenes_basesperfaline except overlapped genes
                geneOverlapList = [gtfList[i]]
                break
        else:
            if pos > gtfList[-1][3]:
                return []
            else:
                print("getNearestGeneOverlapList", pos, gtfList)
                exit(-1)
 
        furthest = gtfList[i][3]
        i += 1
        while len(gtfList) > i and furthest >= gtfList[i][2]:
            if gtfList[i][0] != geneOverlapList[-1][0]:
                geneOverlapList.append(gtfList[i])
            furthest = max(furthest, gtfList[i][3])
            i += 1
        print("getNearestGeneOverlapList", i, geneOverlapList, pos)
        self.lastgenesRearpos = geneOverlapList[-1][3]
        return geneOverlapList
    def getgeneConsensus(self, RefSeqList, idx_RefSeq, VcfList, idx_vcf, depthfile):
        """
        RefSeqList didn't changed
        """
        CodonTable = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
              'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
              'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
              'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
              'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
              'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
              'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
              'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
              'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
              'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
              'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
              'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
              'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
              'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
              'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
              'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
        curpos = RefSeqList[0] + idx_RefSeq
        tscptSeqAllCds_mut = {}
        ref_amino_seq = {}
        mutat_amino_seq = {}
        cns_append = ""
        originallen = {}  # just for test
#        indelatEdgeofCDS=open("indelatEdgeofCDS.txt",'w')
        for gene in self.geneOverlapList:
            genename = gene[0]
#            print(gene, self.cds_frame[genename], sep="\n", file=open("testgeneOverlapList.txt", 'a'))
            
            tscptSeqAllCds_mut[genename] = copy.deepcopy(self.tscptSeqAllCds[genename])
            originallen[genename] = len(tscptSeqAllCds_mut[genename])  # just for test
        
        while idx_vcf != -1 and idx_vcf != len(VcfList) and VcfList[idx_vcf][0] <= self.geneOverlapList[-1][3]:
            vcfpos = VcfList[idx_vcf][0];refalle = VcfList[idx_vcf][1];altalle = VcfList[idx_vcf][2]
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + (vcfpos - curpos)]))
            idx_RefSeq += (vcfpos - curpos)
            curpos = RefSeqList[0] + idx_RefSeq
            
            if re.search(r'[^a-zA-Z]', altalle) != None:  # contain ',' ie. multiple alle
                cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(refalle)]))
                idx_RefSeq += len(refalle)
                curpos = RefSeqList[0] + idx_RefSeq
                idx_vcf += 1
                continue
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(altalle)]))
            idx_RefSeq += len(refalle)  # here should still be refalle
            curpos = RefSeqList[0] + idx_RefSeq
            n_refbases = len(refalle);n_altbases = len(altalle)  # situation TAA     TA;     TTA     TTAAACTTCTATACTA;      C       T;    T       TATA;    ACG     A
# first for every variant making cns_append,and then substitute the seq in the cds seq,and finialy translate to protein 
            for gene in self.geneOverlapList:
                genename = gene[0]
                if gene[2] <= vcfpos and gene[3] >= vcfpos:
                    t4_indx = 3
                    for feature, elemStart, elemEnd, frame in gene[4:]:
                        t4_indx += 1
                        if feature == 'CDS' and vcfpos <= elemEnd and vcfpos >= elemStart:
                            if vcfpos + n_refbases - 1 > elemEnd or vcfpos + n_altbases - 1 > elemEnd:
                                print(VcfList[idx_vcf], genename, "indelatEdgeofCDS")
                                break
                            if n_refbases > n_altbases:  # situation TAA     TA;ACG     A
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases)] = list(altalle)
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)] = [' '] * (n_refbases - n_altbases)
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                            elif n_refbases < n_altbases:  # situation TTA     TTAAACTTCTATACTA;T       TATA;
                                if n_refbases == 1:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                else:
                                    tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1)] = list(altalle[0:(n_refbases - 1)])
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1] = altalle[(n_refbases - 1):]
                            else:  # n_refbases==n_altbases==1
                                try:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                except IndexError:
                                    print(self.cds_frame)
                                    print(genename, vcfpos, t4_indx, altalle, elemStart, feature, len(tscptSeqAllCds_mut[genename]))
                                    exit(-1)
       
            idx_vcf += 1
# 该翻译蛋白了吧 还有 看看长度一样不  将最后一个vcf记录之后的序列加入一致序列字符串
        cns_append += "".join(RefSeqList[idx_RefSeq:idx_RefSeq + (self.geneOverlapList[-1][3] - curpos) + 1])
        for gene in self.geneOverlapList:
            
            genename = gene[0]
#            ref_amino_seq[genename] = []
            mutat_amino_seq[genename] = []
#            mutationTypeList=[]
            if originallen[genename] != len(tscptSeqAllCds_mut[genename]):  # just for test
                print(self.tscptSeqAllCds[genename])
                print(tscptSeqAllCds_mut[genename])
                print("Util getgeneConsensus: length of tscptSeqAllCds changed,so there is a indel in the rearpart of the trcptSeq", genename)

            tscptSeqAllCds_mut_str = "".join(filter(lambda e:e.strip() != "", tscptSeqAllCds_mut[genename]))           
            if gene[1] == '+':
                tscptSeqAllCds_mut[genename] = list(tscptSeqAllCds_mut_str)
                for i in range(self.cds_frame[genename][sorted(self.cds_frame[genename].keys())[0]][0], len(tscptSeqAllCds_mut_str), 3):
#                    codon = "".join(self.tscptSeqAllCds[genename][i:i + 3]).lower()
                    codon_m = tscptSeqAllCds_mut_str[i:i + 3].lower()
#                    ref_amino_seq[genename].append(CodonTable[codon])
                    try:
                        mutat_amino_seq[genename].append(CodonTable[codon_m])
                    except KeyError:
                        mutat_amino_seq[genename].append('X')
            else:  # strand == '-'
                tscptSeqAllCds_Revr_Cmplm = complementary(self.tscptSeqAllCds[genename])
                tscptSeqAllCds_Revr_Cmplm.reverse()                
                tscptSeqAllCds_mut_str = list(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm = complementary(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm.reverse()
                tscptSeqAllCds_mut_str_Revr_Cmplm = "".join(tscptSeqAllCds_mut_str_Revr_Cmplm)
                tscptSeqAllCds_mut[genename] = list(tscptSeqAllCds_mut_str_Revr_Cmplm)
#                 print(genename,sorted(self.cds_frame[genename].keys()),len(tscptSeqAllCds_mut_str_Revr_Cmplm),3)
                for i in range(self.cds_frame[genename][sorted(self.cds_frame[genename].keys())[-1]][0], len(tscptSeqAllCds_mut_str_Revr_Cmplm), 3):
#                    codon = "".join(self.tscptSeqAllCds[genename][i:i+3]).lower()
                    codon_m = tscptSeqAllCds_mut_str_Revr_Cmplm[i:i + 3].lower()
#                    ref_amino_seq[genename].append(CodonTable[codon])
                    try:
                        mutat_amino_seq[genename].append(CodonTable[codon_m])
                    except KeyError:
                        mutat_amino_seq[genename].append('X')
#        testfile = open("animo_acid.txt", 'w')
#        for gene in self.geneOverlapList:#for test only
#            genename = gene[0]
#            print(">" + genename + "\n", file=testfile)
#            print("".join(ref_amino_seq[genename]), file=testfile)
#            print("\n" + "".join(mutat_amino_seq[genename]) + "\n", file=testfile)
#        testfile.close()
        return  tscptSeqAllCds_mut, mutat_amino_seq, cns_append, idx_vcf
def getGtfMap(gtfFileName, elementTypes=["CDS", "stop_codon"]):
    """protein_codingMap={chromNo:[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,],
               chromNo:[],,,,,,,,,,,,,,}
        chrtranscrpitididxMap{chromNo:{transcript_id:ttanscript_id_idx,transcript_id:ttanscript_id_idx,,,,,},
                                chromNo:{},chromNo:{},,,,}
        utrMap={chromNo:{{transcript_id:[("UTR",start,end),(),()],}}}
        allgeneSetMap={chromNo:{transcript_order:[],transcript_id:(strand,startpos,endpos),transcript_id:(strand,startpos,endpos),,,},chromNo:{transcript_id:(strand,startpos,endpos),,,,,},,,,}
    """
    try:
        
        protein_codingMap=pickle.load(open(gtfFileName+".protein_codingMap.landmine","rb"))
        utrMap=pickle.load(open(gtfFileName+".utrMap.landmine","rb"))
        allgeneSetMap=pickle.load(open(gtfFileName+".allgeneSetMap.landmine","rb"))
        return protein_codingMap,utrMap,allgeneSetMap
    except IOError:
        print("getGtfMap")
    gtfFileHandler = open(gtfFileName, 'r')
    protein_codingMap = {}
    chrtranscrpitididxMap = {}
    utrMap={}
    allgeneSetlist={}
#     gtfline = gtfFileHandler.readline()
    jumpout = False
    for getfirstcds in gtfFileHandler:
        if re.search(r"^#", getfirstcds) != None:
            print(getfirstcds)
            continue
        gtfColList = re.split(r'\s+', getfirstcds)
        chromNo = gtfColList[0].strip()
        protein_codingMap[chromNo] = []
        if "transcript_id" in gtfColList:
            transcript_id_idx = gtfColList.index("transcript_id") + 1
            gene_id = gtfColList.index("gene_id")
        else:
            continue
        print("transcript_id_idx", transcript_id_idx)
        transcript_id = re.search(r'\"(.*)\";', gtfColList[transcript_id_idx].strip()).group(1)
        countInChrom = 0
        if gtfColList[2].strip()=="transcript":
            if chromNo in allgeneSetlist:
                allgeneSetlist[chromNo].append((transcript_id,gtfColList[6],int(gtfColList[3]), int(gtfColList[4])))
            else:
                allgeneSetlist[chromNo]=[(transcript_id,gtfColList[6],int(gtfColList[3]), int(gtfColList[4]))]
        for elementType in elementTypes:
            if elementType == gtfColList[2].strip():
                jumpout = True
                protein_codingMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
        else:
            if gtfColList[2].strip()=="UTR" or "utr" in gtfColList[2].strip().lower():
                if int(gtfColList[3])!=int(gtfColList[4]):#ensembl's bug
                    utrMap[chromNo]={transcript_id:[("UTR",int(gtfColList[3]), int(gtfColList[4]))]}
            print(getfirstcds)
        if jumpout:
            break
    
    chrtranscrpitididxMap[chromNo] = {transcript_id:0}
    for gtfline in gtfFileHandler:
        gtfColList = re.split(r'\s+', gtfline)
        if gtfColList[2].strip()=="gene":continue
        transcript_id = re.search(r'\"(.*)\";', gtfColList[transcript_id_idx].strip()).group(1)
        chromNo = gtfColList[0].strip()
        if gtfColList[2].strip()=="transcript":
            if chromNo in allgeneSetlist:
                allgeneSetlist[chromNo].append((transcript_id,gtfColList[6],int(gtfColList[3]), int(gtfColList[4])))
            else:
                allgeneSetlist[chromNo]=[(transcript_id,gtfColList[6],int(gtfColList[3]), int(gtfColList[4]))]
        for elementType in elementTypes:
            if elementType == gtfColList[2].strip():
                break
        else:
            if (gtfColList[2].strip()=="UTR"  or "utr" in gtfColList[2].strip().lower()) and int(gtfColList[3])!=int(gtfColList[4]):#ensembl's bug
                if chromNo in utrMap:
                    if transcript_id in utrMap[chromNo]:
                        utrMap[chromNo][transcript_id].append(("UTR",int(gtfColList[3]), int(gtfColList[4])))
                    else:
                        utrMap[chromNo][transcript_id]=[("UTR",int(gtfColList[3]), int(gtfColList[4]))]
                else:
                    utrMap[chromNo]={transcript_id:[("UTR",int(gtfColList[3]), int(gtfColList[4]))]}
            continue
        
        if chromNo in protein_codingMap:
            if transcript_id in chrtranscrpitididxMap[chromNo].keys():
                tanscript_id_idx = chrtranscrpitididxMap[chromNo][transcript_id]
                protein_codingMap[chromNo][tanscript_id_idx].append((gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7]))
                protein_codingMap[chromNo][tanscript_id_idx][2] = min(protein_codingMap[chromNo][tanscript_id_idx][2], int(gtfColList[3]))
                protein_codingMap[chromNo][tanscript_id_idx][3] = max(protein_codingMap[chromNo][tanscript_id_idx][3], int(gtfColList[4]))
            else:
                protein_codingMap[chromNo].append([transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])])
                chrtranscrpitididxMap[chromNo][transcript_id] = len(protein_codingMap[chromNo]) - 1
        else:
             protein_codingMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
             chrtranscrpitididxMap[chromNo] = {transcript_id:0}
    else:
        pass                 
    gtffilepath = re.search(r"^.*[/]", gtfFileName).group(0)
    testfile = open(gtffilepath + "protein_codingMap.sort.txt", 'w')
    allgeneSetMap={}
    for chromNo in allgeneSetlist.keys():
        allgeneSetMap[chromNo]={"transcript_order":[]}
        allgeneSetlist[chromNo].sort(key=lambda listRec:listRec[2])
        for tp_id,strand,startpos,endpos in allgeneSetlist[chromNo]:
            allgeneSetMap[chromNo]["transcript_order"].append(tp_id)
            allgeneSetMap[chromNo][tp_id]=(strand,startpos,endpos)
        
    for chromNo in protein_codingMap.keys():
        protein_codingMap[chromNo].sort(key=lambda listRec:listRec[2])
        # 先按照转录本起始坐标排序，下面是对转录本内元件排序，不过是什么排序方法忘记了，仔细读一下吧
        for j in range(len(protein_codingMap[chromNo])):
            for t4_indx in range(5, len(protein_codingMap[chromNo][j])):
                t4_key = protein_codingMap[chromNo][j][t4_indx]
                t4_indxp = t4_indx - 1
                while t4_indxp >= 4 and protein_codingMap[chromNo][j][t4_indxp][1] > t4_key[1]:
                    protein_codingMap[chromNo][j][t4_indxp + 1] = protein_codingMap[chromNo][j][t4_indxp]
                    t4_indxp = t4_indxp - 1
                else:
                    protein_codingMap[chromNo][j][t4_indxp + 1] = t4_key
#            print(protein_codingMap[chromNo][j],protein_codingMap[chromNo][j][t4_indxp][1] > t4_key[1],t4_indxp)
        print(chromNo, "num of transcrpit:", len(protein_codingMap[chromNo]), file=testfile)
        for i in range(len(protein_codingMap[chromNo])):
#            print(protein_codingMap[chromNo][i][0],protein_codingMap[chromNo][i][1],protein_codingMap[chromNo][i][2],protein_codingMap[chromNo][i][3])
            for k in range(len(protein_codingMap[chromNo][i])):
                print(protein_codingMap[chromNo][i][k], file=testfile)
    testfile.close()
    gtfFileHandler.close()
#     testutr=open("testURT","w")
    for chrom in utrMap.keys():
#         print(chrom,file=testutr)
        for tpid in utrMap[chrom].keys():
            utrMap[chrom][tpid].sort(key=lambda listRec:listRec[1])
#             print(tpid,file=testutr)
#             print(utrMap[chrom][tpid],file=testutr)
#     testutr.close()
    pickle.dump(protein_codingMap,open(gtfFileName+".protein_codingMap.landmine", 'wb'))
    pickle.dump(utrMap,open(gtfFileName+".utrMap.landmine", 'wb'))
    pickle.dump(allgeneSetMap,open(gtfFileName+".allgeneSetMap.landmine", 'wb'))
    return protein_codingMap,utrMap,allgeneSetMap
def getNearestGenegroup(gtfList, pos):
    """
    input:for a chrom,contain all transcript of this chrom
    gtfList    =    [[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
    return: the first gene that after the pos and the genes contain in or overlap with or contact with this gene indirect
    geneOverlapList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
    order by "start"
    """
    
    if gtfList == None:
        return []
    for i in range(len(gtfList)):
        if gtfList[i][2] > pos: 
            geneOverlapList = [gtfList[i]]
            break
    else:
        if pos > gtfList[-1][3]:
            print("what's wrong")
            return []
        else:
            exit(-1)

    furthest = gtfList[i][3]
    i += 1
    while len(gtfList) > i and furthest >= gtfList[i][2]:
        if gtfList[i][0] != geneOverlapList[-1][0]:  # this judgement maybe not need
            geneOverlapList.append(gtfList[i])
        furthest = max(furthest, gtfList[i][3])
        i += 1
    return furthest, geneOverlapList
def getSNPrecInCDS(codon_idx_in_tscptSeqAllCds, lenOftscpt, codon, codon_m, cds_frame_ONETRSCPT, gtftrscpt_list):
    for i in range(3):
        if codon[i] != codon_m[i]:
            snp_idx_in_tscptALLCds = i + codon_idx_in_tscptSeqAllCds
            ref_base = codon[i].upper()
            alt_base = codon_m[i].upper()
    if gtftrscpt_list[1] == "+":
        for cds_idx in sorted(cds_frame_ONETRSCPT.keys()):
            if cds_frame_ONETRSCPT[cds_idx][1] > snp_idx_in_tscptALLCds:
                cds_idx -= 1
                break
        snp_pos = snp_idx_in_tscptALLCds - cds_frame_ONETRSCPT[cds_idx][1] + gtftrscpt_list[cds_idx][1]
    elif gtftrscpt_list[1] == "-":
        for cds_idx in sorted(cds_frame_ONETRSCPT.keys()):
            if cds_frame_ONETRSCPT[cds_idx][1] > lenOftscpt - snp_idx_in_tscptALLCds - 1:
                cds_idx -= 1
                break
        snp_pos = lenOftscpt - snp_idx_in_tscptALLCds - 1 - cds_frame_ONETRSCPT[cds_idx][1] + gtftrscpt_list[cds_idx][1]
        ref_base = complementary(ref_base).upper()
        alt_base = complementary(alt_base).upper()
    return snp_pos, ref_base, alt_base
def complementary(seq):
    """
    ['tg', 'a', 't', 'g', 'c', 'acacacgatg', 'ctttttcccccccc', 'c', 'c', 'a', 'a', 'aaagagagagacagaaaaaggc', 'atatcgactg', 'catcga']
    reverse to 
    ['ac', 't', 'a', 'c', 'g', 'tgtgtgctac', 'gaaaaagggggggg', 'g', 'g', 't', 't', 'tttctctctctgtctttttccg', 'tatagctgac', 'gtagct']
    """
    newseq = []
    for i in range(0, len(seq)):
        if seq[i].lower() == 'a':
            newseq.insert(i, 't')
        elif seq[i].lower() == 't':
            newseq.insert(i, 'a')
        elif seq[i].lower() == 'c':
            newseq.insert(i, 'g')
        elif seq[i].lower() == 'g':
            newseq.insert(i, 'c')
        elif len(seq[i]) > 1:
            newseq.insert(i, complementary(seq[i]))
        else:
            newseq.insert(i, seq[i])
    if isinstance(seq, str):
        newseq = "".join(newseq)
    return newseq
def getGeneGrouplist(gtfList):
    """
        geneGrouplist=[genegroup1,genegroup2,genegroup3,....]
                     =[[furthest,geneOverlapLists1],[furthest,geneOverlapLists2],[furthest,geneOverlapLists3],.....]
                     =[[furthest,[],[],[],[],....],[furthest,[],[],[],[],....],[furthest,[],[],[],[],....],.....]
                     =[[furthest,[transcript_id,strand,start,end,(),(),(),...],[transcript_id,strand,start,end,(),(),(),...],....],[furthest,[transcript_id,strand,start,end,(),(),(),...],[transcript_id,strand,start,end,(),(),(),...],....],[furthest,[],[],...],[],........]
    """
    newgtfList = copy.deepcopy(gtfList)
    geneGrouplist = []
    curpos = newgtfList[0][2] - 1
    while curpos < newgtfList[-1][2]:
        furthest, geneGroup = getNearestGenegroup(newgtfList, curpos)
        curpos = furthest + 1
        geneGroup.insert(0, furthest)
        geneGrouplist.append(geneGroup)
        
    return geneGrouplist
def make_getElemBed(elementfold,targetseqnamesubstr,pathtoblastn,reffa):
    """
    targetseqnamesubstr is the str before the first space ,after the >
    """
    allseqtobed={}#{chrID:[(sstart,send,elem,qstart,qend,revcom,len),(sstart,send,elem,qstart,qend,revcom,len),,,],,,,}
    if elementfold.endswith("/") or elementfold.endswith("\\"):
        elementfold=elementfold[:-1]
    if os.path.isfile(elementfold+"/"+targetseqnamesubstr+".bed"):
        bedfile=open(elementfold+"/"+targetseqnamesubstr+".bed","r")
        bedfile.readline()#title
        for bedline in  bedfile:
            bedlinelist=re.split(r"\t+",bedline)
            if bedlinelist[0].strip() in allseqtobed:
                allseqtobed[bedlinelist[0].strip()].append((int(bedlinelist[1]),int(bedlinelist[2]),bedlinelist[3],int(bedlinelist[4]),int(bedlinelist[5]),bedlinelist[6],int(bedlinelist[7]),int(bedlinelist[8])))
            else:
                allseqtobed[bedlinelist[0].strip()]=[(int(bedlinelist[1]),int(bedlinelist[2]),bedlinelist[3],int(bedlinelist[4]),int(bedlinelist[5]),bedlinelist[6],int(bedlinelist[7]),int(bedlinelist[8]))]
        bedfile.close()
        return allseqtobed
    randomstr=config.random_str()
    targetseqnamesubstr_lenmap={}
    if targetseqnamesubstr=="none":
        shellstatment=pathtoblastn+" -query "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".fa"+" -task blastn -db "+reffa+" -out "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".blastout -outfmt 7 -num_alignments 10 -num_threads 6"
    queryfafile=open(elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".collectionfas",'w')
    i=0
    for elem in os.listdir(path=elementfold):
        path = elementfold + "/" + elem
        
        if (not os.path.isdir(path)) and (path.endswith("fa") or path.endswith("fasta")):#True is fa file
            print(path,i)
            i+=1
            
            if targetseqnamesubstr.lower().strip()=="none":
                pathfile=open(path,"r")
                for line in pathfile:
                    print(line.strip(),file=queryfafile)
                    if line.startswith(">"):
                        seqname=line.strip()
                    else:
                        targetseqnamesubstr_lenmap[seqname[1:]]=len(line.strip())
#                 print(targetseqnamesubstr_lenmap)
                pathfile.close()
            else:
                muscleout_seqgenerator=SeqIO.parse(path,"fasta")
                for seq_rec in muscleout_seqgenerator:
                    if seq_rec.id==targetseqnamesubstr:
                        seqstr="".join(seq_rec.seq).replace("-", "")
                        print(">"+elem,file=queryfafile)
    #                     allseqtobed[elem]=[]
                        targetseqnamesubstr_lenmap[elem]=len(seqstr)
                        print(seqstr,file=queryfafile)
                        break
                else:
                    print(targetseqnamesubstr,"dosenot exist",elem)
    queryfafile.close()
    shellstatment=pathtoblastn+" -query "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".collectionfas"+" -task blastn -db "+reffa+" -out "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".blastout -outfmt 7 -num_alignments 10 -num_threads 6"
    print(shellstatment)
    a=os.system(shellstatment)
    if a!=0:
        print("error")
        exit(-1)
    blastout=open(elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".blastout","r")
    for line in blastout:
        if re.search(r"^#",line)!=None:
            lastblastlen=None
            continue
        linelist=re.split(r"\s+",line)
        blastlen=int(linelist[3])
        if lastblastlen==None or (blastlen>lastblastlen-10 or  blastlen*0.95>=lastblastlen):
            fafilename=linelist[0]
            chrom = linelist[1]
            sstartpos = int(linelist[8])
            sendpos = int(linelist[9])
            revcom="forward"
            if sstartpos > sendpos:
                temp = sstartpos
                sstartpos = sendpos
                sendpos = temp
                revcom="revcom"
            qstartpos=int(linelist[6])
            qendpos=int(linelist[7])
            total_bases=targetseqnamesubstr_lenmap[fafilename]
            gap_open=int(linelist[5])
            if chrom in allseqtobed:
                allseqtobed[chrom].append((sstartpos,sendpos,fafilename,qstartpos,qendpos,revcom,total_bases,gap_open))
            else:
                allseqtobed[chrom]=[(sstartpos,sendpos,fafilename,qstartpos,qendpos,revcom,total_bases,gap_open)]
            lastblastlen=blastlen
    bedfile=open(elementfold+"/"+targetseqnamesubstr+".bed","w")
    print("chrNo","Region_start","Region_end","fastafilename","startbase","endbase","revcom_forward","total_bases","gap_open",sep="\t",file=bedfile)
    for chrom in allseqtobed.keys():
        allseqtobed[chrom].sort(key=lambda listRec:listRec[1])
        for startpos,endpos,fafilename,qs,qe,revcom,total_bases,gap_open in allseqtobed[chrom]:
            print(chrom,startpos,endpos,fafilename,qs,qe,revcom,total_bases,gap_open,sep="\t",file=bedfile)
    blastout.close()
    bedfile.close()
    return allseqtobed
    os.system("rm "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".fa "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".blastout")
    