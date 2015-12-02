import copy
import os
import re
import time

from Bio import SeqIO
from scipy import stats

from NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


SLEEP_FOR_NEXT_TRY=3
def GOenrichment(gotablefile,outpre,genelist=None,trscptlist=None):
    gotablefile=open(gotablefile,'r')
    title = gotablefile.readline()
    titlelist= [e.strip().lower() for e in re.split(r"\t",title)]
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
def findTrscpt(winfile,outbedfilename,upextend,downextend,winwidth,slideSize,winType,morethan_lessthan,threshold=None,percentage=None,mergeNA=False,extendtodistal=0,found=False):

    if percentage!=None and threshold!=None:
        print("-t conflict with -p")
        exit(-1)
    winFileName7Field = winfile
    re.search(r"[^/]*$",winFileName7Field).group(0)
    if re.search(r'^.*/',outbedfilename)!=None:
        path=re.search(r'^.*/',outbedfilename).group(0)
    else:
        a = os.popen("pwd")
        path=a.readline().strip()+"/"
        a.close()
    if found:
        outfileNameWINwithGENE=path+re.search(r"[^/]*$",winFileName7Field).group(0)+".wincopywithgene"
        return outfileNameWINwithGENE   
    outfile=open(outbedfilename+".bed.selectedgene",'w')
    print("chrNo\tRegion_start\tRegion_end\tNoofWin\textram"+winType+"\ttranscpt\toverlapcode\tgeneID",file=outfile)
    outfileNameWINwithGENE=path+re.search(r"[^/]*$",winFileName7Field).group(0)+".wincopywithgene"
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname) 
    winGenome = Util.WinInGenome(Util.ghostdbname, winFileName7Field)
    time.sleep(SLEEP_FOR_NEXT_TRY)
    winGenome.appendGeneName(Util.TranscriptGenetable, genomedbtools, winwidth, slideSize, outfileNameWINwithGENE,upextend,downextend,(10,morethan_lessthan))
    selectWinNos="threshold method"
    if percentage!=None:
        totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintablewithoutNA)[0][0]
        selectWinNos = int(float(percentage) * totalWin)
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" desc limit 0," + str(selectWinNos))
            print("select * from "+winGenome.wintablewithoutNA + " where 1 order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" asc limit 0," + str(selectWinNos))
            print("select * from " + winGenome.wintablewithoutNA + " where 1 order by "+winType+" asc limit 0," + str(selectWinNos))
    elif threshold!=None:
        if morethan_lessthan=="m" or morethan_lessthan=="M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 and "+winType+">=" + threshold)
        elif morethan_lessthan=="l" or morethan_lessthan=="L":
            print("select", "select * from " + winGenome.wintablewithoutNA + " where "+winType+"!= 'NA' and "+winType+"<=" + threshold)
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 and "+winType+"<=" + threshold)
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
            if int(selectedWinMap[chrom][i-1][1])+1==int(selectedWinMap[chrom][i][1]):#continues win
                mergedRegion.append(selectedWinMap[chrom][i])
            else:#not continues
                #process last region
                Region_start=int(mergedRegion[0][1])*slideSize
                Region_end=int(mergedRegion[-1][1])*slideSize+winwidth
                Nwin=len(mergedRegion)
                extremeValues=[]
                for e in mergedRegion:
                    if winType=="winvalue":
                        extremeValues.append(float(e[5]))
                    elif winType=="zvalue": 
                        extremeValues.append(float(e[6]))
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue=min(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=max(extremeValues)
                selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
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
            for e in mergedRegion:
                if winType=="winvalue":
                    extremeValues.append(float(e[5]))
                elif winType=="zvalue": 
                    extremeValues.append(float(e[6]))
            if morethan_lessthan == "m" or morethan_lessthan == "M":
                extremeValue=min(extremeValues)
            elif morethan_lessthan == "l" or morethan_lessthan == "L":
                extremeValue=max(extremeValues)            
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
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
                    selectedRegion[chrom][i]=(chrom,selectedRegion[chrom][i-1][1],selectedRegion[chrom][i][2],selectedRegion[chrom][i-1][3]+selectedRegion[chrom][i][3]+len(wincount_to_add),extremeValue)
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
    final_table={}
    for chrom in selectedRegion:
        for region in selectedRegion[chrom]:
            if extendtodistal>0:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,Util.TranscriptGenetable,region,upextend,downextend,extendtodistal)
            else:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,Util.TranscriptGenetable,region,upextend,downextend)
#process top outlier values

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
    
    randomstr="UrZTOpSJ"#Util.random_str()
    targetseqnamesubstr_lenmap={}
    queryfafile=open(elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".fa",'w')
    for elem in os.listdir(path=elementfold):
        path = elementfold + "/" + elem
        if (not os.path.isdir(path)) and (path.endswith("fa") or path.endswith("fasta")):#True is fa file
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
    shellstatment=pathtoblastn+" -query "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".fa"+" -task blastn -db "+reffa+" -out "+elementfold+"/"+randomstr+"_"+targetseqnamesubstr+".blastout -outfmt 7 -num_alignments 10 -num_threads 6"
    print(shellstatment)
#     a=os.system(shellstatment)
#     if a!=0:
#         print("error")
#         exit(-1)
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
    