'''
Created on 2018年5月24日

@author: Dr.liu
'''
from multiprocessing.dummy import Pool
import os,re,sys
from src.NGS.BasicUtil import *
SNPIDmissedinextract=[]
def extarctBlastOutvariablelenwaper(a):
    extarctBlastOutvariablelen(**a)
def extarctBlastOutvariablelen(BlastOutFile,faMap,fromposOfBlastf=0):
    global SNPIDmissedinextract
#     me=open(missedinextract,'r')
    ospid=os.getpid()
    f=open(BlastOutFile+"missedinextract"+str(ospid)+".fa",'a')

    for snpid in SNPIDmissedinextract:
        print(">"+snpid+"\n"+faMap[snpid],file=f)

    f.close()
        
def extarctBlastOut(BlastOutFile,queryFaFile,filterlen=100,mismatch=6):
    print(" query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
    ambigousfile=open(BlastOutFile+".ambigous.blastrec",'w')
    qfa=open(queryFaFile,'r')
    faMap={}
    subid=os.getpid()
    for line in qfa:
        if line[0]==">":
            queryID=line[1:].strip()
        elif len(line.strip())>0:
            faMap[queryID]=line.strip()
    allIDf=open('allSNPID'+str(subid),'w')
    for ID in faMap.keys():
        print(ID,file=allIDf)
    allIDf.close()
    a = os.popen("awk '$1!~/^#/ && $4>="+str(filterlen)+" && $5<="+str(mismatch)+"  {print $0}' " + BlastOutFile)
    posfile=open(BlastOutFile+".pos",'w')

    revcom="+"
#    initial
    hit = a.readline()
    hitlist = re.split(r"\s+", hit)

    sendpos = int(hitlist[9])
    sstartpos = int(hitlist[8])
    qstartpos = int(hitlist[6])
    blastlen=int(hitlist[3])
    
    
    if sstartpos > sendpos:
        temp = sstartpos
        sstartpos = sendpos
        sendpos = temp
        revcom="-"
    
    lastsnpID = hitlist[0]
    flanklen=faMap[lastsnpID.strip()].index("N")
    snp_loc_s=sstartpos+flanklen
    chrom= hitlist[1]

    for hit in a:
        hitlist = re.split(r"\s+", hit)
        if lastsnpID==hitlist[0]:
            print(hit.strip(),file=ambigousfile)
            continue#skip this query blast result
        else:
            print(lastsnpID,chrom,snp_loc_s,file=posfile)
            chrom = hitlist[1]
            sstartpos = int(hitlist[8])
            sendpos = int(hitlist[9])
            qstartpos = int(hitlist[6])
            blastlen=int(hitlist[3])
            lastsnpID = hitlist[0]
            flanklen=faMap[lastsnpID.strip()].index("N")
            snp_loc_s=sstartpos+flanklen
            if sstartpos > sendpos:
                temp = sstartpos
                sstartpos = sendpos
                sendpos = temp
                revcom="-"
            else:
                revcom="+"
            
    posfile.close()
    print("postion writing done")
    
    os.system("awk '{print $1}' "+BlastOutFile+".pos > extractedexactly"+str(subid))
    
    os.system("grep -vFf extractedexactly"+str(subid)+" allSNPID"+str(subid)+" > missedinextractmp"+str(subid))
    
#     a = os.popen("less missedinextractmp"+str(subid)+"|wc -l ")
#     totalsnpforcount=int(a.readline().strip())
#     a.close()
    f = open("missedinextractmp"+str(subid),'r')
    for line in f:
        SNPIDmissedinextract.append(line.strip())
    f.close()
    parameterstuples_list=[]
    f=open(BlastOutFile,'r')
    blastoutfilendpos=f.seek(0,os.SEEK_END)
    f.close()
    d=int(blastoutfilendpos/1)+1
    j=0
    for i in range(0,blastoutfilendpos,d):
#         print("sed -n '"+str(i)+","+str(i+d-1)+"p' missedinextractmp"+str(subid)+" >missedinextractmp"+str(subid)+str(j))
#         os.system("sed -n '"+str(i)+","+str(i+d-1)+"p' missedinextractmp"+str(subid)+" >missedinextractmp"+str(subid)+str(j))
        parameterstuples_list.append({"BlastOutFile":BlastOutFile,"faMap":faMap,"fromposOfBlastf":i})
        j+=1
    pool=Pool(1)
    pool.map(extarctBlastOutvariablelenwaper,parameterstuples_list)
    pool.close()
    pool.join()
    ambigousfile.close()
    print("finish")
if __name__ == '__main__':
    print(len(sys.argv))
    infile1map={}
    if len(sys.argv)==2:
        print("print duprecs")
        seqidx=2
        f=open(sys.argv[1],'r');f.readline();ofo=open(sys.argv[1]+"dup",'w')
        for line in f:
            recl=re.split(r"\t",line.strip())
            seql=re.split(r"\[.+\]",recl[seqidx].strip())#############
            flankseq=list(seql[0]+seql[1]);flankseq.reverse()
            compflankseq=geneUtil.complementary(flankseq)
            print(seql[0]+seql[1],compflankseq)
            if seql[0]+seql[1] in infile1map :
                print(*recl,sep="\t",file=ofo)
                print(*infile1map[seql[0]+seql[1]],sep="\t",file=ofo)
            elif ("".join(compflankseq)).upper() in infile1map:
                print(*recl,"recomp",sep="\t",file=ofo)
                print(*infile1map[("".join(compflankseq)).upper()],sep="\t",file=ofo)
            else:
                infile1map[seql[0]+seql[1]]=recl
        f.close();ofo.close()
        print(recl[2].strip(),seql[0]+seql[1])
        print("notice! when you use grep -w/vFf to get the same/different file, notice some ID may not only occur in the ID column,but also duplicates cols ")
        print("some times may print one records multiple times")
        exit()
    ##################################################################################################################################
    dupmap={};dupids=[]
    if len(sys.argv)==3 and os.path.exists(sys.argv[2]):
        print("remove dup,input command: python recswithdup duprecs")
        tilingorder_idx=2
        f=open(sys.argv[1],'r');dupf=open(sys.argv[2],'r');remdupf=open(sys.argv[1]+"redup",'w');print(f.readline().strip(),file=remdupf)
        for line in dupf:
            recl=re.split(r"\t",line.strip())
            if recl[0] in dupmap:
                dupmap[recl[0]].append(recl)
            else:
                dupids.append(recl[0])
                dupmap[recl[0]]=[recl]
        print(len(dupids))
        for line in f:
            recl=re.split(r"\t",line.strip())
            if  recl[0] not in dupids:#(recl[0] in dupmap and recl[-1]=="5") or
                print(*recl,sep="\t",file=remdupf)
            else:
                recl[tilingorder_idx]="10"
                print(*recl,sep="\t",file=remdupf)
        remdupf.close();f.close();dupf.close()
        exit()
    if len(sys.argv)==4 and os.path.exists(sys.argv[1]) and os.path.exists(sys.argv[2]) and os.path.exists(sys.argv[3]):
        allscoredf=open(sys.argv[2],'r')
        lowscoref=open(sys.argv[3],'r')
        lowscores={}
        for line in lowscoref:
            chrNo=re.split(r",",line.strip())[0].split("_")[0];pos=re.split(r",",line.strip())[0].split("_")[-1]
            try:
                lowscores[chrNo].append(pos)
            except:
                lowscores[chrNo]=[pos]
        lowscoref.close()
        idf=open(sys.argv[1],'r')        
        for line in idf:
            if re.search(r"^Locus_Name",line)!=None:
                print(line,"titlepassed")
                break
        IDlist=[]
        for line in idf:
            IDlist.append(re.split(r",",line.strip())[0])
        idf.close()
        fw=open("asdfaf",'w')
        for line in allscoredf:
            if re.search(r"^Locus_Name",line)!=None:
                print(line,"titlepassed")
                break
        cchr="";scorelist=[];poslist=[]
        for line in allscoredf:
            recl=re.split(r",",line.strip())
            if recl[0].split("_")[0]==cchr:
                poslist.append(recl[0].split("_")[-1])
                scorelist.append(float(recl[9]))
            else:
                if cchr in lowscores.keys():
                    lowscores[cchr].sort()
                    for pos in lowscores[cchr]:#for every sites need to sub
                        print(cchr+"_"+pos)
                        print(cchr+"_"+pos,end="\t",file=fw)
                        idx=poslist.index(pos);i=1
                        while idx!=0 or i!=-1:
                            idx+=i
                            if len(poslist)!=idx:
                                if cchr+"_"+poslist[idx] not in IDlist and float(scorelist[idx]) >0.6:
                                    print(cchr+"_"+poslist[idx],file=fw);IDlist.append(cchr+"_"+poslist[idx]);break#sub find
                            else:
                                idx=poslist.index(pos);i=-1
                else:
                    print(cchr,"sdfasdf")
                print("",file=fw)
                scorelist=[];cchr=recl[0].split("_")[0];poslist=[]
        else:
            if cchr in lowscores.keys():
                lowscores[cchr].sort()
                for pos in lowscores[cchr]:#for every sites need to sub
                    print(cchr+"_"+pos)
                    print(cchr+"_"+pos,end="\t",file=fw)
                    idx=poslist.index(pos);i=1
                    while idx!=0 or i!=-1:
                        idx+=i
                        if len(poslist)!=idx:
                            if cchr+"_"+poslist[idx] not in IDlist and float(scorelist[idx]) >0.6:
                                print(cchr+"_"+poslist[idx],file=fw);IDlist.append(cchr+"_"+poslist[idx]);break#sub find
                        else:
                            idx=poslist.index(pos);i=-1                      
        fw.close()                        
        exit()                            
###add priority to scored file################################################################################################################################
    if len(sys.argv)==5:
        print("python chip_design.py scoredfile winsize I/A thresholdfor_I")
        print("add prority accord tiling_order and win; i.e repriority")
        print("thresholdfor_I argv[4] is for I only, but should give any value for A ")
        if sys.argv[3].strip().upper()=="A":
            musctincludes=set()
            """this code is for xiaomi chip, there 417 sites that customer request to include in this chip
            mustincludef=open("lkjltail603uniq417",'r');
            for line in mustincludef:
                mustinl=re.split(r"\t",line.strip())
                musctincludes.add(mustinl[1])
            mustincludef.close()
            """
            
            
            dupseqmap={}
        #     dupf=open("fillgapMergeOldmarker.txtdup",'r');
            """{seq1:[highestpriority,id1,id2,id3],seq2:[highestpriority,id1,id2],seq3:[highestpriority,id1,id2],,}"""
            os.system("sed -n '1,1p' " +sys.argv[1]+" > titlelinetemp")
        #     os.system("sed -n '2,$p' " +sys.argv[1]+"|sort -t$'\t'  -k5,5n -k6,6n - |cat titlelinetemp -  >"+sys.argv[1]+".sorted");sys.argv[1]=sys.argv[1]+".sorted"
            
            print(sys.argv[1])
            f=open(sys.argv[1],'r');title=re.split(r"\t",f.readline().strip())
            tidx=title.index("tiling_order");curchr=title[4];print("curchr:",curchr,"tilingorder idx:",tidx)
            
            """ this code is for xiaomi specified to remove duplicate sites, has been instead by first code segment "len(sys.argv)==2:" 
            
            
            for drec in dupf:
                drecl=re.split(r"\t",drec.strip())
                seql=re.split(r"\[.+\]",drecl[2].strip())
                
                seqmerge=seql[0]+seql[1]
                if seqmerge in dupseqmap:
                    dupseqmap[seqmerge].append(drecl[1])
                    dupseqmap[seqmerge][0]=max(dupseqmap[seql[0]+seql[1]][0],int(drecl[3]))
                else:
                    dupseqmap[seqmerge]=[int(drecl[3]),drecl[1]]
            """    
            
        #     print(dupseqmap);exit()
            win = Util.Window()
            ofo=open(sys.argv[1]+sys.argv[2],'w')
            print(*title,"priority",sep="\t",file=ofo)
            addprortycaculator = Caculators.Caculator_addpriority(of=ofo,tilingorderidx=tidx,best_recommendation=title.index("best_recommendation"),rmdupmap=dupseqmap,mustin=musctincludes)
            addprortycaculator.addVIPidx=title.index("importance")
            recs=[];count=0;tcount=0
            for line in f:
                tcount+=1
                recl=re.split(r"\t",line.strip())
                
                if curchr==recl[4]:
                    currentchrLen=int(recl[5])
        #             print("collect rec in a win")
                    recs.append([int(recl[5])]+recl)
                elif curchr!="cust_chr":
        #             print(recl,"sliding win",len(recs),currentchrLen)
                    
                    count+=len(recs)
                    if count!=tcount-1:
                        print(line,count,tcount);exit(-1)
                    win.slidWindowOverlap(recs, currentchrLen, int(sys.argv[2]), int(sys.argv[2]), addprortycaculator)
        #             print("win",count)
                    recs=[[int(recl[5])]+recl];curchr=recl[4]
                    currentchrLen=int(recl[5])
                else:#first line
        #             print(recl)
                    recs=[[int(recl[5])]+recl];curchr=recl[4];currentchrLen=int(recl[5])
            else:
                win.slidWindowOverlap(recs, currentchrLen, int(sys.argv[2]), int(sys.argv[2]), addprortycaculator)
            ofo.close()
            f.close();addprortycaculator.temp.close()
            #NF==? print $1~$(NF-1),"",$NF need to be variable for some col no best_strand rec , fuck affy
            os.system("""awk 'BEGIN{FS="\t";OFS="\t"}{if(NF==41){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,"",$41}else{print $0}}' KPS_Rice_scored.txt.sortedredup6000 > KPS_Rice_scored.txt.sortedredup6000.modifylastcol""")
        #     os.system("""awk 'BEGIN{FS="\t"}{if($34!=8 && $34!=9){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol|awk 'BEGIN{FS="\t"}{if($32=="neutral"){print $0}}' > Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol1_7netural""")
        #     os.system("""awk 'BEGIN{FS="\t"}{if($34==8){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol|cat - Axiom_KPSmilet_redo_scored.txt.sorted3700withtitle.miodifylastcol1_7netural|sort -t$'\t'  -k5,5 -k6,6n > Axiom_KPSmilet_redo_scored.txt.sorted3700.sorted8M1_7netural""")
        elif sys.argv[3].strip().upper()=="I":
            #sort
#             os.system("sed -n '1,16p' " +sys.argv[1]+" > titlelinetemp")
    #         os.system("sed -n '17,$p' " +sys.argv[1]+"|sort -t$','  -k4,4n -k5,5n - |cat titlelinetemp -  >"+sys.argv[1]+".sorted");
#             sys.argv[1]=sys.argv[1]+".sorted"
    
            win = Util.Window()
            ofo=open(sys.argv[1]+sys.argv[2],'w')
            f=open(sys.argv[1],'r')
            
            recs=[];count=0;tcount=0
            curchr=None
            for line in f:
                tcount+=1
                recl=re.split(r",",line.strip())
    #             print(recl)
                if recl[3]!="Chromosome" and curchr is None:
                    continue
                elif curchr is None:
                    curchr=recl[3];title=recl#"Chromosome"
                    tidx=title.index("Source_Version");print("curchr:",curchr,"tilingorder idx:",tidx);best_recommendation=title.index("Final_Score")
                    addprortycaculator = Caculators.Caculator_selectsitesForillumina(of=ofo,tilingorderidx=tidx,best_recommendation=best_recommendation,T=float(sys.argv[4]),VIPidx=title.index("Source"))
                    addprortycaculator.spf=open(sys.argv[1]+sys.argv[2]+".originscored",'w')
                    continue
                if curchr==recl[3]:
                    currentchrLen=int(recl[4])
                    recs.append([int(recl[4])]+recl)
                elif curchr!="Chromosome":
                    if curchr=="0":
                        for e in recs:
                            if e[title.index("Source")+1]=="VIP" or( int(e[tidx+1])<=3 and float(e[best_recommendation+1])>0.6): 
                                print(*e[1:],sep="\t",file=addprortycaculator.spf)
                                print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],addprortycaculator.m[e[8].upper()],"Soybean","FALSE",sep=",",file=ofo)
                    else:
                        win.slidWindowOverlap(recs, currentchrLen, int(sys.argv[2]), int(sys.argv[2]), addprortycaculator)
                    recs=[[int(recl[4])]+recl];curchr=recl[3]
                    currentchrLen=int(recl[4])
                else:#first line
                    recs=[[int(recl[4])]+recl];curchr=recl[3];currentchrLen=int(recl[4])
            else:
                win.slidWindowOverlap(recs, currentchrLen, int(sys.argv[2]), int(sys.argv[2]), addprortycaculator)
    
            #slide window
        exit()        
    ###other firt function##############################################
    if len(sys.argv)<4:
        print("python chip_design.py chrchangemapfile blastout.pos out.pos")
        exit(-1)
    f=open(sys.argv[1],'r')
    changemap={}
    for line in f:
        linelist=re.split(r"\s+",line.strip())
        changemap[linelist[0]]=linelist[1].strip()
    f.close()
    f=open(sys.argv[2],'r')
    fo=open(sys.argv[3],'w')
    for line in f:
        linelist=re.split(r"\s+",line.strip())
        print(changemap[linelist[1].strip()],linelist[2].strip(),linelist[0].strip(),sep="\t",file=fo)
    fo.close();f.close()