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
    ###add priority to scored file
    if len(sys.argv)!=3:
        print("python chip_design.py scoredfile winsize")
    
    os.system("sort -t$'\t'  -k5,5 -k6,6n "+sys.argv[1]+">"+sys.argv[1]+".sorted")
    f=open(sys.argv[1],'r');title=re.split(r"\t",f.readline().strip())
    tidx=title.index("tiling_order");curchr=title[4]
    win = Util.Window()
    hp_caculator = Caculators.Caculate_Hp(SeqMethodlist=methodlist,minsnps=10,depth=int(options.mindepth))
    for line in f:
        recl=re.split(r"\t",line.strip())
        recs=[]
        if curchr==recl[4]:
            currentchrLen=int(recl[5])
            print("collect rec in a win")
            recs.append(recl)
        elif recl[4]!="cust_chr":
            print("sliding win")
            win.slidWindowOverlap(recs, currentchrLen, 800, 800, hp_caculator)
            
    ###
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