'''
Created on 2018年5月24日

@author: Dr.liu
'''
from multiprocessing.dummy import Pool
import os,re,sys
SNPIDmissedinextract=[]
def extarctBlastOutvariablelenwaper(a):
    extarctBlastOutvariablelen(**a)
def extarctBlastOutvariablelen(BlastOutFile,faMap,fromposOfBlastf="missedinextracttemp2"):
    global SNPIDmissedinextract
    posfile=open(BlastOutFile+".pos",'a')
#     me=open(missedinextract,'r')
    ospid=os.getpid()
    f=open(BlastOutFile+"missedinextract"+str(ospid)+".fa",'a')
    bo=open(BlastOutFile,"r")
    bo.seek(fromposOfBlastf[0])
    for outline in bo:
        for snpid in SNPIDmissedinextract:
            if snpid in outline.strip() and outline[:len(snpid.strip())]==snpid.strip():
                hitlist = re.split(r"\s+", outline)
                if int(hitlist[3])>=len(faMap[IDline.strip()])-10 and int(hitlist[4])<=3 and lastsnpID!=hitlist[0]:
                    sendpos = int(hitlist[9])
                    sstartpos = int(hitlist[8])
                    qstartpos = int(hitlist[6])
                    blastlen=int(hitlist[3])
                    flanklen=faMap[IDline.strip()].index("N")
                    snp_loc_s=sstartpos+flanklen
                    snpindex = flanklen+1 - qstartpos
                    if sstartpos > sendpos:
                        temp = sstartpos
                        sstartpos = sendpos
                        sendpos = temp
                        revcom="-"
                    lastsnpID = hitlist[0]
                    chrom= hitlist[1]
                    hitlist = re.split(r"\s+", a.readline().strip())
                    print(lastsnpID,chrom,snp_loc_s,file=posfile)
                    break
                else:#same snpID means not the first record, or not blast records (just #query)
                    print(">"+IDline+faMap[IDline.strip()],file=f)
                    sys.stdout.flush()
            elif snpid in outline.strip():#
                print(">"+IDline+faMap[IDline.strip()],file=f)
                break#
        if IDline.strip() in 
           
    for IDline in me:
        
            
                
                
        a=os.popen("grep "+IDline.strip()+" "+BlastOutFile+"|awk '$1!~/^#/ && $4>="+str(len(faMap[IDline.strip()])-10)+" && $5<=3  {print $0}' " )
        hit = a.readline()
        if hit.strip()=="":
            
            continue
        hitlist = re.split(r"\s+", hit)
        sendpos = int(hitlist[9])
        sstartpos = int(hitlist[8])
        qstartpos = int(hitlist[6])
        blastlen=int(hitlist[3])
        flanklen=faMap[IDline.strip()].index("N")
        snp_loc_s=sstartpos+flanklen
        snpindex = flanklen+1 - qstartpos
        if sstartpos > sendpos:
            temp = sstartpos
            sstartpos = sendpos
            sendpos = temp
            revcom="-"
        lastsnpID = hitlist[0]
        chrom= hitlist[1]
        hitlist = re.split(r"\s+", a.readline().strip())
        if lastsnpID!=hitlist[0]:
            print(lastsnpID,chrom,snp_loc_s,file=posfile)
        else:
            print(">"+IDline+faMap[IDline.strip()],file=f)
        sys.stdout.flush()
    posfile.close()
    me.close()
    f.close()
def extarctBlastOut(BlastOutFile,queryFaFile,flanklen=60,mismatch=6):
    print(" query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
    ambigousfile=open(BlastOutFile+".ambigous.fa",'w')
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
    a = os.popen("awk '$1!~/^#/ && $4>="+str(2*flanklen+1)+" && $5<="+str(mismatch)+"  {print $0}' " + BlastOutFile)
    posfile=open(BlastOutFile+".pos",'w')

    revcom="+"
#    initial
    hit = a.readline()
    hitlist = re.split(r"\s+", hit)

    sendpos = int(hitlist[9])
    sstartpos = int(hitlist[8])
    qstartpos = int(hitlist[6])
    blastlen=int(hitlist[3])
    
    snpindex = flanklen+1 - qstartpos
    if sstartpos > sendpos:
        temp = sstartpos
        sstartpos = sendpos
        sendpos = temp
        revcom="-"
    snp_loc_s=sstartpos+flanklen
    lastsnpID = hitlist[0]
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
            snp_loc_s=sstartpos+flanklen
            snpindex = flanklen+1 - qstartpos
            if sstartpos > sendpos:
                temp = sstartpos
                sstartpos = sendpos
                sendpos = temp
                revcom="-"
            else:
                revcom="+"
            lastsnpID = hitlist[0]
    posfile.close()
    print("postion writing done")
    
    os.system("awk '{print $1}' "+BlastOutFile+".pos > extractedexactly"+str(subid))
    
    os.system("grep -vFf extractedexactly"+str(subid)+" allSNPID"+str(subid)+" > missedinextractmp"+str(subid))
    
#     a = os.popen("less missedinextractmp"+str(subid)+"|wc -l ")
#     totalsnpforcount=int(a.readline().strip())
#     a.close()
    parameterstuples_list=[]
    blastoutfilendpos=f.seek(0,os.SEEK_END)
    d=int(blastoutfilendpos/36)+1
    j=0
    for i in range(0,blastoutfilendpos,d):
#         print("sed -n '"+str(i)+","+str(i+d-1)+"p' missedinextractmp"+str(subid)+" >missedinextractmp"+str(subid)+str(j))
#         os.system("sed -n '"+str(i)+","+str(i+d-1)+"p' missedinextractmp"+str(subid)+" >missedinextractmp"+str(subid)+str(j))
        parameterstuples_list.append({"BlastOutFile":BlastOutFile,"faMap":faMap,"fromposOfBlastf":i})
        j+=1
    pool=Pool(36)
    pool.map(extarctBlastOutvariablelenwaper,parameterstuples_list)
    pool.close()
    pool.join()
    ambigousfile.close()
    print("finish")
if __name__ == '__main__':
    if len(sys.argv)<3:
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