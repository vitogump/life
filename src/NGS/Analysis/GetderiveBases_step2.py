import re, sys, time, pickle, os
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm
'''
Created on 2013-9-24

@author: liurui
'''
if len(sys.argv) < 7:
    print("python GetderiveBases_step2.py [duckref] [originalspeciesref] [winFile1] [winFile2]... [winFileN] [chromtable] [tempwinDBName] [percentage] [outfilename] [m/l]")
    exit(-1)
duckref = sys.argv[1]
originalspeciesref = sys.argv[2]
winFileName6Fields = sys.argv[3:-5]
chromtable = sys.argv[-5]
tempwinDBName = sys.argv[-4]
percentage = float(sys.argv[-3])
outfilename = sys.argv[-2]
morethan_lessthan = sys.argv[-1].strip()
trscptable = "transcript"

snptables = ["yingtaogusnp", "fanyasnp", "gaoyousnp", "jindingsnp", "kangbeiersnp", "lianchengbaisnp", "shanmasnp"]
fintableNamelist = re.split(r"\.", outfilename)
blastoutName = outfilename.replace("fa", "blast")

finaltable = "_".join(fintableNamelist) + "_allselectedSNP"
#gene_sample_venn="gene_sample_venn"
vcftable = None
duckreffile = open(duckref, 'r')
originspeciesfile = open(originalspeciesref, 'r')
#outfile = open(outfilename, 'w')
#filepos = int(outfile.tell())
selectedWins = {}
if __name__ == '__main__':
    try:
        duckrefindex = pickle.load(open(duckref + ".myindex", 'rb'))
        originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
    except IOError:
        Util.generateIndexByChrom(duckref, duckref + ".myindex")
        Util.generateIndexByChrom(originalspeciesref, originalspeciesref + ".myindex")
        duckrefindex = pickle.load(open(duckref + ".myindex", 'rb'))
        originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
    
    dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
    dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", finaltable, "chicken", "tinytext", "default null"))
    a = os.popen("awk '$1!~/^#/ && $5==1 && $4>30 && $6==0 {print $0}' " + blastoutName)
#    hits=a.readlines()

    lastbasesAccur = {}
    revcom=False
#    initial
    hit = a.readline()
    hitlist = re.split(r"\s+", hit)

    sendpos = int(hitlist[9])
    sstartpos = int(hitlist[8])
    qstartpos = int(hitlist[6])
    snpindex = 26 - qstartpos
    if sstartpos > sendpos:
        temp = sstartpos
        sstartpos = sendpos
        sendpos = temp
        revcom=True
    lastsnpID = hitlist[0]
    chrom= hitlist[1]
    RefSeqMap = Util.getRefSeqBypos(refFastahander=originspeciesfile, refindex=originalspeciesindex, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
    if revcom:
        tempStr=RefSeqMap[chrom][1:]
        tempStr.reverse()
        RefSeqMap[chrom][1:]=Util.complementary(tempStr)
        revcom=False
        
    lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
    for hit in a:
        print(hit)
        hitlist = re.split(r"\s+", hit)
        chrom = hitlist[1]
        sstartpos = int(hitlist[8])
        sendpos = int(hitlist[9])
        qstartpos = int(hitlist[6])
        snpindex = 26 - qstartpos
        if sstartpos > sendpos:
            temp = sstartpos
            sstartpos = sendpos
            sendpos = temp
            revcom=True
        if lastsnpID == hitlist[0]:
            RefSeqMap = Util.getRefSeqBypos(refFastahander=originspeciesfile, refindex=originalspeciesindex, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
            if revcom:
                tempStr=RefSeqMap[chrom][1:]
                tempStr.reverse()
                RefSeqMap[chrom][1:]=Util.complementary(tempStr)
                revcom=False            
            print(lastsnpID,"".join(RefSeqMap[chrom][1:]),file=open("chickensnpflank.fa",'a'))
            if RefSeqMap[chrom][snpindex + 1] in lastbasesAccur:
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]].append((chrom, sstartpos, sendpos))
            else:
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
        else:
            if len(lastbasesAccur.keys()) == 1:
                for bases in lastbasesAccur:#only once
                    dbtools.operateDB("update", "update " + finaltable + " set chicken='" + bases + "' where snpID='" + lastsnpID + "'")
            elif len(lastbasesAccur.keys()) == 0:
                exit(-1)
            RefSeqMap = Util.getRefSeqBypos(refFastahander=originspeciesfile, refindex=originalspeciesindex, currentChromNO=chrom, startpos=sstartpos, endpos=sendpos)
            if revcom:
                tempStr=RefSeqMap[chrom][1:]
                tempStr.reverse()
                RefSeqMap[chrom][1:]=Util.complementary(tempStr)
                revcom=False            
            print(hitlist[0],"".join(RefSeqMap[chrom][1:]),file=open("chickensnpflank.fa",'a'))
#            dbtools.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][snpindex + 1] + "' where snpID='" + hitlist[0] + "'")
            lastsnpID = hitlist[0]
            
            lastbasesAccur.clear()
            if RefSeqMap[chrom][snpindex + 1] in lastbasesAccur:
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]].append((chrom, sstartpos, sendpos))
            else:
                lastbasesAccur[RefSeqMap[chrom][snpindex + 1]] = [(chrom, sstartpos, sendpos)]
                          
    dbtools.disconnect()
    originspeciesfile.close()
    a.close()

    
    
    
    
    
        
