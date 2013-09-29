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
    a = os.popen("awk '$1!~/^#/ && $5==1 && $4==51 && $6==0 {print $0}' " + blastoutName)
#    hits=a.readlines()
    lastsnpID = None
    basesAccur = {}
    for hit in a:
        print(hit)
        hitlist = re.split(r"\s+", hit)
        chrom = hitlist[1]
        startpos = int(hitlist[8])
        endpos = int(hitlist[9])
        if startpos>endpos:
            temp=startpos
            startpos=endpos
            endpos=temp
        if lastsnpID == hitlist[0]:
            RefSeqMap, lastchromNo = Util.getRefSeqBypos(refFastahander=originspeciesfile,refindex=originalspeciesindex,currentChromNO= chrom, startpos=startpos,endpos= endpos)
            if RefSeqMap[chrom][25] in basesAccur:
                basesAccur[RefSeqMap[chrom][25]].append((chrom, startpos, endpos))
            else:
                basesAccur[RefSeqMap[chrom][25]] = [(chrom, startpos, endpos)]
        else:
            if len(basesAccur.keys()) == 1:
                for bases in basesAccur:#only once
                    dbtools.operateDB("update", "update " + finaltable + " set chicken='" + bases + "' where snpID='" + lastsnpID + "'")
            elif len(basesAccur.keys()) == 0:
                RefSeqMap, lastchromNo = Util.getRefSeqBypos(refFastahander=originspeciesfile,refindex=originalspeciesindex,currentChromNO= chrom, startpos=startpos,endpos= endpos)
                dbtools.operateDB("update", "update " + finaltable + " set chicken='" + RefSeqMap[chrom][25] + "' where snpID='" + hitlist[0] + "'")
            lastsnpID = hitlist[0]
            basesAccur.clear()
    dbtools.disconnect()
    originspeciesfile.close()
    a.close()

    
    
    
    
    
        
