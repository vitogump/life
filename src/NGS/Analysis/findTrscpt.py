import re, sys,time
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm
SLEEP_FOR_NEXT_TRY=3
'''
Created on 2013-9-3

@author: liurui
'''
if len(sys.argv) != 6:
    print("python findTrscpt.py [winFile1] [tempwinDBName] [threshold] [outfilename] [m/l]")
    exit(-1)
winFileName6Field = sys.argv[1]
tempwinDBName = sys.argv[2]
threshold = sys.argv[3]
outfilename=sys.argv[4]
morethan_lessthan=sys.argv[5].strip()
trscptable="transcript"
gene_sample_venn="gene_sample_venn"
vcftable=None
outfile=open(outfilename,'w')
if __name__ == '__main__':
    dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
    dbtools.operateDB("alter","alter table "+gene_sample_venn+" add "+outfilename+" smallint(3) default 0") 
    winGenome = Util.WinInGenome(tempwinDBName, winFileName6Field)
    time.sleep(SLEEP_FOR_NEXT_TRY)
    if morethan_lessthan=="m" or morethan_lessthan=="M":
        selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue>=" + threshold)
    elif morethan_lessthan=="l" or morethan_lessthan=="L":
        selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue<=" + threshold)

    for win in selectedWins:
        winRegion=(win[0],win[1],40000,20000,win[5])
        winGenome.collectTrscptInWin(dbtools, trscptable, vcftable, winRegion)
    for win in sorted(winGenome.winContainTrscptMap.keys()):
        print("\t".join(map(str,win)),winGenome.winContainTrscptMap[win],sep="\t",file=outfile)
        for gene in winGenome.winContainTrscptMap[win]:
            print(gene)
            print("update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
            dbtools.operateDB("update","update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
    outfile.close()
#    winGenome.windbtools.drop_table(winGenome.wintable)
#    winGenome.windbtools.disconnect()
