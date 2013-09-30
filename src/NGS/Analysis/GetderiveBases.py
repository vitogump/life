import re, sys, time, pickle
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm
SLEEP_FOR_NEXT_TRY = 3
'''
Created on 2013-9-10

@author: liurui
'''
winwidth = 40000
slidesize = 20000

if len(sys.argv) < 7:
    print("python GetderiveBases.py [duckref] [originalspeciesref] [winFile1] [winFile2]... [winFileN] [chromtable] [tempwinDBName] [percentage] [outfilename] [m/l]")
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
#snptables=[]
#for winFileName6Field in winFileName6Fields:
#    pureName = re.search(r"[^/]*$", winFileName6Field).group(0)  # for linux
#    pureName = re.split(r'\.',pureName)[0]
#    snptables.append(pureName)
snptables = ["yingtaogusnp", "fanyasnp", "gaoyousnp", "jindingsnp", "kangbeiersnp", "lianchengbaisnp", "shanmasnp"]
fintableNamelist = re.split(r"\.", outfilename)

finaltable = "_".join(fintableNamelist) + "_allselectedSNP"
#gene_sample_venn="gene_sample_venn"
vcftable = None
duckreffile = open(duckref, 'r')
originspeciesfile = open(originalspeciesref, 'r')
outfile = open(outfilename, 'w')
filepos = int(outfile.tell())
selectedWins = {}
testName=Util.random_str()
testfile=open(testName+"testsnpfile.txt",'w')
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
    TABLES = {}
    TABLES[finaltable] = (
        "CREATE TABLE " + finaltable + " ("
        " `snpID` varchar(128) NOT NULL,"
        " `chrID` varchar(128) NOT NULL DEFAULT '',"
        " `snp_start_pos` bigint(20) NOT NULL DEFAULT '0',"
        " `comefrom` text DEFAULT NULL,"
        " `ref_bases` tinytext,"
        " PRIMARY KEY (`snpID`) "
        ")ENGINE=InnoDB DEFAULT CHARSET=utf8"
        )
    dbtools.drop_table(finaltable)
    dbtools.create_table(TABLES)
    #add column for every tables
    for snptable in snptables:
        dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", finaltable, snptable, "varchar(128)", "default null"))
    dbtools.operateDB("callproc", "mysql_sp_add_column", data=("life_pilot", finaltable, "fafilepos", "bigint(128)", "default null"))
    for winFileName6Field in winFileName6Fields:
        winGenome = Util.WinInGenome(tempwinDBName, winFileName6Field)
        pureName = re.search(r"[^/]*$", winFileName6Field).group(0)  # for linux
        time.sleep(SLEEP_FOR_NEXT_TRY)
        totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintable)[0][0]
        selectWinNos = int(percentage * totalWin)
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins[winFileName6Field] = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue != 'NA' order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins[winFileName6Field] = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue != 'NA' order by zvalue asc limit 0," + str(selectWinNos))
        print(str(selectWinNos), selectedWins[winFileName6Field][-1], winFileName6Field)
        print("selecting " + winFileName6Field + " wins is finished")
        for win in selectedWins[winFileName6Field]:
            for processedWinFile in selectedWins.keys():
                if winFileName6Field == processedWinFile:
                    continue
                for processedWin in selectedWins[processedWinFile]:
                    if win[0] == processedWin[0] and win[1] == processedWin[1]:
                        snps = dbtools.operateDB("update", "update " + finaltable + " set comefrom=concat(comefrom,' " + pureName + "') where  snp_start_pos>=" + str(int(win[1]) * slidesize) + " and snp_start_pos<=" + str(int(win[1]) * slidesize + winwidth))
                        continue
#            winRegion=(win[0],win[1],winwidth,slidesize,win[5])
            for snptable in snptables:
#                print("select","select * from " +snptable+" where chrID="+win[0]+" and snp_start_pos>="+str(int(win[1])*slidesize)+" and snp_start_pos<="+str(int(win[1])*slidesize+winwidth))
                snps = dbtools.operateDB("select", "select * from " + snptable + " where chrID='" + win[0] + "' and snp_start_pos>=" + str(int(win[1]) * slidesize) + " and snp_start_pos<=" + str(int(win[1]) * slidesize + winwidth))
                for snp in snps:
                    dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", snp[5])
                    alt_dp4 = snp[4] + ":" + dp4.group(1) + "," + dp4.group(2) + "," + dp4.group(3) + "," + dp4.group(4)
#                    print("insert into "+finaltable+"(snpID,chrID,snp_start_pos,comefrom,ref_bases,"+snptable+") values(%s,%s,%s,%s,%s,%s) on duplicate key update comefrom=concat(comefrom,"+"' "+snptable+"'),"+snptable+"= '"+alt_dp4+"'",(snp[0],snp[1],snp[2],snptable,snp[3],alt_dp4))
                    # comefrom should not be null in mysql database
                    dbtools.operateDB("insert", "insert into " + finaltable + "(snpID,chrID,snp_start_pos,comefrom,ref_bases," + snptable + ") values(%s,%s,%s,%s,%s,%s) on duplicate key update " + snptable + "= '" + alt_dp4 + "'", data=(snp[0], snp[1], snp[2], pureName, snp[3], alt_dp4))
        print("filled all snp in the wins of " + winFileName6Field + " into the mysql table " + finaltable)
    print("all snps were insert into " + finaltable + " already . \n going to extract snp flanks seqs")
    
    totalChroms = dbtools.operateDB("select", "select count(*) from " + chromtable)[0][0]
    for i in range(0, totalChroms, 20):
        currentsql = "select * from " + chromtable + " order by chrID limit " + str(i) + ",20"
        result = dbtools.operateDB("select", currentsql)

        for row in result:
            currentchrID = row[0]
            currentchrLen = int(row[2])
            totalsnpsInchr = dbtools.operateDB("select", "select count(*) from " + finaltable + " where chrID ='" + currentchrID + "'")[0][0]
            duckreffile.seek(duckrefindex[currentchrID])
            RefSeqMap, lastchromNo = Util.getRefSeqMap(refFastafilehander=duckreffile, currentChromNO=currentchrID)
            for j in range(0, totalsnpsInchr, 1000):
                print("select * from " + finaltable + " where chrID='" + currentchrID + "' order by snp_start_pos limit " + str(j) + ",1000")
                snps = dbtools.operateDB("select", "select * from " + finaltable + " where chrID='" + currentchrID + "' order by snp_start_pos limit " + str(j) + ",1000")
#                int(snps[-1][2])
                if snps == None:
                    print("no snp in " + currentchrID + " in table " + finaltable)
                    continue
                for snp in snps:
                    currentsnpID = snp[0]
                    currentsnpChrId = snp[1]
                    currentsnpPos = int(snp[2])
                    if len(snp[4]) != 1:
#                        print(snp[4])
                        continue# skip indel
                    if currentsnpPos + 25 <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos - 25 > RefSeqMap[lastchromNo][0] :
                        snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - 25 - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                        print(currentsnpID,snpflankseq[25],file=testfile)
                        snpflankseq=snpflankseq[0:25]+'N'+snpflankseq[26:]
                        
                    elif currentsnpPos <= RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1 and currentsnpPos + 25 > RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1:
                        snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - 25 - RefSeqMap[currentsnpChrId][0]):(currentsnpPos - RefSeqMap[currentsnpChrId][0] + 1)])
                        print(currentsnpID,snpflankseq[25],file=testfile)
                        snpflankseq=snpflankseq[0:25]+'N'
                        
                    elif currentsnpPos - 25 <= RefSeqMap[lastchromNo][0]:
                        snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                        print(currentsnpID,snpflankseq[0],file=testfile)
                        snpflankseq = 'N'+snpflankseq[1:26]
                        
                    elif currentsnpPos>RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1:
                        RefSeqMap, lastchromNo = Util.getRefSeqMap(duckreffile, currentChromNO=currentchrID, preBaseTotal=RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1)
                        snpflankseq = ''.join(RefSeqMap[currentsnpChrId][(currentsnpPos - RefSeqMap[currentsnpChrId][0]):(currentsnpPos + 25 - RefSeqMap[currentsnpChrId][0] + 1)])
                        print(currentsnpID,snpflankseq[0],file=testfile)
                        snpflankseq = 'N'+snpflankseq[1:26]
                        
                    else:
                        print("what's wrong?")
                    print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=outfile)
#                    print("update "+finaltable+" set fafilepos="+str(filepos)+" where snpID='"+currentsnpID+"'")
                    dbtools.operateDB("update", "update " + finaltable + " set fafilepos=" + str(filepos) + " where snpID='" + currentsnpID + "'")
                    filepos = int(outfile.tell())
    dbtools.disconnect()               
    originspeciesfile.close()                
    duckreffile.close()
    outfile.close()
