import re, sys,time
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm
from optparse import OptionParser
SLEEP_FOR_NEXT_TRY=3
'''
Created on 2013-9-3

@author: liurui
'''
parser = OptionParser()
parser.add_option("-i", "--winfile", dest="winfileName",
                  help="reference.fa", metavar="FILE")
parser.add_option("-D", "--tempDBname", dest="tempdbname", help="dbname")
parser.add_option("-t", "--threshold", dest="threshold", help="conflict with -p")
parser.add_option("-p", "--percentage", dest="percentage", help="conflict with -t")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-x", "--morethan_lessthan", dest="morethan_lessthan", help="morethan or lessthan")
parser.add_option("-T", "--trscptable", dest="trscptable", help="trscptable")
parser.add_option("-u", "--upextend", dest="upextend", help="upextend")
parser.add_option("-d", "--downextend", dest="downextend", help="downextend")
parser.add_option("-s","--slideSize",dest="slideSize",default="20000",help="win slide size")
parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
#if len(sys.argv) != 6:
#    print("python findTrscpt.py [winFile1] [tempwinDBName] [threshold] [outfilename] [m/l]")
#    exit(-1)
upextend=int(options.upextend);slideSize=int(options.slideSize);winWidth=int(options.winWidth)
downextend=int(options.downextend)
winFileName6Field = options.winfileName
path=re.search(r'^.*/',options.winfileName).group(0)
tempwinDBName = options.tempdbname
threshold = options.threshold
percentage = options.percentage
outfilename=path+options.outfileprename
morethan_lessthan=options.morethan_lessthan
TranscriptGenetable=options.trscptable
if percentage!=None and threshold!=None:
    print("-t conflict with -p")
    exit(-1)
#gene_sample_venn="gene_sample_venn"
vcftable=None
outfile=open(outfilename,'w')
outfileNameWINwithGENE=winFileName6Field+".wincopywithgene"
if __name__ == '__main__':
    dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
#    dbtools.operateDB("alter","alter table "+gene_sample_venn+" add "+outfilename+" smallint(3) default 0") 
    winGenome = Util.WinInGenome(tempwinDBName, winFileName6Field)
    time.sleep(SLEEP_FOR_NEXT_TRY)
    
    selectWinNos="threshold method"
    if percentage!=None:
        totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintable)[0][0]
        selectWinNos = int(float(percentage) * totalWin)
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue != 'NA' order by zvalue desc limit 0," + str(selectWinNos))
            print("select * from "+winGenome.wintable + " where zvalue != 'NA' order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue != 'NA' order by zvalue asc limit 0," + str(selectWinNos))
            print("select * from " + winGenome.wintable + " where zvalue != 'NA' order by zvalue asc limit 0," + str(selectWinNos))
    elif threshold!=None:
        if morethan_lessthan=="m" or morethan_lessthan=="M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue!= 'NA' and zvalue>=" + threshold)
        elif morethan_lessthan=="l" or morethan_lessthan=="L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintable + " where zvalue!= 'NA' and zvalue<=" + threshold)
        print("select * from " + winGenome.wintable + " where zvalue!= 'NA' and zvalue>=" + threshold)
    selectedWins.sort(key=lambda listRec:float(listRec[5]))
    print(outfilename,selectWinNos,selectedWins[0],selectedWins[-1],len(selectedWins))
    selectedWinMap={}
    for win in selectedWins:
        if win[0] in selectedWinMap:
            selectedWinMap[win[0]].append(win)
        else:
            selectedWinMap[win[0]]=[win]
    #selectedWins is a list :[('KB743038.1', '9', '181586', '219606', '0.3816832053195056', '-0.00013080719016'),(),(),(),(),()]
    #selectedWinMap {chrom1:[(chrom1, '9', '181586', '219606', '0.3816832053195056', '-0.00013080719016'),(),(),()],chrom2:[],}
    #mergedRegion [(chrom1, '9', '181586', '219606', '0.3816832053195056', '-0.00013080719016'),(),(),()] continues
    selectedRegion={}
    #fill selectedRegion map
    #selectedRegion {chrom:[chrom,Region_start,Region_end,Nwin,extremeValue],chrom:[],,,,}
    #merge continues win into a region
    for chrom in selectedWinMap:
        selectedWinMap[chrom].sort(key=lambda listRec: int(listRec[1]))
        selectedRegion[chrom]=[]
        mergedRegion=[selectedWinMap[chrom][0]]
        for i in range(1,len(selectedWinMap[chrom])):
            if int(selectedWinMap[chrom][i-1][1])+1==int(selectedWinMap[chrom][i][1]):#continues win
                mergedRegion.append(selectedWinMap[chrom][i])
            else:#not continues
                #process last region
                Region_start=int(mergedRegion[0][1])*slideSize-upextend
                Region_end=int(mergedRegion[-1][1])*slideSize+winWidth+downextend
                Nwin=len(mergedRegion)
                extremeValues=[]
                for e in mergedRegion:
                    extremeValues.append(float(e[5]))
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue=max(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=min(extremeValues)
                selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
                #process this win
                mergedRegion=[selectedWinMap[chrom][i]]
        else:
            Region_start=int(mergedRegion[0][1])*slideSize-upextend
            Region_end=int(mergedRegion[-1][1])*slideSize+winWidth+downextend
            Nwin=len(mergedRegion)
            extremeValues=[]
            for e in mergedRegion:
                extremeValues.append(float(e[5]))
            if morethan_lessthan == "m" or morethan_lessthan == "M":
                extremeValue=max(extremeValues)
            elif morethan_lessthan == "l" or morethan_lessthan == "L":
                extremeValue=min(extremeValues)            
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
#    get final table
    final_table={}
    for chrom in selectedRegion:
        for region in selectedRegion[chrom]:
            final_table[region]=winGenome.collectTrscptInWin(dbtools,TranscriptGenetable,vcftable,region)
#    for win in selectedWins:
#        winRegion=(win,upextend,downextend)
#        winGenome.collectTrscptInWin(dbtools, TranscriptGenetable, vcftable, winRegion)
    for region in sorted(final_table.keys()):
        tcpts=""
        for tcpt in final_table[region]:
            tcpts+=(tcpt[0]+"\t")
        print("\t".join(map(str,region)),tcpts,sep="\t",file=outfile)
    winGenome.appendGeneName(TranscriptGenetable, dbtools, winWidth, slideSize, outfileNameWINwithGENE)   
    winGenome.windbtools.drop_table(winGenome.wintable)        
#        for gene in winGenome.winContainTrscptMap[win]:
#            print(gene)
#            print("update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
#            dbtools.operateDB("update","update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
    outfile.close()
#    winGenome.windbtools.drop_table(winGenome.wintable)
#    winGenome.windbtools.disconnect()
