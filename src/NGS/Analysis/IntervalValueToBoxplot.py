# -*- coding: UTF-8 -*-
'''
Created on 2014-7-13

@author: liurui
'''
from optparse import OptionParser
import re

from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-i", "--interval", dest="interval", nargs=3,
                  help="minvalue maxvalue breaks", metavar="FILE")
parser.add_option("-K", "--dndsfile", dest="dndsfile", help="dndsfile,dnds_col_idx(begin as 1)", nargs=2)
parser.add_option("-D", "--tempDBname", dest="tempdbname", help="dbname")
parser.add_option("-x", "--morethan_lessthan", dest="morethan_lessthan", default="l", help="morethan or lessthan")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")

parser.add_option("-T", "--trscptable", dest="trscptable", help="trscptable")
parser.add_option("-u", "--upextend", dest="upextend", help="upextend")
parser.add_option("-d", "--downextend", dest="downextend", help="downextend")
parser.add_option("-s", "--slideSize", dest="slideSize", default="20000", help="win slide size")
parser.add_option("-b", "--winfileName", dest="winfileName", help="winfileName ")
parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
parser.add_option("-X", "--winType", dest="winType", default="zvalue", help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()

TranscriptGenetable = options.trscptable
print(options.interval, options.dndsfile, options.outfileprename,)
tempwinDBName = options.tempdbname
winFileName6Field = options.winfileName
print("sssssssssss",options.winfileName)
path=re.search(r'^.*/',options.winfileName).group(0)
outfilename=path+options.outfileprename

minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
dincrease = (maxvalue - minvalue) / breaks
intervalmap = {}  # {(minvalue,maxvalue):selectedRegion,(minvalue,maxvalue):selectedRegion,,,,,,}
upextend = int(options.upextend);slideSize = int(options.slideSize);winWidth = int(options.winWidth);winsizebykb=int(winWidth/1000)
downextend = int(options.downextend)
morethan_lessthan = options.morethan_lessthan


dndsfile = open(options.dndsfile[0], 'r')
outfile=open(outfilename,'w')
dndscolidx = int(options.dndsfile[1]) - 1
dndstitle = dndsfile.readline()
if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    winGenome = Util.WinInGenome(tempwinDBName, winFileName6Field)
    while minvalue <= maxvalue - dincrease:
        intervalmap[minvalue,minvalue + dincrease]=[]
        selectedWinMap = {}  # selectedWinMap {chrom1:[(chrom1, '9', '181586', '219606', '0.3816832053195056', '-0.00013080719016'),(),(),()],chrom2:[],} derived from findTrscpt.py
        print(minvalue, minvalue + dincrease)
        selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where " + options.winType + "!= 'NA' and " + options.winType + "/"+str(winsizebykb)+">=" + str(minvalue) + " and " + options.winType+ "/"+str(winsizebykb) + "<" + str(minvalue + dincrease))
        selectedWins.sort(key=lambda listRec:float(listRec[5]))
        for win in selectedWins:
            if win[0] in selectedWinMap:
                selectedWinMap[win[0]].append(win)
            else:
                selectedWinMap[win[0]] = [win]
        # selectedRegion {chrom:[chrom,Region_start,Region_end,Nwin,extremeValue],chrom:[],,,,}
        # mergedRegion [(chrom1, '9', '181586', '219606', '0.3816832053195056', '-0.00013080719016'),(),(),()] continues
        selectedRegion = {}
        for chrom in selectedWinMap:
            selectedWinMap[chrom].sort(key=lambda listRec: int(listRec[1]))
            selectedRegion[chrom] = []
            mergedRegion = [selectedWinMap[chrom][0]]
            i = 1
            while i < len(selectedWinMap[chrom]):
    #             print(chrom,selectedWinMap[chrom][i])
    #             try:
                if int(selectedWinMap[chrom][i - 1][1]) + 1 == int(selectedWinMap[chrom][i][1]):  # continues win
                    mergedRegion.append(selectedWinMap[chrom][i])
                else:  # not continues
                    # process last region
                    Region_start = int(mergedRegion[0][1]) * slideSize - upextend
                    Region_end = int(mergedRegion[-1][1]) * slideSize + winWidth + downextend
                    Nwin = len(mergedRegion)
                    extremeValues = []
                    for e in mergedRegion:
                        if options.winType == "winvalue":
                            extremeValues.append(float(e[4]))
                        elif options.winType == "zvalue": 
                            extremeValues.append(float(e[5]))
                    if morethan_lessthan == "m" or morethan_lessthan == "M":
                        extremeValue = max(extremeValues)
                    elif morethan_lessthan == "l" or morethan_lessthan == "L":
                        extremeValue = min(extremeValues)
                    selectedRegion[chrom].append((chrom, Region_start, Region_end, Nwin, extremeValue))
                    # process this win
                    if i != len(selectedWinMap[chrom]):
                        mergedRegion = [selectedWinMap[chrom][i]]
                        i += 1
                i += 1
    #             except IndexError:
    #                 print(i,len(selectedWinMap[chrom]),selectedWinMap[chrom])
    #                 exit(-1)
            else:
                Region_start = int(mergedRegion[0][1]) * slideSize - upextend
                Region_end = int(mergedRegion[-1][1]) * slideSize + winWidth + downextend
                Nwin = len(mergedRegion)
                extremeValues = []
                for e in mergedRegion:
                    if options.winType == "winvalue":
                        extremeValues.append(float(e[4]))
                    elif options.winType == "zvalue": 
                        extremeValues.append(float(e[5]))
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue = max(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue = min(extremeValues)            
                selectedRegion[chrom].append((chrom, Region_start, Region_end, Nwin, extremeValue))        
        final_table = {}
        # final_table={region:[tp2,tp2,],region:}
        for chrom in selectedRegion:
            for region in selectedRegion[chrom]:
                intervalmap[minvalue, minvalue + dincrease] += winGenome.collectTrscptInWin(dbtools, TranscriptGenetable, region)
        minvalue += dincrease
    kaksMap={}
    for line in dndsfile:
        linelist=re.split(r'\s+',line)
        try:
            kaksMap[linelist[0].strip()]=float(linelist[dndscolidx])
        except:
            print(line)
    meankaksmap={}
    for a,b in sorted(intervalmap.keys()):
        meankaksmap[a,b]=0
        n=0
        for trscptrec in intervalmap[a,b]:
            if trscptrec[0].strip() in kaksMap:
                if int(kaksMap[trscptrec[0].strip()])<2:
                    n+=1
                    meankaksmap[a,b]+=kaksMap[trscptrec[0].strip()]
                               
                print(str(a)+str(b),str(kaksMap[trscptrec[0].strip()]),trscptrec[0],sep="\t",file=outfile)
        else:
            print(str(a)+str(b),str(meankaksmap[a,b]/n))
    winGenome.windbtools.drop_table(winGenome.wintabletextvalueallwin)
    winGenome.windbtools.drop_table(winGenome.wintablewithoutNA)    
    
    
    
    dndsfile.close()
    outfile.close()
