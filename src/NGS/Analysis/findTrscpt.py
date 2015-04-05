from optparse import OptionParser
import os
import re, sys, time

from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm


SLEEP_FOR_NEXT_TRY=3
'''
Created on 2013-9-3

@author: liurui
'''
parser = OptionParser()
parser.add_option("-i", "--winfile", dest="winfileName",
                  help="winfileName", metavar="FILE")
parser.add_option("-1", "--tempDBname", dest="tempdbname", help="dbname")
parser.add_option("-t", "--threshold", dest="threshold", help="conflict with -p")
parser.add_option("-p", "--percentage", dest="percentage",default=None, help="conflict with -t")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-x", "--morethan_lessthan", dest="morethan_lessthan", help="m:morethan or l:lessthan")
parser.add_option("-T", "--trscptable", dest="trscptable", help="trscptable")
parser.add_option("-u", "--upextend", dest="upextend", help="upextend")
parser.add_option("-d", "--downextend", dest="downextend", help="downextend")
parser.add_option("-s","--slideSize",dest="slideSize",default="20000",help="win slide size")
parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
parser.add_option("-X","--winType",dest="winType",default="zvalue",help="winvalue or zvalue")
parser.add_option("-N","--mergeNAorNOT",dest="mergeNAorNOT",action="store_true",default=False,help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
#if len(sys.argv) != 6:
#    print("python findTrscpt.py [winFile1] [tempwinDBName] [threshold] [outfilename] [m/l]")
#    exit(-1)
upextend=int(options.upextend);slideSize=int(options.slideSize);winWidth=int(options.winWidth)
downextend=int(options.downextend)
winFileName7Field = options.winfileName
if re.search(r'^.*/',options.winfileName)!=None:
    path=re.search(r'^.*/',options.winfileName).group(0)
else:
    a = os.popen("pwd")
    path=a.readline().strip()+"/"
    a.close()

tempwinDBName = options.tempdbname
threshold = options.threshold
percentage = options.percentage
outfilename=path+options.outfileprename
morethan_lessthan=options.morethan_lessthan
TranscriptGenetable=options.trscptable.strip()
mergeNAorNOT=options.mergeNAorNOT
print(mergeNAorNOT)
if percentage!=None and threshold!=None:
    print("-t conflict with -p")
    exit(-1)
#gene_sample_venn="gene_sample_venn"
vcftable=None
outfile=open(outfilename,'w')
print("chrNo\tRegion_start\tRegion_end\tNoofWin\textram"+options.winType+"\ttranscpt\tgeneID",file=outfile)
outfileNameWINwithGENE=winFileName7Field+".wincopywithgene"
if __name__ == '__main__':
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
#    dbtools.operateDB("alter","alter table "+gene_sample_venn+" add "+outfilename+" smallint(3) default 0") 
    winGenome = Util.WinInGenome(tempwinDBName, winFileName7Field)
    time.sleep(SLEEP_FOR_NEXT_TRY)
    winGenome.appendGeneName(TranscriptGenetable, genomedbtools, winWidth, slideSize, outfileNameWINwithGENE)
    selectWinNos="threshold method"
    if percentage!=None:
        totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintablewithoutNA)[0][0]
        selectWinNos = int(float(percentage) * totalWin)
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+" != 'NA' order by "+options.winType+" desc limit 0," + str(selectWinNos))
            print("select * from "+winGenome.wintablewithoutNA + " where zvalue != 'NA' order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+" != 'NA' order by "+options.winType+" asc limit 0," + str(selectWinNos))
            print("select * from " + winGenome.wintablewithoutNA + " where "+options.winType+" != 'NA' order by "+options.winType+" asc limit 0," + str(selectWinNos))
    elif threshold!=None:
        if morethan_lessthan=="m" or morethan_lessthan=="M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+"!= 'NA' and "+options.winType+">=" + threshold)
        elif morethan_lessthan=="l" or morethan_lessthan=="L":
            print("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+"!= 'NA' and "+options.winType+"<=" + threshold)
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+"!= 'NA' and "+options.winType+"<=" + threshold)
        selectWinNos=len(selectedWins)
    selectedWins.sort(key=lambda listRec:float(listRec[5]))
    if selectWinNos==0:
        outfile.close()
        print("selectWinNos==0")
        exit(0)
    print(outfilename,selectWinNos,"~=",len(selectedWins),selectedWins[0],selectedWins[-1])
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
        i=1
        while i < len(selectedWinMap[chrom]):
#             print(chrom,selectedWinMap[chrom][i])
#             try:
            if int(selectedWinMap[chrom][i-1][1])+1==int(selectedWinMap[chrom][i][1]):#continues win
                mergedRegion.append(selectedWinMap[chrom][i])
            else:#not continues
                #process last region
                Region_start=int(mergedRegion[0][1])*slideSize
                Region_end=int(mergedRegion[-1][1])*slideSize+winWidth
                Nwin=len(mergedRegion)
                extremeValues=[]
                for e in mergedRegion:
                    if options.winType=="winvalue":
                        extremeValues.append(float(e[5]))
                    elif options.winType=="zvalue": 
                        extremeValues.append(float(e[6]))
                if morethan_lessthan == "m" or morethan_lessthan == "M":
                    extremeValue=max(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=min(extremeValues)
                selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
                #process this win
                mergedRegion=[selectedWinMap[chrom][i]]
            i+=1
#             except IndexError:
#                 print(i,len(selectedWinMap[chrom]),selectedWinMap[chrom])
#                 exit(-1)
        else:
            Region_start=int(mergedRegion[0][1])*slideSize
            Region_end=int(mergedRegion[-1][1])*slideSize+winWidth
            Nwin=len(mergedRegion)
            extremeValues=[]
            for e in mergedRegion:
                if options.winType=="winvalue":
                    extremeValues.append(float(e[5]))
                elif options.winType=="zvalue": 
                    extremeValues.append(float(e[6]))
            if morethan_lessthan == "m" or morethan_lessthan == "M":
                extremeValue=max(extremeValues)
            elif morethan_lessthan == "l" or morethan_lessthan == "L":
                extremeValue=min(extremeValues)            
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
    if mergeNAorNOT:
        for chrom in selectedRegion:
            selectedRegion[chrom].sort(key=lambda listRec: int(listRec[1]))
            i=1
            idxlist_to_pop=[]
            while i <len(selectedRegion[chrom]):
                winNo_end=str(int(selectedRegion[chrom][i][1]/slideSize))
                winNo_start=str(int((selectedRegion[chrom][i-1][2]-winWidth)/slideSize))
                print("select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and  winNo<"+winNo_end)
                wincount_to_determine=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintablewithoutNA + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
                wincount_to_add=winGenome.windbtools.operateDB("select","select * from "+ winGenome.wintabletextvalueallwin + " where "+" chrID='"+chrom+"' and winNo>"+winNo_start+" and winNo<"+winNo_end)
                if len(wincount_to_determine)==0:
                    if morethan_lessthan == "m" or morethan_lessthan == "M":
                        extremeValue=max(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    elif morethan_lessthan == "l" or morethan_lessthan == "L":
                        extremeValue=min(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    selectedRegion[chrom][i]=(chrom,selectedRegion[chrom][i-1][1],selectedRegion[chrom][i][2],selectedRegion[chrom][i][3]+len(wincount_to_add),extremeValue)
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
            final_table[region]=winGenome.collectTrscptInWin(genomedbtools,TranscriptGenetable,region,upextend,downextend)
#    for win in selectedWins:
#        winRegion=(win,upextend,downextend)
#        winGenome.collectTrscptInWin(dbtools, TranscriptGenetable, vcftable, winRegion)
    for chrom in winGenome.chromOrder:
        if chrom not in selectedRegion:
            continue
        for region in selectedRegion[chrom]:
            if chrom.strip()==region[0].strip():
                tcpts=""
                gnames=""
                for tcpt in final_table[region]:
                    tcpts+=(tcpt[0]+",")
                    if tcpt[2].strip()!="":
                        gnames+=(tcpt[2]+",")
                print("\t".join(map(str,region)),tcpts[:-1],gnames[:-1],sep="\t",file=outfile)                  
#     for region in sorted(final_table.keys()):
#         tcpts=""
#         for tcpt in final_table[region]:
#             tcpts+=(tcpt[0]+"\t")
#         print("\t".join(map(str,region)),tcpts,sep="\t",file=outfile)
       
    winGenome.windbtools.drop_table(winGenome.wintabletextvalueallwin)
    winGenome.windbtools.drop_table(winGenome.wintablewithoutNA)
#        for gene in winGenome.winContainTrscptMap[win]:
#            print(gene)
#            print("update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
#            dbtools.operateDB("update","update "+gene_sample_venn+" set "+outfilename+"=1 where geneID='"+gene[0]+"'")
    outfile.close()
#    winGenome.windbtools.drop_table(winGenome.wintable)
#    winGenome.windbtools.disconnect()
