from optparse import OptionParser
import os,numpy
import re, sys, time,math

from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm


SLEEP_FOR_NEXT_TRY=3
'''
Created on 2013-9-3

@author: liurui
'''
parser = OptionParser()
parser.add_option("-i", "--winfile", dest="winfileName",
                  help="winfileName only one", metavar="FILE")

parser.add_option("-t", "--threshold", dest="threshold", help="conflict with -p")
parser.add_option("-p", "--percentage", dest="percentage",default=None, help="conflict with -t")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-x", "--morethan_lessthan", dest="morethan_lessthan", help="m:morethan or l:lessthan")
parser.add_option("-u", "--upextend", dest="upextend", help="upextend")
parser.add_option("-d", "--downextend", dest="downextend", help="downextend")
parser.add_option("-s","--slideSize",dest="slideSize",default="20000",help="win slide size")
parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
parser.add_option("-X","--winType",dest="winType",default="zvalue",help="winvalue or zvalue")
parser.add_option("-N","--mergeNA",dest="mergeNA",default=False,help="winvalue or zvalue")
parser.add_option("-n","--numberofoutlier_to_NearestGene",dest="numberofoutlier_to_NearestGene",default=0,help="number of outlier value")
parser.add_option("-a","--anchorfile",dest="anchorfile",default=None,help="winvalue or zvalue")
                                                                                                                                                          
(options, args) = parser.parse_args()

#if len(sys.argv) != 6:
#    print("python findTrscpt.py [winFile1] [tempwinDBName] [threshold] [outfilename] [m/l]")
#    exit(-1)
upextend=int(options.upextend);slideSize=int(options.slideSize);winWidth=int(options.winWidth)
downextend=int(options.downextend)
winFileName7Field = options.winfileName
if options.anchorfile!=None:
    anchorDATASTRUCTURE,reverseAnchorDATASTRUCTURE=Util.loadAnchorFile(options.anchorfile)
    ##############
    winfile=open(winFileName7Field,'r')
    title=winfile.readline()
    winMap={}#{scaffold:[(startpos,endpos,noofsnp,winvalue,zvalue),(),(),,,]}
    for line in winfile:
        linelist=re.split(r"\s+",line.strip())
        if linelist[0].strip()  in winMap:
            winMap[linelist[0].strip()].append((int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6]))
        else:
            winMap[linelist[0].strip()]=[(int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6])]
    winfile.close()
    ##################winfile has been loaded into memonery
    ##################reZ-transform the winvalue by seperate the autochromosome and sex chromosome
    
    winCrossGenomeMap={"autosome":[],"Z":[],"W":[],"X":[],"Y":[]}
    winFileName7Field=winFileName7Field+"sexchromseperatestandard"
    
    
    f=open(winFileName7Field,'w')
    print(title,end="",file=f)
    for scaffold in winMap.keys():
        for startpos,endpos,noofsnp,winvalue,zvalue in winMap[scaffold]:
            if  re.search(r"^[1234567890\.e-]+$",winvalue)==None:
                continue
            if scaffold not in reverseAnchorDATASTRUCTURE:
                winCrossGenomeMap["autosome"].append(float(winvalue))
            elif "Z" in reverseAnchorDATASTRUCTURE[scaffold] or "z" in reverseAnchorDATASTRUCTURE[scaffold]:
                winCrossGenomeMap["Z"].append(float(winvalue))
            elif "W" in reverseAnchorDATASTRUCTURE[scaffold] or "w" in reverseAnchorDATASTRUCTURE[scaffold]:
                winCrossGenomeMap["W"].append(float(winvalue))
            elif "X" in reverseAnchorDATASTRUCTURE[scaffold] or "x" in reverseAnchorDATASTRUCTURE[scaffold]:
                winCrossGenomeMap["X"].append(float(winvalue))
            elif "Y" in reverseAnchorDATASTRUCTURE[scaffold] or "y" in reverseAnchorDATASTRUCTURE[scaffold]:
                winCrossGenomeMap["Y"].append(float(winvalue)) 
            else:
                winCrossGenomeMap["autosome"].append(float(winvalue))
    autoexception=numpy.mean(winCrossGenomeMap["autosome"])
    autostd1=numpy.std(winCrossGenomeMap["autosome"],ddof=1)
    sexexception=numpy.mean(winCrossGenomeMap["Z"]+winCrossGenomeMap["W"]+winCrossGenomeMap["X"]+winCrossGenomeMap["Y"])
    sexstd1=numpy.std(winCrossGenomeMap["Z"]+winCrossGenomeMap["W"]+winCrossGenomeMap["X"]+winCrossGenomeMap["Y"],ddof=1)
    for scaffold in sorted(winMap.keys()):
        winNo=0
        for startpos,endpos,noofsnp,winvalue,zvalue in winMap[scaffold]:
            if re.search(r"^[1234567890\.e-]+$",winvalue)!=None:
                if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                    zscore=(float(winvalue)-autoexception)/autostd1
                else:
                    zscore=(float(winvalue)-sexexception)/sexstd1
                print(scaffold,winNo,startpos,endpos,noofsnp,winvalue,zscore,sep="\t",file=f)
            else:
                print(scaffold,winNo,startpos,endpos,noofsnp,winvalue,zvalue,sep="\t",file=f)
            winNo+=1
    f.close()
re.search(r"[^/]*$",winFileName7Field).group(0)
if re.search(r'^.*/',options.outfileprename)!=None:
    path=re.search(r'^.*/',options.outfileprename).group(0)
else:
    a = os.popen("pwd")
    path=a.readline().strip()+"/"
    a.close()

total_outliers=int(options.numberofoutlier_to_NearestGene)
threshold = options.threshold
percentage = options.percentage
outfilename=options.outfileprename
morethan_lessthan=options.morethan_lessthan.lower()
mergeNA=options.mergeNA
print(mergeNA)
if percentage!=None and threshold!=None:
    print("-t conflict with -p")
    exit(-1)
#gene_sample_venn="gene_sample_venn"ninglabvariantdata_tmp
vcftable=None
outfile=open(outfilename,'w')
print("chrNo\tRegion_start\tRegion_end\tNoofWin\textram"+options.winType+"\ttranscpt\tgeneID",file=outfile)
outfileNameWINwithGENE=path+re.search(r"[^/]*$",winFileName7Field).group(0)+".wincopywithgene"

    
if __name__ == '__main__':
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname) 
    winGenome = Util.WinInGenome(Util.ghostdbname, winFileName7Field)
    time.sleep(SLEEP_FOR_NEXT_TRY)
    winGenome.appendGeneName(Util.TranscriptGenetable, genomedbtools, winWidth, slideSize, outfileNameWINwithGENE,upextend,downextend,(total_outliers,morethan_lessthan))
    selectWinNos="threshold method"
    if percentage!=None:
        totalWin = winGenome.windbtools.operateDB("select", "select count(*) from " + winGenome.wintablewithoutNA)[0][0]
        selectWinNos = int(float(percentage) * totalWin)
        if morethan_lessthan == "m" or morethan_lessthan == "M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+options.winType+" desc limit 0," + str(selectWinNos))
            print("select * from "+winGenome.wintablewithoutNA + " where 1 order by zvalue desc limit 0," + str(selectWinNos))
        elif morethan_lessthan == "l" or morethan_lessthan == "L":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 order by "+options.winType+" asc limit 0," + str(selectWinNos))
            print("select * from " + winGenome.wintablewithoutNA + " where 1 order by "+options.winType+" asc limit 0," + str(selectWinNos))
    elif threshold!=None:
        if morethan_lessthan=="m" or morethan_lessthan=="M":
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 and "+options.winType+">=" + threshold)
        elif morethan_lessthan=="l" or morethan_lessthan=="L":
            print("select", "select * from " + winGenome.wintablewithoutNA + " where "+options.winType+"!= 'NA' and "+options.winType+"<=" + threshold)
            selectedWins = winGenome.windbtools.operateDB("select", "select * from " + winGenome.wintablewithoutNA + " where 1 and "+options.winType+"<=" + threshold)
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

    selectedRegion={}

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
                    extremeValue=min(extremeValues)
                elif morethan_lessthan == "l" or morethan_lessthan == "L":
                    extremeValue=max(extremeValues)
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
                extremeValue=min(extremeValues)
            elif morethan_lessthan == "l" or morethan_lessthan == "L":
                extremeValue=max(extremeValues)            
            selectedRegion[chrom].append((chrom,Region_start,Region_end,Nwin,extremeValue))
    if mergeNA!=False and int(mergeNA)>0:
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
                if len(wincount_to_determine)==0 and len(wincount_to_add)<= int(mergeNA):
                    if morethan_lessthan == "m" or morethan_lessthan == "M":
                        extremeValue=min(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    elif morethan_lessthan == "l" or morethan_lessthan == "L":
                        extremeValue=max(selectedRegion[chrom][i][4],selectedRegion[chrom][i-1][4])
                    selectedRegion[chrom][i]=(chrom,selectedRegion[chrom][i-1][1],selectedRegion[chrom][i][2],selectedRegion[chrom][i-1][3]+selectedRegion[chrom][i][3]+len(wincount_to_add),extremeValue)
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
            if total_outliers>0:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,Util.TranscriptGenetable,region,upextend,downextend,True)
            else:
                final_table[region]=winGenome.collectTrscptInWin(genomedbtools,Util.TranscriptGenetable,region,upextend,downextend)
#process top outlier values

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

    winGenome.windbtools.drop_table(winGenome.wintabletextvalueallwin)
    winGenome.windbtools.drop_table(winGenome.wintablewithoutNA)
    outfile.close()
