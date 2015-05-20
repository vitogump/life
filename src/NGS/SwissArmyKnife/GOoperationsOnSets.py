'''
Created on 2014-5-10

@author: liurui
'''
from optparse import OptionParser
import re

#if len(sys.argv) < 4:
#    print("python GOoperationsOnSets.py [selectedgenefileName] [selectedgenefileName] [selectedgenefileName]....-g [gotablefile] -o [outfileprename] ")
#    exit(-1)
parser = OptionParser()

parser.add_option("-g", "--gotablefile", dest="gotablefile", help="gtffile")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
gotablefile=open(options.gotablefile,'r')
if len(args)==1:
    intersectionSetOutfile=open(options.outfileprename,'w')
else:
    intersectionSetOutfile=open(options.outfileprename+".intersection_set",'w')
    unionSetOutfile=open(options.outfileprename+".union_set",'w')

if __name__ == '__main__':
    print(args[:])
    tablegeneSet={}
    """
    tablegeneSet={selectedgenefileName:{tpID1,tpID2,,,,,,}}
    """
    filetrsptmap={}
    """
    filetrsptmap={selectedgenefileName:{tpID1:[(selectedRegion,extremeValue)]},selectedgenefileName:{tpID2:[(selectedRegion,extremeValue)]},,,,,,}
    """
    for selectedgenefileName in args[:]:
        selectedfile=open(selectedgenefileName,'r')
        selectedgenefileNamewithoutpath = re.search(r'[^/]*$',selectedgenefileName).group(0)
        trsptMap={}
        for line in selectedfile:
            linelist=re.split(r"\s+",line.strip())
            if len(linelist)<6:
                continue
            else:
                chrom=linelist[0]
                startPos=linelist[1]
                endPos=linelist[2]
                extremeValue=linelist[4]
                tableRegion=chrom+":"+startPos+"-"+endPos
                tpIDlist=re.search(r",",linelist[5].strip())
                for tpID in tpIDlist:
                    if tpID not in trsptMap:
                        trsptMap[tpID]=[(tableRegion,extremeValue)]
                    else:
                        if (float(extremeValue)<0 and float(extremeValue)<float(trsptMap[tpID][0][1])) or (float(extremeValue)>0 and float(extremeValue)>float(trsptMap[tpID][0][1])):
                            trsptMap[tpID].insert(0,(tableRegion,extremeValue))
                        else:
                            trsptMap[tpID].append((tableRegion,extremeValue))
        filetrsptmap[re.search(r"^.*\.",selectedgenefileNamewithoutpath).group(0)[:-1]]=trsptMap
        tablegeneSet[re.search(r"^.*\.",selectedgenefileNamewithoutpath).group(0)[:-1]]=set(trsptMap.keys())
        selectedfile.close()
        
        
    intersectionSet = tablegeneSet[re.search(r"^.*\.",selectedgenefileNamewithoutpath).group(0)[:-1]]
    unionSet = tablegeneSet[re.search(r"^.*\.",selectedgenefileNamewithoutpath).group(0)[:-1]]

    for selectedgenetable in tablegeneSet:
         intersectionSet=intersectionSet & tablegeneSet[selectedgenetable]
         unionSet=unionSet | tablegeneSet[selectedgenetable]
        
#================================================ operations On Sets finished,begin GO annotion
    title = gotablefile.readline()
    titlelist= [e.strip().lower() for e in re.split(r",",title)]
    geneididx=titlelist.index("ensembl gene id")
    tpididx=titlelist.index("ensembl transcript id")
    gotermaccessionidx=titlelist.index("go term accession")
    gotermNameidx=titlelist.index("go term name")
    godomainidx=titlelist.index("go domain")
    genenameidx=titlelist.index("associated gene name")
    gotable={}
    """
    {tp_id1:(geneID,geneName,bp,cc,mf),tp_id2:(geneID,geneName,bp,cc,mf),,,,,,}
    """
    bp="";cc="";mf="";geneName="";geneID=""
    for termline in  gotablefile:
        termlist=re.split(r",",termline)
        if termlist[tpididx].strip() in gotable:
            if termlist[godomainidx].lower().strip()=="biological_process":
                bp+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
            elif termlist[godomainidx].lower().strip()=="cellular_component":
                cc+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"                 
            elif termlist[godomainidx].lower().strip()=="molecular_function":
                mf+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
           
        else:#new gene start
            gotable[termlist[tpididx].strip()]=(geneID,geneName,bp,cc,mf)#both ok to first time or second time
            geneName=termlist[genenameidx]
            geneID=termlist[geneididx]
            bp="";cc="";mf=""
            if termlist[godomainidx].lower().strip()=="biological_process":
                bp+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"
            elif termlist[godomainidx].lower().strip()=="cellular_component":
                cc+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"            
            elif termlist[godomainidx].lower().strip()=="molecular_function":
                mf+=termlist[gotermaccessionidx]+";"+termlist[gotermNameidx]+";"

# write output file 
    title="ensembl transcript id\t"
    for selectedgenetablefile in sorted(filetrsptmap.keys()):
        title+=selectedgenetablefile+"Region\textremeValue\t"
    title+="geneID\tgeneName\tbiological_process domain\tcellular_component domain\tmolecular_function domain"
    print(title,file=intersectionSetOutfile)
            
    for intersectionTrspt in intersectionSet:
        if intersectionTrspt not in gotable:
            print(intersectionTrspt,"don't have go annotion")
            continue
        print(intersectionTrspt,end="\t",file=intersectionSetOutfile)
        for selectedgenetablefile in sorted(filetrsptmap.keys()):
            print(filetrsptmap[selectedgenetablefile][intersectionTrspt][0][0],filetrsptmap[selectedgenetablefile][intersectionTrspt][0][1],sep="\t",end="\t",file=intersectionSetOutfile)
        print(gotable[intersectionTrspt][0],gotable[intersectionTrspt][1],gotable[intersectionTrspt][2],gotable[intersectionTrspt][3],gotable[intersectionTrspt][4],sep="\t",file=intersectionSetOutfile)
    if len(args)>1:
        print(title,file=unionSetOutfile)
        for unionTrspt in unionSet:
            if unionTrspt not in gotable:
                print(unionTrspt,"don't have go annotion")
                continue
            print(unionTrspt,end="\t",file=unionSetOutfile)
            for selectedgenetablefile in sorted(filetrsptmap.keys()):
                if unionTrspt in filetrsptmap[selectedgenetablefile]:
                    print(filetrsptmap[selectedgenetablefile][unionTrspt][0][0],filetrsptmap[selectedgenetablefile][unionTrspt][0][1],sep="\t",end="\t",file=unionSetOutfile)
                else:
                    print("NA\tNA\t",end="",file=unionSetOutfile)
            print(gotable[unionTrspt][0],gotable[unionTrspt][1],gotable[unionTrspt][2],gotable[unionTrspt][3],gotable[unionTrspt][4],sep="\t",file=unionSetOutfile)    
        unionSetOutfile.close()
    intersectionSetOutfile.close()