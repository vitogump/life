'''
Created on 2016-12-13

@author: liurui
'''
import math
from optparse import OptionParser
import re,copy


parser = OptionParser()
parser.add_option("-m", "--mapfile", dest="mapfile", help="optional")
parser.add_option("-a", "--anchorfile", dest="anchorfile", help="")
parser.add_option("-i", "--winfile", dest="winfile", help="")
parser.add_option("-w", "--winwidth", dest="winwidth", help="")
parser.add_option("-s", "--slidesize", dest="slidesize", help="")
(options, args) = parser.parse_args()
newanchorfilehandler=options.anchorfile
if options.mapfile:
    scaffoldmap={}
    mapfile=open(options.mapfile,'r')
    for line in mapfile:
        linelist=re.split(r"\s+",line.strip())
        scaffoldmap[linelist[0].strip().lower()]=linelist[1].strip()
    oldanchorfilehandler=open(options.anchorfile,'r')
    newanchorfilehandler=open(options.anchorfile+"changed",'w')
    for line in oldanchorfilehandler:
        linelist=re.split(r"\s+",line.strip())
        if linelist[3] in scaffoldmap:
            linelist[3]=scaffoldmap[linelist[3].strip().lower()]
        print(*linelist,sep="\t",file=newanchorfilehandler)
    oldanchorfilehandler.close()
    newanchorfilehandler.close()
        
    
anchorDATASTRUCTURE={}
"""
{chr1:[(53353,53806,scaffold451,558997,558537,-),(57200,62371,scaffold451,553669,548504,-),(),,],chr2:[],,,,}
"""
reverseAnchorDATASTRUCTURE={}
"""
{scaffold451:{chr1:[0,1,2,,,,]},C17734302:{chr1:[idx]}}  idx is idx in the list of anchorDATASTRUCTURE[chr1] 
"""
if options.mapfile:
    
    newanchorfilehandler=open(options.anchorfile+"changed",'r')
else:
    newanchorfilehandler=open(options.anchorfile,'r')
for line in newanchorfilehandler:
    linelist=re.split(r"\s+",line.strip())
    if linelist[0].strip() in anchorDATASTRUCTURE:
        anchorDATASTRUCTURE[linelist[0].strip()].append((int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip()))
    else:
        anchorDATASTRUCTURE[linelist[0].strip()]=[(int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip())]
    #fill reverseAnchorDATASTRUCTURE
    if linelist[3].strip() in reverseAnchorDATASTRUCTURE:
        if linelist[0].strip() in reverseAnchorDATASTRUCTURE[linelist[3].strip()]:
            reverseAnchorDATASTRUCTURE[linelist[3].strip()][linelist[0].strip()].append(len(anchorDATASTRUCTURE[linelist[0].strip()])-1)
        else:
            reverseAnchorDATASTRUCTURE[linelist[3].strip()]={linelist[0].strip():[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]}
    else:
        reverseAnchorDATASTRUCTURE[linelist[3].strip()]={linelist[0].strip():[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]}
newanchorfilehandler.close()
if __name__ == '__main__':
    winfile=open(options.winfile,'r')
    title=winfile.readline()
    winMap={}#{scaffold:[(startpos,endpos,noofsnp,winvalue,zvalue),(),(),,,]}
    for line in winfile:
        linelist=re.split(r"\s+",line.strip())
        if linelist[0].strip()  in winMap:
            winMap[linelist[0].strip()].append((int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6]))
        else:
            winMap[linelist[0].strip()]=[(int(linelist[2]),int(linelist[3]),int(linelist[4]),linelist[5],linelist[6])]
    winfile.close()
    #winfile has been loaded into memonery
    #print a new winfile mark auto or sex chromosome
    f=open(options.winfile)
    winMapMarked=copy.deepcopy(winMap)

    outwinfile=open(options.winfile+"arrangemented",'w')
    print(title.strip(),file=outwinfile)
    for chrom in sorted(anchorDATASTRUCTURE.keys()):
        idx=0#it seems not useful
        if anchorDATASTRUCTURE[chrom][idx][5]=="-":
            regionstart=anchorDATASTRUCTURE[chrom][idx][4]
            regionend=anchorDATASTRUCTURE[chrom][idx][3]
        elif anchorDATASTRUCTURE[chrom][idx][5]=="+":
            regionstart=anchorDATASTRUCTURE[chrom][idx][3]
            regionend=anchorDATASTRUCTURE[chrom][idx][4]
        lastscaffold=anchorDATASTRUCTURE[chrom][idx][2]
        if lastscaffold not in winMap:
            for startpos,endpos,scaffold,sstartpos,sendpos,foward_reverse in anchorDATASTRUCTURE[chrom]:
                idx+=1
                if scaffold in winMap:
                    lastscaffold=scaffold
                    break

        for startpos,endpos,scaffold,sstartpos,sendpos,foward_reverse in anchorDATASTRUCTURE[chrom][idx:]:
            if lastscaffold==scaffold:
#                 print("this code block counting the continue stretch . and Should consider wither there are some gap in it")
                if foward_reverse=="-":
                    regionstart=min(regionstart,sendpos)
                    regionend=max(regionend,sstartpos)
                else:
                    regionstart=min(regionstart,sstartpos)
                    regionend=max(regionend,sendpos)
                
            else:
                print("after determining the start or the end of the stretch, arranging the scaffolds which in the stretch")
                if regionstart < int(options.winwidth):
                    winstartNo=0
                else:
                    winstartNo=math.ceil((regionstart-int(options.winwidth))/int(options.slidesize))
                if regionend<int(options.winwidth):
                    winendNo=0
                else:
                    winendNo=math.ceil((regionend-int(options.winwidth))/int(options.slidesize))
                if anchorDATASTRUCTURE[chrom][idx-1][5]=="-":#last scaffold so here is idx-1
                    for i in range(winstartNo,winendNo+1)[::-1]:
                        if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                            winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["auto"])
                        else:
                            winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                        print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
                else:
                    for i in range(winstartNo,winendNo+1):
                        if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                            winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["auto"])
                        else:
                            winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                        print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
                print("new region")
                if foward_reverse=="-":
                    regionstart=sendpos
                    regionend=sstartpos
                elif foward_reverse=="+":
                    regionstart=sstartpos
                    regionend=sendpos
                if scaffold not in winMap:
                    continue
                lastscaffold=scaffold
            idx+=1
        else:
            if regionstart < int(options.winwidth):
                winstartNo=0
            else:
                winstartNo=math.ceil((regionstart-int(options.winwidth))/int(options.slidesize))
            if regionend<int(options.winwidth):
                winendNo=0
            else:
                winendNo=math.ceil((regionend-int(options.winwidth))/int(options.slidesize))
            if anchorDATASTRUCTURE[chrom][idx-1][5]=="-":
                for i in range(winstartNo,winendNo+1)[::-1]:
                    if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["auto"])
                    else:
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                    print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
            else:
                for i in range(winstartNo,winendNo+1):
                    if scaffold not in reverseAnchorDATASTRUCTURE or ("Z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "z" not in reverseAnchorDATASTRUCTURE[scaffold] and  "W" not in reverseAnchorDATASTRUCTURE[scaffold] and "w" not in reverseAnchorDATASTRUCTURE[scaffold] and "X" not in reverseAnchorDATASTRUCTURE[scaffold] and "x" not in reverseAnchorDATASTRUCTURE[scaffold] and "Y" not in reverseAnchorDATASTRUCTURE[scaffold] and "y" not in reverseAnchorDATASTRUCTURE[scaffold]):
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["auto"])
                    else:
                        winMapMarked[lastscaffold][i]=tuple((list(winMapMarked[lastscaffold][i]))+["Z"])
                    print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],lastscaffold,sep="\t",file=outwinfile)
            
#             for i in range(winstartNo,winendNo+1):
#                 print(chrom,i,winMap[lastscaffold][i][0],winMap[lastscaffold][i][1],winMap[lastscaffold][i][2],winMap[lastscaffold][i][3],winMap[lastscaffold][i][4],scaffold,sep="\t",file=outwinfile)
    outwinfile.close()
    outwinfile=open(options.winfile+"marked",'w')
    print(title.strip()+"\tmark",file=outwinfile)
    for scaffold in sorted(winMapMarked.keys()):
        for winNo in range(len(winMapMarked[scaffold])):
            if len(winMapMarked[scaffold][winNo])==5 or winMapMarked[scaffold][winNo][5]=="auto":
                print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"autochromosome",sep="\t",file=outwinfile)
            else:
                print(scaffold,winNo,winMapMarked[scaffold][winNo][0],winMapMarked[scaffold][winNo][1],winMapMarked[scaffold][winNo][2],winMapMarked[scaffold][winNo][3],winMapMarked[scaffold][winNo][4],"sexchromosome",sep="\t",file=outwinfile)
    outwinfile.close()