
#import Make_Picture
'''
Created on 2013-8-11

@author: rui
'''

from optparse import OptionParser
import os,numpy
import re

from NGS.RUtil import *
from src.NGS.BasicUtil import geneUtil


parser = OptionParser()

parser.add_option("-o","--pathoutputfilename",dest="pathoutputfilename",help="default infile2_infile1")#
# parser.add_option("-P","--positive_withgenename",dest="multiple_positive_winfiles_withgenename",action="append",nargs=2,default=[],help="on top,file name and threshold")#
# parser.add_option("-N","--negtive_withgenename",dest="multiple_negtive_winfiles_withgenename",action="append",nargs=2,default=[],help="at bottom,file name and threshold")#

parser.add_option("-p","--positive",dest="multiple_positive_winfiles",action="append",nargs=3,default=[],help="on top,filename ond threshold and outpre")#
parser.add_option("-n","--negtive",dest="multiple_negtive_winfiles",action="append",nargs=3,default=[],help="at bottom")#
parser.add_option("-a","--allvalue",dest="multiple_allvalue_winfiles",action="append",nargs=3,default=[],help="at bottom")#

parser.add_option("-A","--anchorfile",dest="anchorfile",default=None,help="winvalue or zvalue")
parser.add_option("-g", "--gotablefile", dest="gotablefile", help="gotable title with :Ensembl Gene ID    Ensembl Transcript ID    GO Term Accession    GO Term Evidence Code    GO domain    GO Term Name    GO Term Definition,order and upper/lower case is arbitrarily")

parser.add_option("-x", "--threshold_percentage", dest="threshold_percentage",help="t / p", metavar="FILE")
parser.add_option("-e", "--distalextend", dest="distalextend",default="180000",help="t / p", metavar="FILE")
parser.add_option("-f", "--trscptfound", dest="trscptfound",action="store_true",default=False, help="outfileprename")
parser.add_option("-S", "--splitintopart", dest="splitintopart",default=1, help="split winfile into part")
parser.add_option("-u", "--upextend", dest="upextend", help="upextend")
parser.add_option("-d", "--downextend", dest="downextend", help="downextend")
parser.add_option("-s","--slideSize",dest="slideSize",default="20000",help="win slide size")
parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
parser.add_option("-X","--winType",dest="winType",default="zvalue",help="winvalue or zvalue")
parser.add_option("-N","--mergeNA",dest="mergeNA",default=False,help="winvalue or zvalue")
parser.add_option("-m", "--minmaxSNP", dest="minmaxSNP",default="0", help="upextend")
parser.add_option("-r", "--removegenelistfile", dest="removegenelistfile", help="upextend")
# parser.add_option("-t","--numberofoutlier_to_NearestGene",dest="numberofoutlier_to_NearestGene",default=0,help="number of outlier value,for example top 10")
(options, args) = parser.parse_args()
columnname=options.winType.strip()

def standardseparately(anchorfile,winfilein,upextend,downextend,winWidth,slideSize):
    anchorDATASTRUCTURE={}
    """
    {chr1:[(53353,53806,scaffold451,558997,558537,-),(57200,62371,scaffold451,553669,548504,-),(),,],chr2:[],,,,}
    """
    reverseAnchorDATASTRUCTURE={}
    """
    {scaffold451:{chr1:[0,1,2,,,,]},C17734302:{chr1:[idx]}}  idx is idx in the list of anchorDATASTRUCTURE[chr1] 
    """
    newanchorfilehandler=open(anchorfile,'r')
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
    ##############
    winfile=open(winfilein,'r')
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
    winFileName7Field=winfilein+"sexchromseperatestandard"
    
    
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
    return winFileName7Field
    
if __name__ == '__main__':

    makeMhtGraph = Make_Picture.MakeMhtGraph()
    outfileNameWINwithGENE_Plist=[];outfileNameWIN_Plist=[]
    outfileNameWINwithGENE_Nlist=[];outfileNameWIN_Nlist=[]
    uniontpidlist=[];intertpidset=set()
    removed=[]
    if options.splitintopart==1:
        if options.multiple_positive_winfiles!=[]:
            for p_inputfileName,threshold_title,outbedfilename in options.multiple_positive_winfiles[:]:
                if options.anchorfile!=None:
                    p_inputfileName=standardseparately(options.anchorfile,p_inputfileName,int(options.upextend),int(options.downextend),int(options.winWidth),int(options.slideSize))
                outfileNameWIN_Plist.append(p_inputfileName)
                threshold_title_list=re.split(r"_",threshold_title.strip())
                outfileNameWINwithGENE_Plist.append((geneUtil.findTrscpt(p_inputfileName, outbedfilename, int(options.upextend), int(options.downextend), int(options.winWidth), int(options.slideSize), options.winType, "m", threshold_title_list[0], None, options.mergeNA, int(options.distalextend),options.trscptfound),threshold_title,outbedfilename))
                makeMhtGraph.makeHistonPicture(p_inputfileName, "Fst")#,"c(0,2000)","c(0,45)"
                makeMhtGraph.makeHistonPicture(outfileNameWINwithGENE_Plist[-1][0], "Fst")
                print("awk 'NR!=1{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk'$0~/^ENS/{print $0}' >"+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+options.minmaxSNP +"{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk '$0~/^ENS/{print $0}' >"+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+options.minmaxSNP +"{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk '$0!~/^ENS/{print $0}' >"+outbedfilename+".miRNAlist")
                if options.removegenelistfile!=None:
                    f=open(outbedfilename+".ENStrscptIDlist",'r')
                    mylist=[];
                    for line in f:

    #                     uniontpidlist.append(line.strip())
                        mylist.append(line.strip())
                    f.close()
                    ff=open(outbedfilename+".ENStrscptIDlist",'w')
                    removelist=[]    
############################
                    f=open(options.removegenelistfile,'r')
                    
                    for line in f:
                        removelist.append(line.strip())
                    f.close()
#########################
                    finallist=list(set(mylist)-set(removelist))
                    removed+=list(set(mylist).intersection(set(removelist)))
                    for e in finallist:
                        print(e,file=ff)
                    ff.close()                
                f=open(outbedfilename+".ENStrscptIDlist",'r')
                curset=set()
                mylist=[]
                for line in f:
                    curset.add(line.strip())
                    uniontpidlist.append(line.strip())
                    mylist.append(line.strip())
                f.close()
                if intertpidset:
                    intertpidset=intertpidset.intersection(curset)
                else:
                    intertpidset=curset

                geneUtil.GOenrichment(options.gotablefile,outbedfilename,None,list(set(mylist)),None)

                print("grep -wFf "+outbedfilename+".ENStrscptIDlist"+""" /home/bioinfo/databases/ensembleIDconvert.txt|awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d'>"""+outbedfilename+""".Homologs_human""")
                os.system("grep -wFf "+outbedfilename+".ENStrscptIDlist"+""" /home/bioinfo/databases/ensembleIDconvert.txt|awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d'>"""+outbedfilename+""".Homologs_human""")
                print("grep -wFf "+outbedfilename+".Homologs_human /home/bioinfo/databases/humangenesymbl.txt|awk '{print $3}'|sort|uniq|sed '/^$/d'>"+outbedfilename+"Homologs_human_genesymbl")
                os.system("grep -wFf "+outbedfilename+".Homologs_human /home/bioinfo/databases/humangenesymbl.txt|awk '{print $3}'|sort|uniq|sed '/^$/d'>"+outbedfilename+"Homologs_human_genesymbl")
                print("""grep -wFf """+outbedfilename+""".Homologs_human /home/bioinfo/databases/humanGO.table |awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d' > """+outbedfilename+".Homologs_humanEntrezGeneID")
                os.system("""grep -wFf """+outbedfilename+""".Homologs_human /home/bioinfo/databases/humanGO.table |awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d' > """+outbedfilename+".Homologs_humanEntrezGeneID")
        intersectionlist=[]
        if options.multiple_negtive_winfiles!=[]:
            for n_inputfileName,threshold_title,outbedfilename in options.multiple_negtive_winfiles[:]:
                if options.anchorfile!=None:
                    n_inputfileName=standardseparately(options.anchorfile,n_inputfileName,int(options.upextend),int(options.downextend),int(options.winWidth),int(options.slideSize))
                threshold_title_list=re.split(r"_",threshold_title.strip())
                outfileNameWIN_Nlist.append(n_inputfileName)
                outfileNameWINwithGENE_Nlist.append((geneUtil.findTrscpt(n_inputfileName,outbedfilename, int(options.upextend), int(options.downextend), int(options.winWidth), int(options.slideSize), options.winType, "l", threshold_title_list[0], None, options.mergeNA, int(options.distalextend),options.trscptfound),threshold_title,outbedfilename))
                makeMhtGraph.makeHistonPicture(n_inputfileName, "Hp")#,"c(0,2000)","c(0,45)"
                makeMhtGraph.makeHistonPicture(outfileNameWINwithGENE_Nlist[-1][0], "Hp")#,"c(0,2000)","c(0,45)"
#                     print("awk 'NR!=1{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|grep --wFf """+options.removegenelistfile + """ - > """+outbedfilename+".trscptIDlist")
                print("awk 'NR!=1{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq>"""+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+options.minmaxSNP +"{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|sed '/^$/d'|awk '$0~/^ENS/{print $0}'>"""+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+options.minmaxSNP +"{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|sed '/^$/d'|awk '$0!~/^ENS/{print $0}'>"""+outbedfilename+".miRNAlist")
                if options.removegenelistfile!=None:
                    f=open(outbedfilename+".ENStrscptIDlist",'r')
                    mylist=[];
                    for line in f:

    #                     uniontpidlist.append(line.strip())
                        mylist.append(line.strip())
                    f.close()
                    ff=open(outbedfilename+".ENStrscptIDlist",'w')
                    removelist=[]    
############################
                    f=open(options.removegenelistfile,'r')
                    for line in f:
                        removelist.append(line.strip())
                    f.close()
#########################
                    finallist=list(set(mylist)-set(removelist))
                    removed+=list(set(mylist).intersection(set(removelist)))
                    for e in finallist:
                        print(e,file=ff)
                    ff.close()
                f=open(outbedfilename+".ENStrscptIDlist",'r')
                curset=set()
                mylist=[];
                for line in f:
                    curset.add(line.strip())
#                     uniontpidlist.append(line.strip())
                    intersectionlist.append(line.strip())
                    mylist.append(line.strip())
                f.close()
                removelist=[]    
                if options.removegenelistfile!=None:
                    f=open(options.removegenelistfile,'r')
                    
                    for line in f:
                        removelist.append(line.strip())
                    f.close()
                geneUtil.GOenrichment(options.gotablefile,outbedfilename,None,list(set(mylist)),None)
                if intertpidset:
                    intertpidset=intertpidset.intersection(curset)
                else:
                    intertpidset=curset
                print("grep -wFf "+outbedfilename+".ENStrscptIDlist"+""" /home/bioinfo/databases/ensembleIDconvert.txt|awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d'>"""+outbedfilename+""".Homologs_human""")
                os.system("grep -wFf "+outbedfilename+".ENStrscptIDlist"+""" /home/bioinfo/databases/ensembleIDconvert.txt|awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d'>"""+outbedfilename+""".Homologs_human""")
                print("grep -wFf "+outbedfilename+".Homologs_human /home/bioinfo/databases/humangenesymbl.txt|awk '{print $3}'|sort|uniq|sed '/^$/d'>"+outbedfilename+"Homologs_human_genesymbl")
                os.system("grep -wFf "+outbedfilename+".Homologs_human /home/bioinfo/databases/humangenesymbl.txt|awk '{print $3}'|sort|uniq|sed '/^$/d'>"+outbedfilename+"Homologs_human_genesymbl")
                print("""grep -wFf """+outbedfilename+""".Homologs_human /home/bioinfo/databases/humanGO.table |awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d' > """+outbedfilename+".Homologs_humanEntrezGeneID")
                os.system("""grep -wFf """+outbedfilename+""".Homologs_human /home/bioinfo/databases/humanGO.table |awk '{FS="\t";print $3}'|sort|uniq|sed '/^$/d' > """+outbedfilename+".Homologs_humanEntrezGeneID")            
        print("outfileNameWINwithGENE_Plist",outfileNameWINwithGENE_Plist)
        print("outfileNameWINwithGENE_Nlist",outfileNameWINwithGENE_Nlist)
        genelist,interlist=makeMhtGraph.makeMhtplots_compareInOnePicture_withgeneName(options.pathoutputfilename+".withgene", outfileNameWINwithGENE_Plist, outfileNameWINwithGENE_Nlist, 0,columnname)
        f=open(options.pathoutputfilename+".u."+str(len(genelist)),"w")
        for gene in genelist:
            print(gene,file=f)
        f.close()
        f=open(options.pathoutputfilename+".i."+str(len(interlist)),"w")
        for gene in interlist:
            print(gene,file=f)
        f.close()
        
        geneUtil.GOenrichment(options.gotablefile,options.pathoutputfilename+"u",None,list(set(list(intersectionlist)+uniontpidlist)),None)
        
        intersectionlist=set(intersectionlist).intersection(set(uniontpidlist))
        if intersectionlist !=set():
            geneUtil.GOenrichment(options.gotablefile,options.pathoutputfilename+"i",None,list(intersectionlist),None)
        splitinto=int(options.splitintopart)
    else:
        splitinto=int(options.splitintopart)
        outfileNameWIN_Nlist=[];outfileNameWIN_Plist=[];outfileNameWIN_Alist=[]
        for p_inputfileName,threshold,outbedfilename in options.multiple_positive_winfiles[:]:
            outfileNameWIN_Plist.append((p_inputfileName,threshold,outbedfilename))
        for n_inputfileName,threshold,outbedfilename in options.multiple_negtive_winfiles[:]:
            outfileNameWIN_Nlist.append((n_inputfileName,threshold,outbedfilename))
        for a_inputfileName,threshold,outbedfilename in options.multiple_allvalue_winfiles[:]:
            outfileNameWIN_Alist.append((a_inputfileName,threshold,outbedfilename))
        makeMhtGraph.makeMhtplots_compareInOnePicture(options.pathoutputfilename, outfileNameWIN_Plist, outfileNameWIN_Nlist,outfileNameWIN_Alist, 0,columnname,splitinto)
    print("removed genes:",removed)