
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
parser.add_option("-M", "--mapfile", dest="mapfile",default=None, help="optional")
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
# parser.add_option("-m", "--minmaxSNP", dest="minmaxSNP",default="0", help="upextend")
parser.add_option("-r", "--removegenelistfile", dest="removegenelistfile", help="upextend")
# parser.add_option("-t","--numberofoutlier_to_NearestGene",dest="numberofoutlier_to_NearestGene",default=0,help="number of outlier value,for example top 10")
(options, args) = parser.parse_args()
columnname=options.winType.strip()
minmaxSNP="0"#because the roh regions. in early selection . the candidate selection regions are naturally consider the minmaxSNP of 7
if __name__ == '__main__':
    makeMhtGraph = Make_Picture.MakeMhtGraph()
    outfileNameWINwithGENE_Plist=[];outfileNameWIN_Plist=[]
    outfileNameWINwithGENE_Nlist=[];outfileNameWIN_Nlist=[]
    outfileNameWINwithGENE_Alist_unfinished=[];outfileNameWIN_Alist=[]
    uniontpidlist=[];intertpidset=set()
    removed=[]
    if options.splitintopart==1:
        if options.multiple_positive_winfiles!=[]:
            for p_inputfileName,threshold_title,outbedfilename in options.multiple_positive_winfiles[:]:
#                 outfileNameWIN_Plist.append(p_inputfileName)
                threshold_title_list=re.split(r"_",threshold_title.strip())
                outfileNameWIN_Plist.append(((p_inputfileName+"arrangemented","none",outbedfilename)))
                outfileNameWINwithGENE_Plist.append((geneUtil.findTrscpt(p_inputfileName, outbedfilename, int(options.upextend), int(options.downextend), int(options.winWidth), int(options.slideSize), options.winType, "m", threshold_title_list, None, options.mergeNA, int(options.distalextend),options.anchorfile,options.trscptfound,options.mapfile),threshold_title,outbedfilename))
#                 makeMhtGraph.makeHistonPicture(p_inputfileName, "Fst")#,"c(0,2000)","c(0,45)"
                makeMhtGraph.makeHistonPicture(outfileNameWINwithGENE_Plist[-1][0], "Fst")
                print("awk 'NR!=1{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk'$0~/^ENS/{print $0}' >"+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+minmaxSNP +"{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk '$0~/^ENS/{print $0}' >"+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+minmaxSNP +"{print $8}' "+outbedfilename+".bed.selectedgene"+"|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|awk '$0!~/^ENS/{print $0}' >"+outbedfilename+".miRNAlist")
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

                threshold_title_list=re.split(r"_",threshold_title.strip())
                outfileNameWIN_Nlist.append((n_inputfileName+"arrangemented","none",outbedfilename))
                outfileNameWINwithGENE_Nlist.append((geneUtil.findTrscpt(n_inputfileName,outbedfilename, int(options.upextend), int(options.downextend), int(options.winWidth), int(options.slideSize), options.winType, "l", threshold_title_list, None, options.mergeNA, int(options.distalextend),options.anchorfile,options.trscptfound,options.mapfile),threshold_title,outbedfilename))
                makeMhtGraph.makeHistonPicture(n_inputfileName, "Hp")#,"c(0,2000)","c(0,45)"
                makeMhtGraph.makeHistonPicture(outfileNameWINwithGENE_Nlist[-1][0], "Hp")#,"c(0,2000)","c(0,45)"
#                     print("awk 'NR!=1{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|grep --wFf """+options.removegenelistfile + """ - > """+outbedfilename+".trscptIDlist")
                print("awk 'NR!=1{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq>"""+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+minmaxSNP +"{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|sed '/^$/d'|awk '$0~/^ENS/{print $0}'>"""+outbedfilename+".ENStrscptIDlist")
                os.system("awk 'NR!=1 && $7>="+minmaxSNP +"{print $8}' "+outbedfilename+""".bed.selectedgene"""+"""|sed 's/,/\\n/g' |sed  '/^$/d' |sort|uniq|sed '/^$/d'|awk '$0!~/^ENS/{print $0}'>"""+outbedfilename+".miRNAlist")
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
        if options.multiple_allvalue_winfiles!=[]:
            for a_inputfileName,threshold_title,outbedfilename in options.multiple_allvalue_winfiles[:]:
                threshold_title_list=re.split(r"_",threshold_title.strip())
                outfileNameWIN_Nlist.append((a_inputfileName+"arrangemented","none",outbedfilename))
                outfileNameWINwithGENE_Alist_unfinished.append((geneUtil.findTrscpt(a_inputfileName,outbedfilename, int(options.upextend), int(options.downextend), int(options.winWidth), int(options.slideSize), options.winType, "l", threshold_title_list, None, options.mergeNA, int(options.distalextend),options.anchorfile,options.trscptfound,options.mapfile),threshold_title,outbedfilename))
                
        print("outfileNameWINwithGENE_Plist",outfileNameWINwithGENE_Plist)
        print("outfileNameWINwithGENE_Nlist",outfileNameWINwithGENE_Nlist)
        for a_inputfileName,threshold,outbedfilename in options.multiple_allvalue_winfiles[:]:
            outfileNameWIN_Alist.append((a_inputfileName+"arrangemented","none",outbedfilename))
        makeMhtGraph.makeMhtplots_compareInOnePicture(options.pathoutputfilename, outfileNameWIN_Plist, outfileNameWIN_Nlist,outfileNameWIN_Alist, 0,columnname)
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