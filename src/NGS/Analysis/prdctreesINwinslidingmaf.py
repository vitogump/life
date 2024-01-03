'''
Created on 2022年11月8日

@author: RuiLiu
'''
from optparse import OptionParser
import os
import re,subprocess

from Bio import AlignIO
from Bio.AlignIO import MafIO


# for multiple_alignment in AlignIO.parse("BHG_muscovy_D2B_mallard_shaoxing2BJchr1.maf", "maf"):
#     for seqrec in multiple_alignment:
#         print(seqrec.annotations["start"])
#         print(seqrec.id)
#         if '4926'==seqrec.annotations["start"] or seqrec.annotations["start"]==4926:
#             break
#     else:
#         continue
#     break
# idx=MafIO.MafIndex("test.mafindex","BHG_muscovy_D2B_mallard_shaoxing2BJchr1.maf","ZJU1BJNC_051772.1")
# results=idx.search([4411,4511,4611],[4510,4610,4710])
# for maln in results:
#     for seqrec in maln:
#         print(seqrec)
# fafile=open("mytestmaf.slidebywin.fa","a")       
# maln=idx.get_spliced([4411,4511,4611], [4510,4610,4710], strand=1)
# AlignIO.write(maln, fafile, "fasta")    
parser = OptionParser()
#parser.add_option("-b","--bamlist_eachpop",dest="bamlist_eachpop",action="append",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-a","--alignmentsMAF",dest="alignmentsMAF",help="R(r)/G(g)")
parser.add_option("-r","--refspeciesname",dest="refspeciesname",help="tell program which refspeciesname in the maf used as ref to search interval")
parser.add_option("-q","--optionforiqtree",dest="optionforiqtree",action="append",default="",help="e.g. -q '-o outgroupname' -q '-pre /path/to/outpre'")
parser.add_option("-o","--outpath",dest="outpath",default="")
(options, args) = parser.parse_args()
logf=open("iqtree.log","wb")
stepsize=winsize=100
if (not os.path.isdir(options.outpath)) or (not options.outpath):
    print(options.outpath,"is not a dir");exit(-1)
if __name__ == '__main__':
    #find the first pos and count how many species 
    
    for maln in AlignIO.parse(options.alignmentsMAF,"maf"):
        firstpos=-1;totalspecies=0#unused
        for seqrec in maln:
            print("# find 'firstpos', the first alignments position from alignment block of maf file")
            if seqrec.id.startswith(options.refspeciesname) and int(seqrec.annotations['size'])>0:
                refchrSize=int(seqrec.annotations['srcSize'])
                firstpos=int(seqrec.annotations['start'])
                firstrefblocksize=seqrec.annotations['size']
                print("find first",seqrec,totalspecies,firstpos)
            totalspecies+=1
        if totalspecies>0 and firstpos>=0:
            print(totalspecies)
            break
    a=os.system("touch "+os.path.join(options.outpath,"estimated_gene_trees.tree"));print(a)
    curpos=firstpos
    
    print("#slide window and filter start from 'firstpos'",firstpos);wincount=0
    idx=MafIO.MafIndex(os.path.join(options.outpath,"test.mafindex"),options.alignmentsMAF,options.refspeciesname)
    while curpos < refchrSize-stepsize and a==0:
        alnseqINwin=idx.get_spliced([curpos],[curpos+winsize], 1)
        seqslist=[]
        
        for seqrec in alnseqINwin:
            seqslist.append(str(seqrec.seq))
            if seqslist[-1].count("N")+seqslist[-1].count("-")>0.3*winsize:
                print("filtered 1) go curpos step forwards, see P82",curpos);break#   
        else:
            samecount=0
            for pos in range(len(seqslist[0])):
                siteBases="".join([seq[pos] for seq in seqslist[:]])# [e for e in seqslist[:][pos]] is wrong
                if len(set(siteBases))==1 or siteBases.count("-")+siteBases.count("N")==len(siteBases):
                    samecount+=1
            if len(set(seqslist))<=4:
                print("filtered 2)  go curpos step forwards",curpos)# 
            elif samecount> 0.7 * winsize:
                print("filtered 3)  go curpos step forwards",curpos)# 
            else:
                AlignIO.write(alnseqINwin, os.path.join(options.outpath,options.alignmentsMAF[:-4]+".tempwin.afa"), "fasta")
                for optionforiqtree in options.optionforiqtree:
                    args=re.split(r"\s+",options.optionforiqtree)
                    iqtreeoutpre=args[-1] if "-pre" in optionforiqtree else os.path.join(options.outpath,options.alignmentsMAF[:-4]+".tempwin.afa")

                proc = subprocess.Popen(["iqtree","-s",options.alignmentsMAF[:-4]+".tempwin.afa","-m","MFP","-redo"]+args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                logf.write(stdout)
                logf.write(stderr)
                a=os.system("cat "+iqtreeoutpre+".treefile "+ os.path.join(options.outpath,"estimated_gene_trees.tree")+" > win"+str(wincount))
                a+=os.system("mv win"+str(wincount)+" "+os.path.join(options.outpath,"estimated_gene_trees.tree"))
                           
            
        wincount+=1
        curpos+=stepsize
    print(curpos,a)
    logf.cloes()
    print("then act command before run: java -jar astral.5.6.1.jar -i /home/lrui/temp/tttttt/estimated_gene_trees.tree -o estimated_species_tree.tree")
    print("sed 's/shaoxing[^:]*:/shaoxing:/g' /home/lrui/avian/duckevo/softmaskedTE/estimated_gene_trees.tree|sed 's/D2B[^:]*:/D2B:/g' |sed 's/anserBHG[^:]*:/anserBHG:/g' |sed 's/muscovy[^:]*:/muscovy:/g' |sed  's/mallard[^:]*:/mallard:/g' |sed 's/ZJU1BJ[^:]*:/ZJU1BJ:/g' > estimated_gene_trees.tree")