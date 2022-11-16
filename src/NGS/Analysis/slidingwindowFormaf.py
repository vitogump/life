'''
Created on 2022年11月8日

@author: RuiLiu
'''
from optparse import OptionParser
import os
import subprocess

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
parser.add_option("-r","--refspeciesname",dest="refspeciesname")
parser.add_option("-C","--chromlistfilename",dest="chromlistfilename")
(options, args) = parser.parse_args()
logf=open("iqtree.log","wb")
stepsize=winsize=100
if __name__ == '__main__':
    #find the first pos and count how many species 
    
    for maln in AlignIO.parse(options.alignmentsMAF,"maf"):
        firstpos=-1;totalspecies=0#unused
        for seqrec in maln:
            print("=====================find the first alignments position from alignment block of maf file")
            if seqrec.id.startswith(options.refspeciesname) and int(seqrec.annotations['size'])>0:
                refchrSize=int(seqrec.annotations['srcSize'])
                firstpos=int(seqrec.annotations['start'])
                firstrefblocksize=seqrec.annotations['size']
                print("find first",seqrec,totalspecies,firstpos)
            totalspecies+=1
        if totalspecies>0 and firstpos>=0:
            break
    a=os.system("touch estimated_gene_trees.tree");print(a)
    curpos=firstpos
    #slide window
    idx=MafIO.MafIndex("test.mafindex",options.alignmentsMAF,options.refspeciesname)
    print("slide win and filter");wincount=0
    while curpos < refchrSize-stepsize and a==0:
        alnseqINwin=idx.get_spliced([curpos],[curpos+winsize], 1)
        seqslist=[]
        
        for seqrec in alnseqINwin:
            seqslist.append(str(seqrec.seq))
            if seqslist[-1].count("N")+seqslist[-1].count("-")>0.3*winsize:
                print("filtered 1) go curpos step forwards",curpos);break#   
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
                AlignIO.write(alnseqINwin, options.alignmentsMAF[:-4]+".tempwin.afa", "fasta")
                proc = subprocess.Popen(["iqtree","-s",options.alignmentsMAF[:-4]+".tempwin.afa","-m","MFP","-redo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                logf.write(stdout)
                logf.write(stderr)
                a=os.system("cat "+options.alignmentsMAF[:-4]+".tempwin.afa.treefile estimated_gene_trees.tree > win"+str(wincount))
                a+=os.system("mv win"+str(wincount)+" estimated_gene_trees.tree")  
                           
            
        wincount+=1
        curpos+=stepsize
    print(curpos,a)
    logf.cloes()