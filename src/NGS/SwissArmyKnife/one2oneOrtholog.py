'''
Created on 2022年3月28日

@author: RuiLiu
'''
from optparse import OptionParser
import os, re
import pickle

from NGS.BasicUtil import Util
from os.path import isfile


parser = OptionParser()
parser.add_option("-i", "--groupfile", dest="groupfile", help="input group file")
parser.add_option("-s", "--steps", dest="steps", help="1: selected assigend species 1-to-1 gene family into group file\n,2:extract seq.fa from previous out,3,4")
parser.add_option("-c", "--compliantFasta",dest="compliantFasta")
parser.add_option("-m", "--musclepath",dest="musclepath")
parser.add_option("-p", "--pamlpath",dest="pamlpath")
parser.add_option("-n", "--speciesname", dest="speciesname",action="append",default=[], help="")
(options, args) = parser.parse_args()
if __name__ == '__main__':
    groupfile=open(options.groupfile,"r")
    if options.steps.strip()=="1":
        
        of=open(options.groupfile+".out","w")
        for line in groupfile:
            #linelist=re.split(r"\s+",line.strip())
            if (line.count("Lox|")<2 or line.count("GRCh38|")<2 or line.count("MMUL|")<2 or line.count("PapA|")<=1 or  line.count("Cjac|")<=1 or line.count("CanF|")<=1 or line.count("Pika|")<=1 or line.count("PIG|")<=1 or line.count("turT|")<=1 or line.count("GRCm|")<=1 or line.count("Equ|")<=1) and ((line.count("UMD|")==line.count("GSM|")==line.count("Oar|")==1 ) ):#or line.count("UMD|")==line.count("GSM|")==1 or line.count("GSM|")==line.count("Oar|")==1 or line.count("Oar|")==line.count("UMD|")==1 #or line.count("UMD|")==1 or line.count("Oar|")==1 or  line.count("GSM|")==1
                if line.count("Lox|")+ line.count("GRCh38|") + line.count("MMUL|") + line.count("PapA|") +  line.count("Cjac|") + line.count("CanF|") + line.count("Pika|") + line.count("PIG|")+ line.count("turT|")+ line.count("GRCm|")+ line.count("Equ|")>0:
                    print(line,end="",file=of)
            else:
                print(line[:50])
        of.close()
    elif options.steps.strip()=="2":
        from Bio import Phylo,SeqIO
        
        try:
            os.mkdir(options.groupfile.replace(".","_"))
        except:
            pass
        fileHmap={"Lox":0,"GRCh38":0,"MMUL":0,"PapA":0,"GSM":0,"Cjac":0,"CanF":0, "Pika":0,"PIG":0, "turT":0 , "GRCm":0, "UMD":0, "Equ":0, "Oar":0}
        for fn in os.listdir(options.compliantFasta):
            pathTofn=os.path.join(options.compliantFasta.strip(),fn.strip())
            print("build/read indexfile for:",fn)
            if fn.replace(".fasta","") in fileHmap and os.path.isfile(pathTofn):
                try:
                    duckrefindex = pickle.load(open(os.path.join(options.compliantFasta.strip(),"."+fn.strip()) + ".myfasteridx", 'rb'))
            #             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
                except IOError:
                    Util.generateFasterRefIndex(pathTofn, os.path.join(options.compliantFasta.strip(),"."+fn.strip()) + ".myfasteridx",chrsignal="|")
                    duckrefindex = pickle.load(open(os.path.join(options.compliantFasta.strip(),"."+fn.strip()) + ".myfasteridx", 'rb'))                 
                fileHmap[fn.replace(".fasta","")]=[open(pathTofn,'r'),duckrefindex]
                print()
            else:
                print("pasd",fn,pathTofn)

        for line in groupfile:
            linelist= re.split(r"\s+",line.strip())
            #for every gene family
            genefamilyfaname=os.path.join(options.groupfile.replace(".","_"),linelist[0][:-1]+".fa")
            genefamilyalnname=os.path.join(options.groupfile.replace(".","_"),linelist[0][:-1]+"_align.fa")
            of_tmp=open(genefamilyfaname,"w")
            for nameOf_s_p in linelist[1:]:
                # genename and speciesname  for each specie
                s_p=re.search(r"(^[\d\w]+)\|(.*)$",nameOf_s_p)
                if s_p.group(1) not in fileHmap:
                    print("skip seq:",s_p.group(1),s_p.group(2))
                    continue
                #extract seq and print in one family file for each gene
                RefSeqMap=Util.getRefSeqBypos_faster(fileHmap[s_p.group(1)][0], fileHmap[s_p.group(1)][1], s_p.group(2), 1, 9999999999)
                print(">"+s_p.group(1)+"\n"+"".join(RefSeqMap[s_p.group(2)][1:]),file=of_tmp)
            of_tmp.close()
            #call muscle
            stat = os.system(options.musclepath.strip() + " -align " + genefamilyfaname + " -output " + genefamilyalnname)
            if stat!=0:
                print("Error muscle call");exit(-1)
            muscleout_seqgenerator=SeqIO.parse(genefamilyalnname,"fasta")
            species_names=[];aaseqs=[]
            for seq_rec in muscleout_seqgenerator:
                species_names.append(seq_rec.id)
                aaseqs.append("".join(seq_rec.seq))
            #trans aligned fa file into phylip for paml
            print(Util.encode_phyliplines(species_names, aaseqs,30), file=open(genefamilyalnname+".phy","w"))
            #config pamle and run condeml
            #stat=os.system(options.pamlpath.strip())
        for fn in fileHmap.keys():
            fileHmap[fn][0].close()
    elif options.steps=="3":
        step2outpath=options.groupfile.replace(".","_")
        treefilefor14species="/home/lrui/specieslist14example_NAMED.nwk.txt"
        from Bio.Phylo.PAML import codeml
        cml=codeml.Codeml()
        "/home/lrui/paml4.9j/comdel.ctl"
        cml.read_ctl_file(os.path.join(options.pamlpath,"codeml.ctl"))
        
        cml.tree=treefilefor14species#"/home/lrui/specieslist14example_NAMED.nwk.txt"
        #cml.ctl_file="codeml.ctl"
        #cml.working_dir="/home/lrui/paml4.9j"
        cml.set_options(seqtype=2)# 2:AAs;
        cml.set_options(model=3)#Empirical+F
        cml.set_options(alpha=0.1)#Empirical+F
        #cml.set_options(NSsites=[0])
        cml.set_options(ncatG=8)#seem not used
        cml.set_options(RateAncestor=1)
        cml.set_options(aaRatefile=os.path.join(options.pamlpath,"dat/wag.dat"))
        
        print("")        
        for fn in os.listdir(step2outpath):
            pathTofn=os.path.join(step2outpath,fn.strip())
            if os.path.isfile(pathTofn) and pathTofn.endswith("_align.fa.phy"):
                cml.alignment=pathTofn
                cml.write_ctl_file()
                stat = os.system(os.path.join(options.pamlpath,"codeml"))
                #copy rst and mcl and rename
                stat += os.system("cp rst "+os.path.join(step2outpath,fn.strip().replace("_align.fa.phy",".rst")))
                stat += os.system("cp mlc "+os.path.join(step2outpath,fn.strip().replace("_align.fa.phy",".mlc")))
                if stat!=0:
                    print("Error paml call");exit(-1)
                #extract 

    print("done")
    groupfile.close()