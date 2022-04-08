'''
Created on 2022年3月28日

@author: RuiLiu
'''
"""
Each function accepts either a file name or an open file handle, so data can be also loaded from compressed files, StringIO objects, and so on. 
If the file name is passed as a string, the file is automatically closed when the function finishes; 
"""

import copy
from io import StringIO
from optparse import OptionParser
import os, re
from os.path import isfile
import pickle

from Bio import AlignIO

from NGS.BasicUtil import Util


parser = OptionParser()
parser.add_option("-i", "--groupfile", dest="groupfile", help="input group file")
parser.add_option("-s", "--steps", dest="steps", help="1: selected assigend species 1-to-1 gene family into group file\n,2:extract seq.fa from previous out,3,4")
parser.add_option("-c", "--compliantFasta",dest="compliantFasta")
parser.add_option("-m", "--musclepath",dest="musclepath")
parser.add_option("-p", "--pamlpath",dest="pamlpath")
parser.add_option("-t", "--targetSpeciesname", dest="targetSpeciesname",action="append",default=[], help="")
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
        from Bio import Phylo,SeqIO
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
        
        finalfile=open(options.groupfile+"".join(options.targetSpeciesname)+".final",'w')        
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
                rst=open(os.path.join(step2outpath,fn.strip().replace("_align.fa.phy",".rst")),'r')
                mlc=open(os.path.join(step2outpath,fn.strip().replace("_align.fa.phy",".mlc")),'r')
                rstContent=rst.readlines();rstLine_idx=0
                while rstLine_idx < len(rstContent):
                    
                    if "tree with node labels for Rod Page's TreeView" in rstContent[rstLine_idx]:
                        while True:# for skip the empty line but It is not actually useful
                            rstLine_idx+=1
                            try:
                                sf = StringIO(rstContent[rstLine_idx])
                                treeWithNodeLables=Phylo.read(sf, "newick")#rstContent[rstLine_idx]
                                tlist=[];sf.close()
                                for targetSpeciesName in options.targetSpeciesname:
                                    tlist.append(list(treeWithNodeLables.find_clades(name=r"\d*_"+targetSpeciesName))[0])
                                break
                            except:
                                print("search Tree Failed? ",rstContent[rstLine_idx],"continue searching the tree")
                                
                        t_mrca = treeWithNodeLables.common_ancestor(*tlist)
                    if "List of extant and reconstructed sequences" in rstContent[rstLine_idx]:
                        sf=StringIO()
                        sf.write(rstContent[rstLine_idx+2]);nodes_bases=re.split(r"\s+",rstContent[rstLine_idx+2].strip())# #nodes #bases
                        rstLine_idx+=4;
                        while re.search(r"[\w\W]{1,30}\s+([\w]{10}\s*)*",rstContent[rstLine_idx])!=None:
                            sf.write(rstContent[rstLine_idx].strip()+"\n")
                            rstLine_idx+=1
                        sf.seek(0)
                        """ no gap in align seq already"""
                        aln_anc_seqgenerator=AlignIO.parse(sf,'phylip')#<class 'generator'>
                        species_names=[];aaalignseqs=[]
                        seq_rec=aln_anc_seqgenerator.__next__()
                        if len(treeWithNodeLables.get_terminals()+treeWithNodeLables.get_nonterminals())!=len(seq_rec):
                            print("Error somewhere,inconstant amount of tree not and seqs")
                        print("loading ancestrall stat/info finished")
                        parallel_sites=[];convergent_sites=[]
                        #1. collect c==d
                        CD_sites_idx=[]
                        for i in range(int(nodes_bases[1])):
                            pos_i_targetsAA=[seqrecord_obj[i] for seqrecord_obj in seq_rec if seqrecord_obj.id in options.targetSpeciesname]
                            if pos_i_targetsAA.count(pos_i_targetsAA[0])==len(options.targetSpeciesname):
                                CD_sites_idx.append(i)
                        #2. search changes along each target to root,A!=C,B!=D
                        
                        for i in CD_sites_idx:
                            Targets_ANCPATH_info={}
                            for t_clade in tlist:
                                foroneTargetTmrca=treeWithNodeLables.trace(t_mrca,t_clade)
                                foroneTargetTmrca_ids=["node #"+ab.confidence for ab in foroneTargetTmrca[:-1]]+[re.search(r"\d+_(\w+)",foroneTargetTmrca[-1].name).group(1)]

                                pos_i_pathMrcaToTargetACorBD=[seqrecord_obj[i] for seqrecord_obj in seq_rec if seqrecord_obj.id in foroneTargetTmrca_ids]
                                if pos_i_pathMrcaToTargetACorBD.count(pos_i_pathMrcaToTargetACorBD[-1])<len(pos_i_pathMrcaToTargetACorBD):#A is not equal to C
                                    print("A != B || C!=D")
                                    Targets_ANCPATH_info[t_clade.id]=[tuple(foroneTargetTmrca_ids),copy.deepcopy(pos_i_pathMrcaToTargetACorBD)]


                            if len(Targets_ANCPATH_info.keys())==len(tlist):
                                print("find A != B || C!=D for every target; now judge A!=B, next")
                                """
                                path example
                                [Clade(confidence=16), Clade(confidence=22), Clade(confidence=24), Clade(confidence=25), Clade(confidence=26), Clade(name='12_GSM')]
                                """
                                #3. search for A!=B
                                for AA in Targets_ANCPATH_info[tlist[0].id][1:Targets_ANCPATH_info[tlist[0].id].index(Targets_ANCPATH_info[tlist[0].id][-1])]:#not include the target AA, find the index of first target AA along path 
                                    for D_targes in tlist[1:]:
                                        for AA2 in Targets_ANCPATH_info[D_targes.id][1:Targets_ANCPATH_info[D_targes.id].index(Targets_ANCPATH_info[D_targes.id][-1])]:# for one path, path to D  for example
                                            if AA2!=AA:
                                                print(" find convergent ")
                                                convergent_sites.append(i)
                                                break
                                        else:
                                            print("find parallel")
                                            parallel_sites.append(i)

                        
                        print("ssss")
                    rstLine_idx+=1        
                print("end while",parallel_sites,convergent_sites)
                if convergent_sites!=[] and parallel_sites!=[]:
                    print(pathTofn,":\nconvergent sites",convergent_sites,parallel_sites,file=finalfile)
                #t_mrca = treeWithNodeLables.common_ancestor(*tlist).confidence
                
                  
        finalfile.close()
    print("done")
    groupfile.close()