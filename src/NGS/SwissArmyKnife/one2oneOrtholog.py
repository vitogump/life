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
import sys

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
    print("Warning:make sure single programm is running, or you may get wrong result because different program may access the same files")
    fileHmap={"Lox":0,"GRCh38":0,"MMUL":0,"PapA":0,"GSM":0,"Cjac":0,"CanF":0, "Pika":0,"PIG":0, "turT":0 , "GRCm":0, "UMD":0, "Equ":0, "Oar":0}
    groupfile=open(options.groupfile,"r")
    if options.steps.strip()=="1":
        
        of=open(options.groupfile+".step1out","w")
        for line in groupfile:
            #linelist=re.split(r"\s+",line.strip())
            if (line.count("Lox|")<2 or line.count("GRCh38|")<2 or line.count("MMUL|")<2 or line.count("PapA|")<=1 or  line.count("Cjac|")<=1 or line.count("CanF|")<=1 or line.count("Pika|")<=1 or line.count("PIG|")<=1 or line.count("turT|")<=1 or line.count("GRCm|")<=1 or line.count("Equ|")<=1) and ((line.count("UMD|")==line.count("GSM|")==line.count("Oar|")==1 ) or line.count("UMD|")==line.count("GSM|")==1 or line.count("Oar|")==line.count("GSM|")==1 ):#or line.count("UMD|")==line.count("GSM|")==1 or line.count("GSM|")==line.count("Oar|")==1 or line.count("Oar|")==line.count("UMD|")==1 #or line.count("UMD|")==1 or line.count("Oar|")==1 or  line.count("GSM|")==1
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
                #skip dup species that is not one-to-one ortholog
                if line.count(s_p.group(1)+"|")>=2:
                    print("skip specie",s_p.group(1)+"|")
                    continue

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
        from ete3 import  Tree
        cml=codeml.Codeml()
        "/home/lrui/paml4.9j/comdel.ctl"
        cml.read_ctl_file(os.path.join(options.pamlpath,"codeml.ctl"))
        
        
        #cml.ctl_file="codeml.ctl"
        #cml.working_dir="/home/lrui/paml4.9j"
        cml.set_options(seqtype=2)# 2:AAs;
        cml.set_options(model=3)#Empirical+F
        cml.set_options(alpha=0.1)#Empirical+F
        #cml.set_options(NSsites=[0])
        cml.set_options(ncatG=8)#seem not used
        cml.set_options(RateAncestor=1)
        cml.set_options(aaRatefile=os.path.join(options.pamlpath,"dat/wag.dat"))
        
        finalfile=open(options.groupfile+"".join(options.targetSpeciesname)+"step3.final",'w')
        tfinalascii_file=open(options.groupfile+"".join(options.targetSpeciesname)+"step3.finaltreefig",'w')        
        for fn in os.listdir(step2outpath):
            pathTofn=os.path.join(step2outpath,fn.strip())
            if os.path.isfile(pathTofn) and pathTofn.endswith("_align.fa.phy"):
                cml.alignment=pathTofn
                
                ########### test temp
                with open(pathTofn,'r') as pffffff:
                    physeqmap=Util.decode_phyliplines(pffffff)
                    print("# of species with one to one ortholog gene",len(physeqmap.keys()),len(fileHmap.keys()))
                    if len(physeqmap.keys())!=len(fileHmap.keys()):
                        print("extract subtree and write to subtree#of##.newick file")
                        tfffffffff=Tree(treefilefor14species,format=1)
                        tfffffffff.prune(list(physeqmap.keys()),preserve_branch_length=True)
                        tfffffffff.write(format=1, outfile="subtree"+str(len(physeqmap.keys()))+"of"+str(len(fileHmap.keys()))+".newick")

                        cml.tree="subtree"+str(len(physeqmap.keys()))+"of"+str(len(fileHmap.keys()))+".newick"
                    else:
                        cml.tree=treefilefor14species#"/home/lrui/specieslist14example_NAMED.nwk.txt"
                cml.write_ctl_file()
                #############test temp end
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
                    print("rstLine_idx",rstLine_idx)
                    if "tree with node labels for Rod Page's TreeView" in rstContent[rstLine_idx]:
                        while True:# for skip the empty line but It is not actually useful
                            rstLine_idx+=1;newickTreeString=rstContent[rstLine_idx].strip()
                            print("tree with node labels for Rod Page's TreeView\n",rstContent[rstLine_idx])
                            try:
                                sf = StringIO(rstContent[rstLine_idx])
                                treeWithNodeLables=Phylo.read(sf, "newick")#rstContent[rstLine_idx]
                                tlist=[];sf.close()
                            except:
                                print("search Tree Failed? ",rstContent[rstLine_idx],"continue searching the tree")
                                
                            for targetSpeciesName in options.targetSpeciesname:
                                tlist.append(list(treeWithNodeLables.find_clades(name=r"\d*_"+targetSpeciesName))[0])
                            break                            
                                
                        t_mrca = treeWithNodeLables.common_ancestor(*tlist)
                    if "List of extant and reconstructed sequences" in rstContent[rstLine_idx]:
                        sf=StringIO()
                        sf.write(rstContent[rstLine_idx+2]);nodes_AAs=re.split(r"\s+",rstContent[rstLine_idx+2].strip())# #nodes #bases
                        rstLine_idx+=4;
                        while re.search(r"[\w\W]{1,30}\s+([\w]{10}\s*)*",rstContent[rstLine_idx])!=None:
                            sf.write(rstContent[rstLine_idx].strip()+"\n")
                            rstLine_idx+=1
                        sf.seek(0)
                        """ no gap in align seq already"""
        
                        aln_anc_seqgenerator=AlignIO.parse(sf,'phylip')#<class 'generator'>
                        #seq_rec=Util.decode_phyliplines(sf.readlines())#can be used instead, but need modify the 10 truncate problem first.
                        
                        seq_rec=aln_anc_seqgenerator.__next__()
                        if len(treeWithNodeLables.get_terminals()+treeWithNodeLables.get_nonterminals())!=len(seq_rec):
                            print("Error somewhere,inconstant amount of tree not and seqs")
                        print("loading ancestrall stat/info finished")
                        parallel_sites={};convergent_sites={}
                        
#                         CD_sites_idx=[]
                        for i in range(int(nodes_AAs[1])):#for every AA position
                            #1. find  c==d sites i.e X3==X4
                            coverent_AA=""
                            pos_i_targetsAA=[seqrecord_obj.seq[i] for seqrecord_obj in seq_rec if seqrecord_obj.id in options.targetSpeciesname]# search each target species'AA for position i
                            if pos_i_targetsAA.count(pos_i_targetsAA[0])==len(options.targetSpeciesname):
                                print("test old version of what pass X3==X4",pos_i_targetsAA,"CD_sites_idx.append(i)",i,fn)                            
                            for AA in pos_i_targetsAA:
                                if pos_i_targetsAA.count(AA)>=2:
                                    print("find X3==X4",pos_i_targetsAA,AA)
                                    tlist_actual=[list(treeWithNodeLables.find_clades(name=r"\d*_"+seqrecord_obj.id))[0] for seqrecord_obj in seq_rec if ((seqrecord_obj.id in options.targetSpeciesname) and seqrecord_obj.seq[i]==AA)]
                                    print([(seqrecord_obj,seqrecord_obj.seq[i]) for seqrecord_obj in seq_rec if ((seqrecord_obj.id in options.targetSpeciesname) and seqrecord_obj.seq[i]==AA)])
                                    break# go step #2
                                    #CD_sites_idx.append(i)
                            else:
                                continue
                            print(tlist_actual)
                            #2. search changes along each target to root,A!=C,B!=D

                            Targets_ANCPATH_info={}
                            print(t_mrca,"may not same as ",treeWithNodeLables.common_ancestor(*tlist_actual),"when less than input","species")
                            t_mrca = treeWithNodeLables.common_ancestor(*tlist_actual)
                            
                            for t_clade in tlist_actual:
                                foroneTargetTmrca=treeWithNodeLables.trace(t_mrca,t_clade)
                                foroneTargetTmrca_ids=["node #"+str(ab.confidence) for ab in foroneTargetTmrca[:-1]]+[re.search(r"\d+_(\w+)",foroneTargetTmrca[-1].name).group(1)]#name_id same as code in  align seq in rst of paml
#检查 确认上一句 名字没有错误
                                pos_i_pathMrcaToTargetACorBD=[]
                                for seq_id in foroneTargetTmrca_ids:
                                    pos_i_pathMrcaToTargetACorBD.append([seqrecord_obj.seq[i] for seqrecord_obj in seq_rec if seqrecord_obj.id==seq_id][0])
                                #pos_i_pathMrcaToTargetACorBD=[seqrecord_obj.seq[i] for seqrecord_obj in seq_rec if seqrecord_obj.id in foroneTargetTmrca_ids]
                                print("check A != C || B!=D along path:",pos_i_pathMrcaToTargetACorBD,pos_i_pathMrcaToTargetACorBD.count(pos_i_pathMrcaToTargetACorBD[-1]),pos_i_pathMrcaToTargetACorBD[-1],pos_i_pathMrcaToTargetACorBD)
                                if pos_i_pathMrcaToTargetACorBD.count(pos_i_pathMrcaToTargetACorBD[-1])<len(pos_i_pathMrcaToTargetACorBD):#A is not equal to C
                                    print("get A != C || B!=D,record",t_clade,pos_i_pathMrcaToTargetACorBD,foroneTargetTmrca_ids)
                                    Targets_ANCPATH_info[t_clade.name]=[tuple(foroneTargetTmrca_ids)]+[*pos_i_pathMrcaToTargetACorBD]

                            print(foroneTargetTmrca_ids,"t_mrca",t_mrca,t_clade.name,type(t_clade.name),tlist_actual,type(tlist_actual[0]))
                            if len(Targets_ANCPATH_info.keys())==len(tlist_actual):
                                print("find A != C && B!=D for every target,ie # of records == # of targets; now judge A!=B, next")
                                """
                                path example
                                [Clade(confidence=16), Clade(confidence=22), Clade(confidence=24), Clade(confidence=25), Clade(confidence=26), Clade(name='12_GSM')]
                                """
                                #3. search for A!=B
                                print(Targets_ANCPATH_info)
                                print("Clades path for first species is :",Targets_ANCPATH_info[tlist_actual[0].name],"common ancestral AA:",Targets_ANCPATH_info[tlist_actual[0].name][1]," may same as target AA:",Targets_ANCPATH_info[tlist_actual[0].name][-1])
                                for X3T_idx in range(len(tlist_actual)):# each Terminal(Target) X3 
                                    X1_idx=-1
                                    for AA in reversed(Targets_ANCPATH_info[tlist_actual[X3T_idx].name][2:]):#not include the target and mrca AA, find the index of AA along PATH from X1  to mcra 
                                        print(Targets_ANCPATH_info[tlist_actual[X3T_idx].name].index(Targets_ANCPATH_info[tlist_actual[X3T_idx].name][-1]))
                                        if AA==Targets_ANCPATH_info[tlist_actual[X3T_idx].name][-1]:
                                            print(Targets_ANCPATH_info[tlist_actual[X3T_idx].name][2:],AA,"=?",Targets_ANCPATH_info[tlist_actual[X3T_idx].name][-1],X1_idx,"not this internal clade")
                                            X1_idx-=1;continue
                                        print("find change in Clade to",tlist_actual[X3T_idx].name)
                                        print(X1_idx,"find Clade  X1:"+Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0][X1_idx]+" in path to X3",tlist_actual[X3T_idx].name,AA)
                                        
                                        for X4T_Clade in tlist_actual[X3T_idx+1:]:# each Terminal(Target) X4 
                                            X2_idx=-1
                                            for AA2 in reversed(Targets_ANCPATH_info[X4T_Clade.name][2:]):# for one path, path to D  for example
                                                if AA2 ==Targets_ANCPATH_info[X4T_Clade.name][-1]:
                                                    X2_idx-=1;continue
                                                print("find Clade of X2 in path to X4",AA,"check X1!=X2")
                                                if Targets_ANCPATH_info[X4T_Clade.name][0][X2_idx]==Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0][X1_idx]:#may be 12 13 has the same internal path to 16 as the example in notebook9
                                                    print("skip to check other species: 12 13 has the same internal path to 16 as the example in notebook9",Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0][X1_idx],Targets_ANCPATH_info[X4T_Clade.name][0][X2_idx])
                                                    break
                                                print("get X1!=X2 Clade")
                                                if AA2!=AA:#
                                                    if i not in convergent_sites:convergent_sites[i]={}
                                                    print(" find convergent ",i,tlist_actual,Targets_ANCPATH_info.keys())
                                                    convergent_sites[i].update({tlist_actual[X3T_idx].name:Targets_ANCPATH_info[tlist_actual[X3T_idx].name][-1],X4T_Clade.name:Targets_ANCPATH_info[X4T_Clade.name][-1],Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0][X1_idx]:AA,Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0]:Targets_ANCPATH_info[tlist_actual[X3T_idx].name][1:],Targets_ANCPATH_info[X4T_Clade.name][0]:Targets_ANCPATH_info[X4T_Clade.name][1:]})#species1:X3_AA
                                                else:# all internal clades are same as the 
                                                    if i not in parallel_sites:parallel_sites[i]={} 
                                                    print("find parallel",i,tlist_actual,Targets_ANCPATH_info.keys())
                                                    parallel_sites[i].update({tlist_actual[X3T_idx].name:Targets_ANCPATH_info[tlist_actual[X3T_idx].name][-1],X4T_Clade.name:Targets_ANCPATH_info[X4T_Clade.name][-1],Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0][X1_idx]:AA,Targets_ANCPATH_info[tlist_actual[X3T_idx].name][0]:Targets_ANCPATH_info[tlist_actual[X3T_idx].name][1:],Targets_ANCPATH_info[X4T_Clade.name][0]:Targets_ANCPATH_info[X4T_Clade.name][1:]})
                                                print(" find convergent/ parallel",convergent_sites,parallel_sites)
                                                break
                                        break

#                                 for AA in Targets_ANCPATH_info[tlist[0].name][2:Targets_ANCPATH_info[tlist[0].name].index(Targets_ANCPATH_info[tlist[0].name][-1])]:#not include the target AA, find the index of AA along first target to mcra path 
#                                     print("get in AA",AA)
#                                     for D_targes in tlist[1:]:# sites for every other target species that convergent/parallel with this target would be recorded
#                                         for AA2 in Targets_ANCPATH_info[D_targes.name][2:Targets_ANCPATH_info[D_targes.name].index(Targets_ANCPATH_info[D_targes.name][-1])]:# for one path, path to D  for example
#                                             if AA2!=AA:
#                                                 print(" find convergent ")
#                                                 convergent_sites.append(i)
#                                                 break
#                                         else:
#                                             print("find parallel")
#                                             parallel_sites.append(i)
                            else:
                                print("skip site not changed in path C or D",i,Targets_ANCPATH_info.keys(),tlist)

                        
                        #print("treeWithNodeLables:",treeWithNodeLables)
                    rstLine_idx+=1        
                print("end while",parallel_sites,convergent_sites)
                if convergent_sites!={} or parallel_sites!={}:
                    print(pathTofn,newickTreeString,treeWithNodeLables,file=finalfile)
                    print("convergent sites:",file=finalfile);Util.mapFormatPrint(convergent_sites, finalfile);print("parallel_sits:",file=finalfile);Util.mapFormatPrint(parallel_sites,finalfile)
                    print(pathTofn,file=tfinalascii_file);Phylo.draw_ascii(treeWithNodeLables.root,tfinalascii_file)
                    sys.stdout.flush()
                #t_mrca = treeWithNodeLables.common_ancestor(*tlist).confidence
                
                  
        finalfile.close()
    print("done")
    groupfile.close()