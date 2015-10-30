# -*- coding: UTF-8 -*-
'''
Created on 2014-4-3

@author: liurui
'''
from NGS.BasicUtil import Util
from optparse import OptionParser
import pickle, sys, os, re, time
from Bio import Phylo,SeqIO
from io import StringIO

parser = OptionParser()
parser.add_option("-m", "--homologousgene", dest="homologousgene",
                  help="homologous file")
parser.add_option("-p", "--proteincdspair", dest="proteincdspairfile", help="proteincdspairfile")
parser.add_option("-c", "--configure", dest="configure")
parser.add_option("-l", "--minlen", dest="minlen")
parser.add_option("-A","--processA",dest="processA",action="store_true",default=False)
parser.add_option("-B","--processB",dest="processB",action="store_true",default=False)
parser.add_option("-C","--processC",dest="processC",action="store_true",default=False)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
minlen = int(options.minlen)
print("processA",options.processA,options.processB)
mytesttempfile=open("mytesttempfile_forprocessA.txt",'w')
chitesttable = {0.01:{2:9.2103,3:11.3449},
                0.05:{2:5.9915,3:7.8147}}

homogenefile = open(options.homologousgene, 'r')
aa_cds_pair_file = open(options.proteincdspairfile, 'r')
configure = open(options.configure, 'r')

cline = configure.readline()
tempPath = re.search(r'tempdir=(.*)', cline).group(1).strip()
cline = configure.readline()
musclePath = re.search(r'musclepath=(.*)', cline).group(1).strip()
cline = configure.readline()
pamlPath = re.search(r'pamlpath=(.*)', cline).group(1).strip()
cline = configure.readline()
outfileNamePre=re.search(r'outfilepre=(.*)', cline).group(1).strip()
cline = configure.readline()
PhyMLpath = re.search(r'PhyMLpath=(.*)', cline).group(1).strip()
# treefile_oringal = open(tempPath+"/", 'r')
# treefileName=(tempPath+"/treefileforonehomogene")
# temp_outtreefileName = treefileName + "_markbranch"


pamlcodeml = pamlPath + "/codeml"
pamlcodemlctl = pamlPath + "/codeml.ctl"
pamlcodemlctlfile = open(pamlcodemlctl, 'r')
MuscleInputFileName = tempPath + "/muscleinseq.fa"
MuscleOutputFileName = tempPath + "/muscleoutseqaln.aln"
pamlInputCDSFileName = tempPath + "/pamlinputfile.phy"
fastalnforcdstree=tempPath+"/pamlinputfile"
if __name__ == '__main__':

# load aa_cds_pair file and cdsfafile aafafile into memory
    aa_cds_filemap = {}
    """
        {species1:[aafafile,cdsfafile,aaindex,cdsfafile],species2:[,,,],,,,,}
    """
    for aa_cds_line in aa_cds_pair_file:
        speciesname = re.split(r':', aa_cds_line)[0].strip()
        aa_cds_pair = re.split(r':', aa_cds_line)[1]
        aa_cds_list = re.split(r';', aa_cds_pair)
        aafafileName = aa_cds_list[0].strip()
        cdsfafileName = aa_cds_list[1].strip()
        print("aafafileName",aafafileName,"cdsfafileName",cdsfafileName)
        if aafafileName != None and cdsfafileName != None:
            aa_cds_filemap[speciesname] = [open(aafafileName, 'r'), open(cdsfafileName, 'r')]
            aaindex = {}
            cdsindex = {}
            try:
                aa_cds_filemap[speciesname].append(pickle.load(open(aafafileName + ".myindex", 'rb')))
                aa_cds_filemap[speciesname].append(pickle.load(open(cdsfafileName + ".myindex", 'rb')))
            except IOError:
                print("generateIndexByChrom",speciesname)
                Util.generateIndexByChrom(aafafileName, aafafileName + ".myindex", "transcript:")
                Util.generateIndexByChrom(cdsfafileName, cdsfafileName + ".myindex")
                aa_cds_filemap[speciesname].append(pickle.load(open(aafafileName + ".myindex", 'rb')))
                aa_cds_filemap[speciesname].append(pickle.load(open(cdsfafileName + ".myindex", 'rb')))
            stat = os.system("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex")
            if stat != 0:
                print("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex" + " os.system return not 0")
                exit(-1)
            print("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex OK", stat)
# run muscle and paml loop
    finalkakslist = []
    """
    [(species1,species2,....,dn/ds,dn,ds),(pig,human,....,,,),,,,]
    
    """
    homotrscpttitle = re.split(r'~', homogenefile.readline())
    homotrscpttitle = [e.strip() for e in homotrscpttitle]
    finalkakslist.append(tuple(homotrscpttitle + ["dn/ds", "dn", "ds"]))
    processB_result_collection = {}
    
    processA_result_collection = {}
    skipthishomotrscptline = False  
    for homotrscptline in homogenefile:
        print("process:",homotrscptline)
        homotrscptlist = re.split(r'~', homotrscptline.strip())
        i = 0
        lenofhomeAA = []
        # make the aa fa file as the input of the muscle and run muscle        
        muscleinfile = open(MuscleInputFileName, 'w')
        firstofhomotrscpts = re.search(r'transcript:(.*)', homotrscptlist[0]).group(1).strip()
        processB_result_collection[firstofhomotrscpts] = []
        for trscpt in homotrscptlist:
            
            homotrscptlist[i] = re.search(r'transcript:(.*)', trscpt).group(1).strip()
            curspecies = homotrscpttitle[i].strip()
            print("collect aa ",trscpt,curspecies)
            if homotrscptlist[i] not in aa_cds_filemap[curspecies][2]:
                print(homotrscptlist[i],curspecies, "not in aa file or cds file")
                skipthishomotrscptline = True
                break
            aa_cds_filemap[curspecies][0].seek(aa_cds_filemap[curspecies][2][homotrscptlist[i]])
            # homotrscptlist[i]==currentChromNO
            print(trscpt,aa_cds_filemap[curspecies][2][homotrscptlist[i]],homotrscptlist[i])
            refSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(aa_cds_filemap[curspecies][0], homotrscptlist[i], mapname="transcript:")
            lenofhomeAA.append(len(refSeqMap[homotrscptlist[i]]))

            if "".join(refSeqMap[homotrscptlist[i]][1:-1]).find("*") == -1:
                if "".join(refSeqMap[homotrscptlist[i]][-1]) == "*":
                    refSeqMap[homotrscptlist[i]] = refSeqMap[homotrscptlist[i]][1:-1]
                print(">" + curspecies, file=muscleinfile)
                print("".join(refSeqMap[homotrscptlist[i]][1:]), file=muscleinfile)
            else:
                skipthishomotrscptline=True

            del refSeqMap[homotrscptlist[i]]
            i += 1
        muscleinfile.close()
        if skipthishomotrscptline:
            skipthishomotrscptline = False
            continue
        if max(lenofhomeAA) / min(lenofhomeAA) >= 2 or min(lenofhomeAA) < minlen:
            skipthishomotrscptline = True
        if skipthishomotrscptline:
            skipthishomotrscptline = False
            continue
        stat = os.system(musclePath + " -in " + MuscleInputFileName + " -out " + MuscleOutputFileName)
        if stat != 0 :
            print("Error:" + musclePath + " -in " + MuscleInputFileName + " -out " + MuscleOutputFileName)
            exit(-1)
        print("muscle ok:",stat)
        
        # fill back cds seq
#         muscleoutfile = open(MuscleOutputFileName, 'r')
        pamlinputcdsfile = open(pamlInputCDSFileName, 'w')
        pamlinputcdsheader = [];pamlinputcdsseq = []
        j = 0;curspecies = homotrscpttitle[j].strip()  # same order as specise
        
#         muscleout_seqmap = Util.decode_phyliplines(muscleoutfile)
        muscleout_seqmap={}
        muscleout_seqgenerator=SeqIO.parse(MuscleOutputFileName,"fasta")
        for seq_rec in muscleout_seqgenerator:
            muscleout_seqmap[seq_rec.id]=seq_rec.seq
#         muscleoutfile.close()
        maxlenlist=[]
        for species_and_trscpt_idx in range(len(homotrscpttitle)):
            curspecies = homotrscpttitle[species_and_trscpt_idx]
            curtrscpt = homotrscptlist[species_and_trscpt_idx]
            curtrscpt_idx_incdsfile = aa_cds_filemap[curspecies][3][curtrscpt]
            aa_cds_filemap[curspecies][1].seek(curtrscpt_idx_incdsfile)
            cdsSeqMap, cdscurrenttrscpt, cdsnexttrscpt = Util.getRefSeqMap(aa_cds_filemap[curspecies][1], currentChromNO=curtrscpt)
            locofaa = 0
            # fill back cds seq
            cdsseqfillback = ""
            for aa in muscleout_seqmap[curspecies][:]:
                if aa == "-":
                    codon = "---"
                else:
                    codon = "".join(cdsSeqMap[cdscurrenttrscpt][locofaa * 3 + 1:locofaa * 3 + 4])
                    locofaa += 1
                cdsseqfillback += codon
            pamlinputcdsheader.append(curspecies)
            maxlenlist.append(len(curspecies))
            pamlinputcdsseq.append(cdsseqfillback)
        maxlen=max(maxlenlist)
        print(maxlen)
        print(Util.encode_phyliplines(pamlinputcdsheader, pamlinputcdsseq,maxlen+2), file=pamlinputcdsfile)
#         ffff=open(fastalnforcdstree,'w')
#         for i in range(len(pamlinputcdsheader)):
#             print(">"+pamlinputcdsheader[i],file=ffff)
#             print(pamlinputcdsseq[i],file=ffff)
#         ffff.close()
#             
        pamlinputcdsfile.close()
        #make tree
        a=os.system(PhyMLpath+" -i "+pamlInputCDSFileName+" -m GTR -b 100 -t e -a e")
        if a!=0:
            print("error",PhyMLpath+" -i "+pamlInputCDSFileName+" -m GTR -b 100 -t e -a e")
            exit(-1)
#         os.system(PhyMLpath+" -infile="+pamlInputCDSFileName+" -type=DNA -output=FASTA -align")
        print(PhyMLpath+" -i "+pamlInputCDSFileName+" -m GTR -b 100 -t e -a e")
        treefile_oringal=open(pamlInputCDSFileName+"_phyml_tree.txt",'r')
        fixtreetext=treefile_oringal.readline()
        treefile_oringal.close()
        print(re.subn(r"\)[\d\.]+:","):",fixtreetext.strip())[0],end="",file=open(pamlInputCDSFileName+"_phyml_tree.txt",'w'))
        treefile_oringal=open(pamlInputCDSFileName+"_phyml_tree.txt",'r')
        temp_outtreefileName = pamlInputCDSFileName+"_phyml_tree.txt" + "_markbranch"
        if skipthishomotrscptline:
            skipthishomotrscptline = False
            continue
        # finishing fill back the cds seq file,next run codeml and extract ka ks value from mlc file
        
        ##### configure codeml.ctl file to process C
        if options.processC:
            pamlcodemlctlfilelines = pamlcodemlctlfile.readlines()
            pamlcodemlctlfile.close()
            os.system("rm " + pamlcodemlctl)
            pamlcodemlctlfile = open(pamlcodemlctl, 'w')
            for line in pamlcodemlctlfilelines:
                print(line)
                if re.search(r'^\s+seqfile\s*=', line) != None:  # seqfile
                    print("      seqfile = " + pamlInputCDSFileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+outfile\s*=', line) != None:
                    print("      outfile = " + tempPath + "/mlc", file=pamlcodemlctlfile)
                elif re.search(r'^\s+seqtype\s*=', line) != None:
                    print("      seqtype = 1", file=pamlcodemlctlfile)
                elif re.search(r'^\s+model\s*=', line) != None:
                    print("        model = 0", file=pamlcodemlctlfile)
                elif re.search(r'^\s+runmode\s*=', line) != None:
                    print("      runmode = -2", file=pamlcodemlctlfile)
                elif re.search(r'^\s+treefile\s*=', line) != None:
                    print("     treefile = stewart.trees      * tree structure file name",file=pamlcodemlctlfile)
                elif re.search(r'^\s+NSsites\s*=', line) != None:
                    print("      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;",file=pamlcodemlctlfile)
                elif re.search(r'^\s+fix_omega\s*=', line) != None:
                    print("    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate",file=pamlcodemlctlfile)
                elif re.search(r'^\s+omega\s*=',line)!=None:
                    print("        omega = .4 * initial or fixed omega, for codons or codon-based AAs",file=pamlcodemlctlfile)
                else:
                    print(line, file=pamlcodemlctlfile, end="")
            pamlcodemlctlfile.close()        
             
            print(pamlcodeml, pamlcodemlctl)
            stat = os.system(pamlcodeml)
            if stat != 0:
                print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                continue
    #             exit(-1)
            mlcfile = open(tempPath + "/mlc", 'r')
            mlclines = mlcfile.readlines()
            try:
                valuesabj = re.search(r'dN/dS=(.*)dN =(.*)dS =(.*)', mlclines[-1])
            except:
                print("may be cause this Error:Make sure to separate the sequence from its name by 2 or more spaces.")
                mlcfile.close()
                continue
            for line in mlclines:
                if re.search(r'^pairwise comparison,',line):
                    startcollect_pairwise_comparison=True
            dnds = valuesabj.group(1).strip()
            dn = valuesabj.group(2).strip()
            ds = valuesabj.group(3).strip()
            finalkakslist.append(tuple(homotrscptlist + [dnds, dn, ds]))
            mlcfile.close()
        
        
        treefile_oringal.seek(0)
        tree = Phylo.read(treefile_oringal, "newick")
        tree_terminal_list = tree.get_terminals()
        treefile_oringal.close()
        print(tree_terminal_list)      
        if options.processB:

            ##### configure codeml.ctl file to process B
            
            pamlcodemlctlfile = open(pamlcodemlctl, 'r')
            pamlcodemlctlfilelines = pamlcodemlctlfile.readlines()
            pamlcodemlctlfile.close()
            os.system("rm " + pamlcodemlctl)
            pamlcodemlctlfile = open(pamlcodemlctl, 'w')
            # model=2 nssite=0 runmode=0
            for line in pamlcodemlctlfilelines:
                if re.search(r'^\s+seqfile\s*=', line) != None:  # seqfile
                    print("      seqfile = " + pamlInputCDSFileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+outfile\s*=', line) != None:
                    print("      outfile = " + tempPath + "/mlc", file=pamlcodemlctlfile)
                elif re.search(r'^\s+seqtype\s*=', line) != None:
                    print("      seqtype = 1", file=pamlcodemlctlfile)
                elif re.search(r'^\s+model\s*=', line) != None:
                    print("        model = 2", file=pamlcodemlctlfile)
                elif re.search(r'^\s+runmode\s*=', line) != None:
                    print("      runmode = 0", file=pamlcodemlctlfile)
                elif re.search(r'^\s+treefile\s*=', line) != None:
                    print("     treefile = " + temp_outtreefileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+NSsites\s*=', line) != None:
                    print("      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;", file=pamlcodemlctlfile)
                elif re.search(r'^\s+fix_omega\s*=', line) != None:
                    print("    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate", file=pamlcodemlctlfile)
                elif re.search(r'^\s+omega\s*=', line) != None:
                    print("        omega = 1.5 * initial or fixed omega, for codons or codon-based AAs", file=pamlcodemlctlfile)
                else:
                    print(line, file=pamlcodemlctlfile, end="")
            pamlcodemlctlfile.close()
            for tree_terminal in tree_terminal_list:
                
                curspecies = tree_terminal.name
                species_and_trscpt_idx = homotrscpttitle.index(curspecies)
                curtrscpt = homotrscptlist[species_and_trscpt_idx]
                sio_f = StringIO()
                Phylo.write(tree, sio_f, "newick")
                sio_f.seek(0)
                temp_outtreefile = open(temp_outtreefileName, 'w')
                print(re.sub(r"" + curspecies, curspecies + "#1", sio_f.readline()), file=temp_outtreefile)
                temp_outtreefile.close();sio_f.close()
    
                stat = os.system(pamlcodeml)
                if stat != 0:
                    print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                    exit(-1)
                # extract data from mlc,fill data into processB_result_collection
                mlcfile = open(tempPath + "/mlc", 'r')
                mlclines = mlcfile.readlines()
                mlcline_idx=0
                while mlcline_idx < len(mlclines):
    #                 leafname_idx_map={}
                    if re.search(r"w \(dN/dS\) for branches:", mlclines[mlcline_idx]) != None:
                        branch_value_obj = re.search(r"w \(dN/dS\) for branches:\s*([\.\d]+)\s+([\.\d]+)", mlclines[mlcline_idx])
                        background_branch_value = branch_value_obj.group(1)
                        foreground_branch_value = branch_value_obj.group(2)
                        processB_result_collection[firstofhomotrscpts].append(foreground_branch_value)
#                         print(firstofhomotrscpts,processB_result_collection[firstofhomotrscpts],file=open("test.txt",'a'))
                        break
                    mlcline_idx+=1
    #                 if re.search(r"^tree length\s*=",mlclines[mlcline_idx])!=None:
    #                     mlcline_idx+=2
    #                     tree_codeleaf_text=re.sub(r"\s+","",mlclines[mlcline_idx])
    #                     mlcline_idx+=1
    #                     tree_orignalleafname_text=re.sub(r"\s+","",mlclines[mlcline_idx])
    #                     tree_codeleaf=Phylo.read(StringIO(tree_codeleaf_text),"newick")
    #                     tree_orignalleafname=Phylo.read(StringIO(tree_orignalleafname_text),"newick")
    #                     tree_orignalleafname_terminals=tree_orignalleafname.get_terminals()
    #                     tree_codeleaf_terminals=tree_codeleaf.get_terminals()
    #                     for terminal_idx in range(len(tree_codeleaf_terminals)):
    #                         leafname_idx_map[tree_orignalleafname_terminals[terminal_idx].name]=tree_codeleaf_terminals[terminal_idx].name
    #                         
                mlcfile.close()
        if options.processA:
            
            #########configure codeml.ctl file to process A ,first time collect LnL #####################
            pamlcodemlctlfile = open(pamlcodemlctl, 'r')
            pamlcodemlctlfilelines = pamlcodemlctlfile.readlines()
            pamlcodemlctlfile.close()
            os.system("rm " + pamlcodemlctl)
            pamlcodemlctlfile = open(pamlcodemlctl, 'w')
            # model=2 nssite=0 runmode=0
            for line in pamlcodemlctlfilelines:
                if re.search(r'^\s+seqfile\s*=', line) != None:  # seqfile
                    print("      seqfile = " + pamlInputCDSFileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+outfile\s*=', line) != None:
                    print("      outfile = " + tempPath + "/mlc", file=pamlcodemlctlfile)
                elif re.search(r'^\s+seqtype\s*=', line) != None:
                    print("      seqtype = 1", file=pamlcodemlctlfile)
                elif re.search(r'^\s+model\s*=', line) != None:
                    print("        model = 2", file=pamlcodemlctlfile)
                elif re.search(r'^\s+runmode\s*=', line) != None:
                    print("      runmode = 0", file=pamlcodemlctlfile)
                elif re.search(r'^\s+treefile\s*=', line) != None:
                    print("     treefile = " + temp_outtreefileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+NSsites\s*=', line) != None:
                    print("      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;", file=pamlcodemlctlfile)
                elif re.search(r'^\s+fix_omega\s*=', line) != None:
                    print("    fix_omega = 1  * 1: omega or omega_1 fixed, 0: estimate", file=pamlcodemlctlfile)
                elif re.search(r'^\s+omega\s*=', line) != None:
                    print("        omega = 1 * initial or fixed omega, for codons or codon-based AAs", file=pamlcodemlctlfile)
                else:
                    print(line, file=pamlcodemlctlfile, end="")
            pamlcodemlctlfile.close()
            processA_result_collection_lnL = {}
            for tree_terminal in tree_terminal_list:
                curspecies = tree_terminal.name
                species_and_trscpt_idx = homotrscpttitle.index(curspecies)
                curtrscpt = homotrscptlist[species_and_trscpt_idx]
                sio_f = StringIO()
                Phylo.write(tree, sio_f, "newick")
                sio_f.seek(0)
                temp_outtreefile = open(temp_outtreefileName, 'w')
                print(re.sub(r"" + curspecies, curspecies + "#1", sio_f.readline()), file=temp_outtreefile)
                temp_outtreefile.close();sio_f.close()
                stat = os.system(pamlcodeml)
                if stat != 0:
                    print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                    exit(-1)
                # extract data from mlc,fill data into processB_result_collection
                mlcfile = open(tempPath + "/mlc", 'r')
                mlclines = mlcfile.readlines()
                mlcline_idx=0
                while mlcline_idx < len(mlclines):
                    if re.search(r"^lnL", mlclines[mlcline_idx]) != None:
                        lnL_0 = re.search(r":\s+([-\.\d]+)",mlclines[mlcline_idx]).group(1)
                        processA_result_collection_lnL[curspecies] = float(lnL_0)
                    mlcline_idx+=1
    ######### configure codeml.ctl file to process A ,second time collect  #####################
            # model=2 nssite=0 runmode=0
            pamlcodemlctlfile = open(pamlcodemlctl, 'r')
            pamlcodemlctlfilelines = pamlcodemlctlfile.readlines()
            pamlcodemlctlfile.close()
            os.system("rm " + pamlcodemlctl)
            pamlcodemlctlfile = open(pamlcodemlctl, 'w')
            for line in pamlcodemlctlfilelines:
                if re.search(r'^\s+seqfile\s*=', line) != None:  # seqfile
                    print("      seqfile = " + pamlInputCDSFileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+outfile\s*=', line) != None:
                    print("      outfile = " + tempPath + "/mlc", file=pamlcodemlctlfile)
                elif re.search(r'^\s+seqtype\s*=', line) != None:
                    print("      seqtype = 1", file=pamlcodemlctlfile)
                elif re.search(r'^\s+model\s*=', line) != None:
                    print("        model = 2", file=pamlcodemlctlfile)
                elif re.search(r'^\s+runmode\s*=', line) != None:
                    print("      runmode = 0", file=pamlcodemlctlfile)
                elif re.search(r'^\s+treefile\s*=', line) != None:
                    print("     treefile = " + temp_outtreefileName, file=pamlcodemlctlfile)
                elif re.search(r'^\s+NSsites\s*=', line) != None:
                    print("      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;", file=pamlcodemlctlfile)
                elif re.search(r'^\s+fix_omega\s*=', line) != None:
                    print("    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate", file=pamlcodemlctlfile)
                elif re.search(r'^\s+omega\s*=', line) != None:
                    print("        omega = 1.5 * initial or fixed omega, for codons or codon-based AAs", file=pamlcodemlctlfile)
                else:
                    print(line, file=pamlcodemlctlfile, end="")
            pamlcodemlctlfile.close()
            processA_result_collection[firstofhomotrscpts]={}
            for tree_terminal in tree_terminal_list:
                significant=False
                curspecies = tree_terminal.name
                processA_result_collection[firstofhomotrscpts][curspecies]=[]
                species_and_trscpt_idx = homotrscpttitle.index(curspecies)
                curtrscpt = homotrscptlist[species_and_trscpt_idx]
                sio_f = StringIO()
                Phylo.write(tree, sio_f, "newick")
                sio_f.seek(0)
                temp_outtreefile = open(temp_outtreefileName, 'w')
                print(re.sub(r"" + curspecies, curspecies + "#1", sio_f.readline()), file=temp_outtreefile)
                temp_outtreefile.close();sio_f.close()
                stat = os.system(pamlcodeml)
                if stat != 0:
                    print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                    exit(-1)
                # extract data from mlc,fill data into processB_result_collection
                mlcfile = open(tempPath + "/mlc", 'r')
                mlclines = mlcfile.readlines()
                mlcfile.close
                mlcline_idx=0
                print(processA_result_collection)
                while mlcline_idx < len(mlclines):
                    if re.search(r"^lnL", mlclines[mlcline_idx]) != None:
                        lnL_1 = float(re.search(r":\s+([-\.\d]+)",mlclines[mlcline_idx]).group(1))
                        if (processA_result_collection_lnL[curspecies] - lnL_1) * 2>chitesttable[0.05][2]:
                            significant=True
                        else:
                            significant=False
                        print(processA_result_collection_lnL,significant,lnL_1,chitesttable[0.05][2])
                    elif significant and re.search(r"Bayes Empirical Bayes \(BEB\)",mlclines[mlcline_idx])!=None:
                        while True:
                            if re.search(r"\s+\d+\s+\w\s+[\.\d]+",mlclines[mlcline_idx])!=None or mlcline_idx==len(mlclines)-2:
                                break
                            mlcline_idx+=1
                        while  mlclines[mlcline_idx].split():
                            t3=re.search(r"\s+(\d+)\s+(\w)\s+([\.\d]+)([\*]*)",mlclines[mlcline_idx])#t3=re.split(r'\s+',mlclines[mlcline_idx].strip())
                            processA_result_collection[firstofhomotrscpts][curspecies].append((int(t3.group(1)),t3.group(2),float(t3.group(3)),t3.group(4)))
                            mlcline_idx+=1
                        else:
                            print(firstofhomotrscpts,processA_result_collection[firstofhomotrscpts],file=open("test.txt",'a'))
                            break
                    elif re.search(r"Bayes Empirical Bayes \(BEB\)",mlclines[mlcline_idx])!=None:#just for test
                        print(firstofhomotrscpts,curspecies,processA_result_collection[firstofhomotrscpts][curspecies],file=mytesttempfile)
                        print("significant judgement lnL:",lnL_1,"null hypothesis lnL:",processA_result_collection_lnL[curspecies],"(lnL_1 - lnL_1(null)) * 2",(processA_result_collection_lnL[curspecies] - lnL_1) * 2,"chitesttable",chitesttable[0.05][2],file=mytesttempfile)
                        while True:
                            if re.search(r"\s+\d+\s+\w\s+[\.\d]+",mlclines[mlcline_idx])!=None or mlcline_idx==len(mlclines)-2:
                                break
                            mlcline_idx+=1
                        while  mlclines[mlcline_idx].split():
                            print(mlclines[mlcline_idx])
                            t3=re.search(r"\s+(\d+)\s+(\w)\s+([\.\d]+)([\*]*)",mlclines[mlcline_idx])#  re.split(r'\s+',mlclines[mlcline_idx].strip())
                            print(int(t3.group(1)),t3.group(2),float(t3.group(3)),t3.group(4),file=mytesttempfile)
                            processA_result_collection[firstofhomotrscpts][curspecies].append((int(t3.group(1)),t3.group(2),float(t3.group(3)),t3.group(4)))
                            mlcline_idx+=1
                        else:
                            print(firstofhomotrscpts,processA_result_collection[firstofhomotrscpts],file=open("test.txt",'a'))
                            break                        
                    mlcline_idx+=1
    processB_outfile=open(outfileNamePre+"_branch",'w')
    #process A output map
    print("tpIDorGeneID",end="\t",file=processB_outfile)
    outfileMap={};specieslist=[]
    for tree_terminal in tree_terminal_list:
        curspecies = tree_terminal.name
        specieslist.append(curspecies)
        outfileMap[curspecies]= open(outfileNamePre+"_"+curspecies+"_branchsite",'w')
        print(curspecies+"kaks",end="\t",file=processB_outfile)
    for firstofhomotrscpts in sorted(processA_result_collection.keys()):
        for curspecies in processA_result_collection[firstofhomotrscpts].keys():
            print(firstofhomotrscpts,end=":",file=outfileMap[curspecies])
            for pos,animo,pro,sig in processA_result_collection[firstofhomotrscpts][curspecies]:
                print(str(pos),animo,str(pro),sig,sep=" ",end=";",file=outfileMap[curspecies])
            print(file=outfileMap[curspecies])
    for firstofhomotrscpts in sorted(processB_result_collection.keys()):
        print("\n",firstofhomotrscpts,end="\t",file=processB_outfile)
        for idx in range(len(processB_result_collection[firstofhomotrscpts])):
            print(processB_result_collection[firstofhomotrscpts][idx],end="\t",file=processB_outfile)
    processB_outfile.close()
    for s in outfileMap.keys():
        outfileMap[s].close()
#     for t in finalkakslist:
#         print("\t".join(t), file=outfile)
#     outfile.close()
#     treefile_oringal.close()
    configure.close()

    mytesttempfile.close()
    print("finish")
        
    
        

        
        
