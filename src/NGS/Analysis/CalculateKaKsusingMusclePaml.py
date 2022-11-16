# -*- coding: UTF-8 -*-
'''
Created on 2014-4-3

@author: liurui
'''
import sys,inspect,os,configparser,random,string
import platform
cfparser = configparser.ConfigParser()

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
parser.add_option("-t", "--tempdir", dest="tempdir",default=os.getcwd())

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
minlen = int(options.minlen)


chitesttable = {0.01:{2:9.2103,3:11.3449},
                0.05:{2:5.9915,3:7.8147}}

homogenefile = open(options.homologousgene, 'r')
aa_cds_pair_file = open(options.proteincdspairfile, 'r')


cfparser = configparser.ConfigParser()
cfparser.read(options.configure, encoding='utf-8')

itemSoftwarecf=dict(cfparser.items("software"))
itemParametercf=dict(cfparser.items("parameter"))
print("allsections",cfparser.sections(),itemParametercf)
musclePath=cfparser.get("software","musclepath")
pamlPath=cfparser.get("software","pamlpath")
PhyMLpath=cfparser.get("software",'PhyMLpath')

if "fgspecies" in itemParametercf:
    print("fgspecies and outgroup must together")
    fgspecies,outgroup=cfparser.get("parameter",'fgspecies'),cfparser.get("parameter",'outgroup')
    print(fgspecies,outgroup)
else:
    fgspecies,outgroup=None,None

if 'pamltreefile' in itemParametercf:
    fixpamltree=cfparser.get("parameter",'pamltreefile')
    treefile_once=open(fixpamltree,'r')    
else:
    fixpamltree=False
ctlfile=cfparser.get("parameter",'ctltemplatefile')
outfileNamePre=cfparser.get("parameter",'outfilepre')
processType=cfparser.get("parameter","processType")
print("processType",processType,outfileNamePre,fixpamltree,processType)



tempPath = options.tempdir #re.search(r'tempdir=(.*)', cline).group(1).strip()

MuscleInputFileName = os.path.join(tempPath, "muscleinseq.fa")
MuscleOutputFileName = os.path.join(tempPath,"muscleoutseqaln.aln")
pamlInputCDSFileName = os.path.join(tempPath , "pamlinputfile.phy")
fastalnforcdstree=os.path.join(tempPath,"pamlinputfile")
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
                Util.generateIndexByChrom(cdsfafileName, cdsfafileName + ".myindex","transcript:")
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
    skipthishomotrscptline = False  
    outfileMap={};specieslist=[]
    if processType=="B":
        processB_result_collection = {}
        processB_outfile=open(outfileNamePre+"_branch",'w')
    elif processType=="A":
        processA_result_collection = {}
        
        mytesttempfile=open("mytesttempfile_forprocessA.txt",'w')
    print("tpIDorGeneID",end="\t",file=processB_outfile)
    fixtreetext=treefile_once.readline()
    treefile_once.close()
    if not fixpamltree or '#1' not in fixtreetext:
        print("set forebranch for each species")
        for sp in homotrscpttitle:
            print(sp+"fbranchw\tbgw",end="\t",file=processB_outfile)
            if processType=="A":
                specieslist.append(sp)
                outfileMap[sp]= open(outfileNamePre+"_"+sp+"_branchsite",'w')
    elif '#1' in fixtreetext:
        print("use seted froebranch for each homo gene group")
        print("fbranchw\tbgw",end="\t",file=processB_outfile)
            
        
    print("",file=processB_outfile)

    
    for homotrscptline in homogenefile:
        print("process:",homotrscptline)
        homotrscptlist = re.split(r'~', homotrscptline.strip())
        i = 0
        lenofhomeAA = []
        # make the aa fa file as the input of the muscle and run muscle        
        muscleinfile = open(MuscleInputFileName, 'w')
        firstofhomotrscpts = homotrscptlist[0]#re.search(r'transcript:(.*)', ).group(1).strip()
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
        pamlinputcdsphy = open(pamlInputCDSFileName, 'w')
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
                if aa == "-" or locofaa * 3 + 1>=len(cdsSeqMap[cdscurrenttrscpt]):# * STOP CODEN
                    codon = "---"
                else:
                    codon = "".join(cdsSeqMap[cdscurrenttrscpt][locofaa * 3 + 1:locofaa * 3 + 4])
                    locofaa += 1
                cdsseqfillback += codon
            pamlinputcdsheader.append(curspecies)
            maxlenlist.append(len(curspecies))
            pamlinputcdsseq.append(cdsseqfillback)
            sys.stdout.flush()
            #print(curspecies,len(cdsseqfillback),cdsSeqMap[cdscurrenttrscpt],len(cdsSeqMap[cdscurrenttrscpt]),locofaa * 3 + 1,locofaa * 3 + 4)
        maxlen=max(maxlenlist)
        print(Util.encode_phyliplines(pamlinputcdsheader, pamlinputcdsseq,maxlen+2), file=pamlinputcdsphy)
        ffff=open(fastalnforcdstree,'w')
        for i in range(len(pamlinputcdsheader)):
            print(">"+pamlinputcdsheader[i],file=ffff)
            print(pamlinputcdsseq[i],file=ffff)
        ffff.close()
#             
        pamlinputcdsphy.close()
        from Bio.Phylo.PAML import codeml
        cml=codeml.Codeml()
        if not fixpamltree:
            print('make tree and then assign to ctl file')
            a=os.system(PhyMLpath+" -i "+fastalnforcdstree+" -m GTR -b 100 -t e -a e")
            if a!=0:#below code block for clustalw func which is in case of the PhyML not work correct, may not used
                print("error",PhyMLpath+" -i "+pamlInputCDSFileName+" -m GTR -b 100 -t e -a e","\nplease make sure the 'clustalw=path' is indicated in the end of the configure file",)
                clustalw=cfparser.get("software",'clustalw')
                os.system(clustalw+" -infile="+fastalnforcdstree+" -type=DNA -output=FASTA -align")
                print(clustalw+" -infile="+pamlInputCDSFileName+" -type=DNA -output=FASTA -align")
                treefile_oringal=open(fastalnforcdstree+".dnd",'r')
                temp_outtreefileName = fastalnforcdstree+".dnd" + "_markbranch"
                #exit(-1)
            else:
                treefile_oringal=open(pamlInputCDSFileName+"_phyml_tree.txt",'r')
                fixtreetext=treefile_oringal.readline();treefile_oringal.close()
                print(re.subn(r"\)[\d\.]+:","):",fixtreetext.strip())[0],end="",file=open(pamlInputCDSFileName+"_phyml_tree.txt",'w'))
                treefile_oringal=open(pamlInputCDSFileName+"_phyml_tree.txt",'r')
                temp_outtreefileName = pamlInputCDSFileName+"_phyml_tree.txt" + "_markbranch"
    # #         os.system(PhyMLpath+" -infile="+pamlInputCDSFileName+" -type=DNA -output=FASTA -align")
    #         print(PhyMLpath+" -i "+pamlInputCDSFileName+" -m GTR -b 100 -t e -a e")
        else:
            print("use config tree file")
            temp_outtreefileName=fixpamltree
            treefile_oringal=open(fixpamltree,'r')
        if skipthishomotrscptline:
            skipthishomotrscptline = False
            continue
        # finishing fill back the cds seq file,next run codeml and extract ka ks value from mlc file
        pamlcodeml = os.path.join(pamlPath,"codeml")
        fixtreetext=treefile_oringal.readline()
        treefile_oringal.seek(0)
        tree = Phylo.read(treefile_oringal, "newick")
        tree_terminal_list = tree.get_terminals()
        treefile_oringal.close()
        print(tree_terminal_list)
        if processType.upper()=="C":
            cml.read_ctl_file(ctlfile)
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc")
            cml.set_options(seqtype=1);cml.set_options(model=0);cml.set_options(runmode=-2)

            cml.set_options(NSsites=0);cml.set_options(fix_omega=0);cml.set_options(omega = .4)
            cml.write_ctl_file()
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
        

              
        if processType.upper()=="B":
            cml.read_ctl_file(ctlfile)
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc")
            cml.tree=temp_outtreefileName
            cml.set_options(runmode=0);cml.set_options(model=2);cml.set_options(seqtype=1)
            cml.set_options(NSsites=[0]);cml.set_options(fix_omega=0);cml.set_options(omega=1.5)
            cml.write_ctl_file()
            ##### configure codeml.ctl file to process B

            # model=2 nssite=0 runmode=0
            
            if fixpamltree and (("#1" not in fixtreetext) and ('$1' not in fixtreetext)):
                print("use fixpamltree file:",fixpamltree)
                print("use fgspecies,outgroup information, if not exist report error")
            elif fixpamltree and fgspecies==outgroup==None:
                print("calculate kaks for the fg branch(s) according fixpamltree with assigned fg branch species marker #1 ")
                "get tree.terminal"
                "get homo species set"
                "现在还无法做 bg species 有空的或dup的，因为homogene file 已经限定每个物种一个基因"
                stat = os.system(pamlcodeml)
                if stat != 0:
                    print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                    exit(-1)
                processB_result_collection,dscollct=Util.extract_kaks_pamlmlc(tempPath + "/mlc",processB_result_collection,firstofhomotrscpts)

            elif not fixpamltree and fgspecies==outgroup==None:
                print("fixpamltree does not exist, construct tree using sequence, and then calculate kaks for each species as forebranch and other species as background branchs")
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
                    processB_result_collection,dscollct=Util.extract_kaks_pamlmlc(tempPath + "/mlc",processB_result_collection,firstofhomotrscpts)
            print(firstofhomotrscpts,end="\t",file=processB_outfile)
            for idx in range(len(processB_result_collection[firstofhomotrscpts])):
                print(*processB_result_collection[firstofhomotrscpts][idx],end="\t",file=processB_outfile)
            print("",file=processB_outfile)
            processB_outfile.flush()
#             print(bgbranchw,file=processB_outfile)
                    # extract data from mlc,fill data into processB_result_collection
                
                    
        if processType.upper()=="A":
            cml.read_ctl_file(ctlfile)
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc")
            cml.tree=temp_outtreefileName
            cml.set_options(runmode=0);cml.set_options(model=2);cml.set_options(seqtype=1)
            cml.set_options(NSsites=[2]);cml.set_options(fix_omega=1);cml.set_options(omega=1)
            cml.write_ctl_file()
            #########configure codeml.ctl file to process A ,first time collect LnL #####################

            # model=2 nssite=0 runmode=0


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
            cml.read_ctl_file(ctlfile)
            cml.set_options(runmode=0);cml.set_options(model=2);cml.set_options(seqtype=1)
            cml.set_options(NSsites=[2]);cml.set_options(fix_omega=0);cml.set_options(omega=1.5)
            cml.write_ctl_file()

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
            print()
        print("one home genes end:",homotrscptline)
        print("tpid","species1dnds","species2dnds","species1ds")
        
    
    if processType=="B":
        processB_outfile.close()


        
    
    #process A output map
    elif processType=="A":

        for firstofhomotrscpts in sorted(processA_result_collection.keys()):
            for curspecies in processA_result_collection[firstofhomotrscpts].keys():
                print(firstofhomotrscpts,end=":",file=outfileMap[curspecies])
                for pos,animo,pro,sig in processA_result_collection[firstofhomotrscpts][curspecies]:
                    print(str(pos),animo,str(pro),sig,sep=" ",end=";",file=outfileMap[curspecies])
                print(file=outfileMap[curspecies])
    
        for s in outfileMap.keys():
            outfileMap[s].close()
#     for t in finalkakslist:
#         print("\t".join(t), file=outfile)
#     outfile.close()
#     treefile_oringal.close()


        mytesttempfile.close()
    print("finish")
        
    
        

        
        
