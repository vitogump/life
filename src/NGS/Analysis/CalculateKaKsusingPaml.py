# -*- coding: UTF-8 -*-
'''
Created on 2014-4-3

@author: liurui
'''
import re,inspect,os,configparser
import sys
from email.policy import default

from scipy.stats import chisquare,chi2
sys.path.append("/opt/life/src")
from bioinfodevelop.analysisUtils import Utils,genomicCodonTable
from optparse import OptionParser

from Bio import Phylo,SeqIO
from io import StringIO
from ete3 import  Tree

parser = OptionParser()
parser.add_option("-m", "--homologousgene", dest="homologousgene",
                  help="homologous file")
parser.add_option("-p", "--proteincdspair", dest="proteincdspairfile", help="proteincdspairfile")
parser.add_option("-c", "--configure", dest="configure")
parser.add_option("-l", "--minlen", dest="minlen")
parser.add_option("-t", "--tempdir", dest="tempdir",default=os.getcwd())
parser.add_option("-s", "--chrsignal",dest="chrsignal",default="transcript:",help="matchedanyway/>/chr/CHR...")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
minlen = int(options.minlen)


chitesttable = {0.01:{1:6.6349,2:9.2103,3:11.3449},
                0.05:{1:3.8415,2:5.9915,3:7.8147}}

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
    fgspecies,outgroup=re.split(r',',cfparser.get("parameter",'fgspecies')),re.split(r',',cfparser.get("parameter",'outgroup'))
    print(fgspecies,outgroup)
else:
    fgspecies,outgroup=None,None

if 'pamltreefile' in itemParametercf:
    fixpamltree=cfparser.get("parameter",'pamltreefile')
    treefile_once=open(fixpamltree,'r')    
else:
    fixpamltree=False
ctlfile=cfparser.get("parameter",'ctltemplatefile')
outfileNamePre=cfparser.get("parameter",'outfilepre')+re.search(r'[^/]+$',options.homologousgene).group().strip() if cfparser.get("parameter",'outfilepre') else os.path.join(os.path.abspath(os.curdir),re.search(r'[^/]+$',options.homologousgene).group().strip())
if os.path.isfile(outfileNamePre):
    print("check outname,exist will be covered")
    exit(-1)
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
        print("main()","aafafileName",aafafileName,"cdsfafileName",cdsfafileName)
        if aafafileName != None and cdsfafileName != None:
            aa_cds_filemap[speciesname] = [open(aafafileName, 'r'), open(cdsfafileName, 'r')]
            aaindex = {}
            cdsindex = {}
            try:
                aa_cds_filemap[speciesname].append(Utils.fastaidxer(aafafileName, options.chrsignal.strip()))#"matchedanyway"
                aa_cds_filemap[speciesname].append(Utils.fastaidxer(cdsfafileName, options.chrsignal.strip()))#"matchedanyway"
            except:
                print("mind pep/cds fa file, '> transcript:ID' format like this? or try add -s matchedanyway")
                exit(-1)
            stat = os.system("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex")
            if stat != 0:
                print("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex" + " Error")
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
        processB_outfile=open(outfileNamePre+".branch.kaks",'w')
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
                outfileMap[sp]= open(outfileNamePre+"_"+sp+".branchsite.kaks",'w')
    elif '#1' in fixtreetext:
        print("use seted froebranch for each homo gene group")
        print("fbranchw\tbgw",end="\t",file=processB_outfile)
            
        
    print("",file=processB_outfile)

    
    for homotrscptline in homogenefile:
        print("looping homogene:",homotrscptline)
        homotrscptlist = re.split(r'~', homotrscptline.strip())
        i = 0
        lenofhomeAA = []
        # make the aa fa file as the input of the muscle and run muscle        
        muscleinfile = open(MuscleInputFileName, 'w')
        firstofhomotrscpts = re.split(r'[;,]',homotrscptlist[0])[0].strip().strip("*").strip().replace("|transcript:",'_').replace("gene:","") #re.search(r'transcript:(.*)', ).group(1).strip()
        processB_result_collection[firstofhomotrscpts] = []
        for trscpt in homotrscptlist:
            homotrscpts=list(set(re.split(r'[;,]',trscpt.strip())))#for one species
#             homotrscptlist[i] = re.search(r'transcript:(.*)', trscpt).group(1).strip()
            homotrscptlist[i]=curhomotrscpt=re.split(r'[;,]',trscpt)[0].strip().strip("*") if re.search(r'transcript:(.*)', trscpt)==None else re.search(r'transcript:(.*)', trscpt).group(1).strip()
            curspecies = homotrscpttitle[i].strip()
            print("collect aa ",trscpt,curspecies)
            i += 1
            if curhomotrscpt.strip()=="NA" and curhomotrscpt not in aa_cds_filemap[curspecies][2]:#include NA
                print(curhomotrscpt,curspecies, "not in aa file or cds file")
                continue#skipthishomotrscptline = True
#                 break
            try:
                print(trscpt,aa_cds_filemap[curspecies][2][curhomotrscpt],curhomotrscpt)
            except KeyError:
                print("mind the format of -m homogenesfile, or the title of species name between -m and -p")
                exit(-1)
            refSeqMap, currentChromNO, nextChromNO = Utils.getRefSeqMap(aa_cds_filemap[curspecies][0],aa_cds_filemap[curspecies][2], curhomotrscpt)#, mapname=None"matchedanyway"
            lenofhomeAA.append(len(refSeqMap[curhomotrscpt]))

            if "".join(refSeqMap[curhomotrscpt][1:-1]).find("*") == -1:#stop coden
                if "".join(refSeqMap[curhomotrscpt][-1]) == "*":
                    refSeqMap[curhomotrscpt] = refSeqMap[curhomotrscpt][1:-1]
                print(">" + curspecies, file=muscleinfile)
                print("".join(refSeqMap[curhomotrscpt][1:]), file=muscleinfile)
            else:
                continue
                #skipthishomotrscptline=True

            del refSeqMap[curhomotrscpt]
            
        #out of trscpt travel of this homotrscptlist
        muscleinfile.close()
        if not lenofhomeAA or min(lenofhomeAA) < minlen or max(lenofhomeAA) / min(lenofhomeAA) >= 2:
            continue#skipthishomotrscptline = True

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
            curtrscpt = homotrscptlist[species_and_trscpt_idx].strip()
            if curspecies not in muscleout_seqmap: print(curspecies,curtrscpt,homotrscptlist);continue
            cdsSeqMap, cdscurrenttrscpt, cdsnexttrscpt = Utils.getRefSeqMap(aa_cds_filemap[curspecies][1],aa_cds_filemap[curspecies][3], currentChromNO=curtrscpt)
            locofaa = 0
            # fill back cds seq
            cdsseqfillback = ""
            for aa in muscleout_seqmap[curspecies][:]:
                if aa == "-" or locofaa * 3 + 1>=len(cdsSeqMap[cdscurrenttrscpt]):# * STOP CODEN
                    codon = "---"
                else:
                    c="".join(cdsSeqMap[cdscurrenttrscpt][locofaa * 3 + 1:locofaa * 3 + 4])
                    codon = c if c in genomicCodonTable and aa==genomicCodonTable[c.lower()] else c+"N"*(3-len(c))#"NNN"
                    locofaa += 1
                cdsseqfillback += codon #if len(codon)==3 else 
            pamlinputcdsheader.append(curspecies)
            maxlenlist.append(len(curspecies))
            pamlinputcdsseq.append(cdsseqfillback)
            sys.stdout.flush()
            #print(curspecies,len(cdsseqfillback),cdsSeqMap[cdscurrenttrscpt],len(cdsSeqMap[cdscurrenttrscpt]),locofaa * 3 + 1,locofaa * 3 + 4)
        maxlen=max(maxlenlist)
        print(Utils.encode_phyliplines(pamlinputcdsheader, pamlinputcdsseq,maxlen+2), file=pamlinputcdsphy)
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
            if set(fgspecies)& set(pamlinputcdsheader)!=set(fgspecies) or not(set(outgroup)&set(pamlinputcdsheader)):
                skipthishomotrscptline=True
#             print('make tree and then assign to ctl file')
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
#             print("use config tree file")
            temp_outtreefileName=fixpamltree
            treefile_oringal=open(fixpamltree,'r')

        # finishing fill back the cds seq file,next run codeml and extract ka ks value from mlc file
        pamlcodeml = os.path.join(pamlPath,"codeml")
        fixtreetext=treefile_oringal.readline()
        treefile_oringal.seek(0)
        tree = Phylo.read(treefile_oringal, "newick")
        tree_terminal_list = tree.get_terminals()
        treefile_oringal.close();fgcount=0;markedfgbranch1=[];markedfgbranch2=[]
        for tree_terminal in tree_terminal_list:
            if ('#1' in tree_terminal.name or '$1' in tree_terminal.name) :
                markedfgbranch1.append(tree_terminal.name);markedfgbranch2.append(re.sub(r"[#|$]1","",tree_terminal.name))
                if markedfgbranch2[-1] not in pamlinputcdsheader:skipthishomotrscptline=True;fgcount+=1
        if skipthishomotrscptline or len(pamlinputcdsheader)==fgcount:
            skipthishomotrscptline = False
            continue
        cml.read_ctl_file(ctlfile)
        if len(pamlinputcdsheader)<len(tree_terminal_list):
            tfffffffff=Tree(temp_outtreefileName,format=1);pamlinputcdsheadertmp=pamlinputcdsheader
            for tnameidx in range(len(markedfgbranch2)):pamlinputcdsheadertmp.remove(markedfgbranch2[tnameidx]);pamlinputcdsheadertmp.append(markedfgbranch1[tnameidx])
            tfffffffff.prune(pamlinputcdsheadertmp,preserve_branch_length=True)
            tfffffffff.write(format=1, outfile="subtree"+str(len(pamlinputcdsheader))+"of"+str(len(tree_terminal_list))+".newick")
            cml.tree="subtree"+str(len(pamlinputcdsheader))+"of"+str(len(tree_terminal_list))+".newick"#;print(pamlinputcdsheadertmp)
        else:
            cml.tree=temp_outtreefileName
        print(tree_terminal_list)
        def marktreeMKer(curspecies):
            species_and_trscpt_idx = homotrscpttitle.index(curspecies)
            curtrscpt = homotrscptlist[species_and_trscpt_idx]
            sio_f = StringIO();Phylo.write(tree, sio_f, "newick");sio_f.seek(0)
            temp_outtreefile = open(temp_outtreefileName, 'w')
            print(re.sub(r"" + curspecies, curspecies + "#1", sio_f.readline()), file=temp_outtreefile)
            temp_outtreefile.close();sio_f.close()
        if processType.upper()=="C":
#             cml.read_ctl_file(ctlfile)
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc"+firstofhomotrscpts+processType);mclfn="mlc"+firstofhomotrscpts+processType
            cml.set_options(seqtype=1);cml.set_options(model=0);cml.set_options(runmode=-2)

            cml.set_options(NSsites=0);cml.set_options(fix_omega=0);cml.set_options(omega = .4)
            cml.write_ctl_file()
            stat = os.system(pamlcodeml)
            if stat != 0:
                print("call paml maybe call this Error", pamlInputCDSFileName, "The seq file appears to be in fasta format, but not aligned?")
                continue
    #             exit(-1)
            mlcfile = open(os.path.join(tempPath,mclfn), 'r')
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
#             cml.read_ctl_file(ctlfile)
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc"+firstofhomotrscpts+processType);mclfn="mlc"+firstofhomotrscpts+processType
            
            cml.set_options(runmode=0);cml.set_options(model=2);cml.set_options(seqtype=1)
            cml.set_options(NSsites=[0]);cml.set_options(fix_omega=0);cml.set_options(omega=1.5)
            cml.write_ctl_file()
            ##### configure codeml.ctl file to process B
            # model=2 nssite=0 runmode=0
            extract_kaks_pamlmlcGenor=Utils.runpamlAccordingTree(pamlcodeml,processType,[fixpamltree,fixtreetext,tree_terminal_list],[fgspecies,outgroup],Utils.extract_kaks_pamlmlc,marktreeMKer)
            lnlfullist=[];LRTlist=[]
            for extract_kaks_pamlmlc in extract_kaks_pamlmlcGenor:# really execute code in runpamlAccordingTree(){...os.system('paml')...} and return a  function() named extract_kaks_pamlmlc
                processB_result_collection,lnl=extract_kaks_pamlmlc(os.path.join(tempPath,mclfn),processB_result_collection,firstofhomotrscpts,pmodel=2,pnssite=0)
            cml.set_options(model=0)
            extract_kaks_pamlmlcGenor=Utils.runpamlAccordingTree(pamlcodeml,processType,[fixpamltree,fixtreetext,tree_terminal_list],[fgspecies,outgroup],Utils.extract_kaks_pamlmlc,marktreeMKer)
            for extract_kaks_pamlmlc in extract_kaks_pamlmlcGenor:
                temp,lnlnull=extract_kaks_pamlmlc(os.path.join(tempPath,mclfn),processB_result_collection,firstofhomotrscpts,pmodel=0,pnssite=0)
            pvalue=chi2.sf(abs(2*(lnlnull-lnl)),1);statics=chisquare([lnl/(lnl+lnlnull),lnlnull/(lnl+lnlnull)], [0.5,0.5], 1)[0]
            
            print(firstofhomotrscpts,end="\t",file=processB_outfile)
            for idx in range(len(processB_result_collection[firstofhomotrscpts])):
                print(*processB_result_collection[firstofhomotrscpts][idx],sep="\t",end="\t",file=processB_outfile)
            print(pvalue,statics,abs(lnlnull - lnl) * 2>chitesttable[0.05][1],sep="\t",file=processB_outfile)
            processB_outfile.flush()
#             print(bgbranchw,file=processB_outfile)
                    # extract data from mlc,fill data into processB_result_collection
                
                    
        if processType.upper()=="A":
            
            cml.alignment=pamlInputCDSFileName
            cml.out_file=os.path.join(tempPath,"mlc"+firstofhomotrscpts+processType);mclfn="mlc"+firstofhomotrscpts+processType
            
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
                mlcfile = open(tempPath + "/"+mclfn, 'r')
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
                mlcfile = open(tempPath + "/"+mclfn, 'r')
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
            pass
        print("tpid","species1dnds","species2dnds","species1ds")
        print("one home genes end:",homotrscptline)
        
        
    
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
        
    
        

        
        
