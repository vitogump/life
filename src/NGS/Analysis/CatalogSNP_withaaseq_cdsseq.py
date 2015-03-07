'''
Created on 2014-12-13

@author: liurui
'''
import copy
from optparse import OptionParser
import pickle
import re

from src.NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm

# testfile


parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa",
                  help="reference.fa")
parser.add_option("-g", "--gtffile", dest="gtffile", help="gtffile")
parser.add_option("-v", "--variantstable", dest="variantstable", help="variants")
parser.add_option("-b", "--bedfiles", dest="bedfiles", action="append", default=[], help="bedfiles")
parser.add_option("-o", "--outputpath", dest="outputpath", help="default infile1_infile2")
parser.add_option("-c", "--chromtablename", dest="chromtablename", help="")
parser.add_option("-d", "--chromdbname", dest="chromdbname", help="")
parser.add_option("-m", "--minlength", dest="minlength")
parser.add_option("-5", "--TSSregion", dest="TSSregion", default="0", help="")
parser.add_option("-3", "--utr3_region", dest="utr3_region", default="0")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
refFastaName = options.reffa
reffastaidxName = refFastaName + ".myindex"
reffahandler = open(options.reffa, "r")


minlength = options.minlength
chromtable = options.chromtablename
dbchromtools = dbm.DBTools("10.2.48.140", "root", "1234567", options.chromdbname)

variantstablename = options.variantstable.strip()
dbvariantstools = dbm.DBTools("10.2.48.140", "root", "1234567", "ninglabvariantdata")

gtfMap = Util.getGtfMap(options.gtffile)
bedfileNames = options.bedfiles

outputpath = options.outputpath.strip()


TSSregionlen = int(options.TSSregion)
utr3_region = int(options.utr3_region)
if outputpath[-1] != "/":
    outputpath = outputpath + "/"
mutaa = open(outputpath+variantstablename+".mutaa", 'w')
testmutcds = open(outputpath+variantstablename+".mutcds", 'w')
testrefaa = open(outputpath+variantstablename+".refaa", "w")
sql = "select * from " + chromtable + " where chrlength>=" + minlength
primaryID = "chrID"

intergenicVF = open(outputpath + "intergenic.Variantfile", 'w')
cdsVF = open(outputpath + "cds.Variantfile", 'w')
intronVF = open(outputpath + "intron.Variantfile", 'w')
utrVF = open(outputpath + "utr.Variantfile", 'w')
titlelist = [a[0].strip() for a in dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + variantstablename + "'")]
print(*(titlelist + ["trscptID", "geneID", "strand", "cdsidx", "refcodon", "refaa", "altcodon", "altaa"]), sep="\t", file=cdsVF)
print(*(titlelist + ["trscptID", "geneID", "strand", "intronidx"]), sep="\t", file=intronVF)
print(*(titlelist + ["trscptID", "geneID", "strand", "5'/3'"]), sep="\t", file=utrVF)
print(*titlelist, sep="\t", file=intergenicVF)
bedfileVFhandlerlist = []
for bedfile in bedfileNames:
    bedfileName = re.search(r'[^/]*$', bedfile).group(0)
    bedfileVFhandlerlist.append(open(outputpath + "bedfileName.Variantfile", "w"))
    print(*titlelist, sep="\t", file=bedfileVFhandlerlist[-1])

CodonTable = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
      'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
      'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
      'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
      'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
      'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
      'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
      'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
      'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
      'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
      'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
      'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
      'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
      'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
      'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
      'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
if __name__ == '__main__':
    genegrouptest = open("genegroup.txt", 'w')
    try:
        refidxByChr = pickle.load(open(reffastaidxName, 'rb'))
    except IOError:
        Util.generateIndexByChrom(refFastaName, reffastaidxName)
        refidxByChr = pickle.load(open(reffastaidxName, 'rb'))
        
    totalChroms = dbchromtools.operateDB("select", "select count(*) from " + chromtable + " where chrlength>=" + minlength)[0][0]
    for i in range(0, totalChroms, 20):
        currentsql = sql + " order by " + primaryID + " limit " + str(i) + ",20"
        result = dbchromtools.operateDB("select", currentsql)
        for row in result:
            seektuple = ()
            currentchrID = row[0]
            currentchrLen = int(row[1])
            if currentchrID not in gtfMap:
                snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(0) + " and snp_pos<" + str(currentchrLen) + " order by snp_pos")
                for snp in snps:
                    print(*snp, sep="\t", file=intergenicVF)            
                print("no gene in the chrom:", currentchrID)
                continue
            snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(0) + " and snp_pos<" + str(gtfMap[currentchrID][0][2]) + " order by snp_pos")
            for snp in snps:
                print(*snp, sep="\t", file=intergenicVF)
            GeneGrouplist = Util.getGeneGrouplist(gtfMap[currentchrID])
            print(currentchrID, "\n", len(GeneGrouplist), GeneGrouplist, file=genegrouptest)
            for geneGroup in GeneGrouplist:
                RefSeqMap = Util.getRefSeqBypos(reffahandler, refidxByChr, currentchrID, geneGroup[1][2], geneGroup[0], currentchrLen, seektuple)
                seektuple = (reffahandler.tell(), len(RefSeqMap[currentchrID]) + RefSeqMap[currentchrID][0] - 1)
                
                tscptSeqAllCds = {};tscptSeqAllCds_mut = {};cds_frame = {};mutat_amino_seq = {};ref_amino_seq = {}
                for gene in geneGroup[1:]:
                    tscptID = gene[0]
                    tscptSeqAllCds[tscptID] = []
                    cds_frame[tscptID] = {}  # {cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}
                    cdsidx = 3

                    for feature, elemStart, elemEnd, frame in gene[4:]:
                        cdsidx += 1
                        if feature == "CDS" or feature == "stop_codon":
                            cds_frame[tscptID][cdsidx] = (int(frame), len(tscptSeqAllCds[tscptID]))
                            tscptSeqAllCds[tscptID] += RefSeqMap[currentchrID][(elemStart - RefSeqMap[currentchrID][0]):(elemEnd - RefSeqMap[currentchrID][0] + 1)]
                    tscptSeqAllCds_mut[tscptID] = copy.deepcopy(tscptSeqAllCds[tscptID])                   
                    
                print("geneGroup", len(geneGroup), geneGroup, file=genegrouptest)
                for gene_idx in range(1, len(geneGroup)):
                    tscptID = geneGroup[gene_idx][0]
                    print(tscptID, file=genegrouptest)
                    linetoCDSMap = {}
                    linetoIntronMap = {}
                    if geneGroup[gene_idx][1] == '+':
                        snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]) + " and snp_pos<" + str(geneGroup[gene_idx][3]) + " order by snp_pos")  # extend UTR here
                        for snp in snps:
                            snppos = snp[1]
                            refbase = snp[3]
                            altbase = snp[4]
                            cdsidx = 3
                            Intron_idx = -1
                            if re.search(r'[^a-zA-Z]', altbase) != None:  # contain ',' ie. multiple alle
                                continue  # go to the next snp
                            for feature, elemStart, elemEnd, frame in geneGroup[gene_idx][4:]:
                                
                                cdsidx += 1
                                if snppos <= elemEnd and snppos >= elemStart:
                                    if feature == 'CDS' or feature == 'stop_codon' :
                                        
                                        if snppos + len(refbase) - 1 > elemEnd or snppos + len(altbase) - 1 > elemEnd:
                                            print(snp, tscptID, "indelatEdgeofCDS")
                                            break  # go to the next snp
###################################################
                                        if len(refbase) > len(altbase):  # situation TAA     TA;ACG     A
                                            tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1]):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(altbase))] = list(altbase)
                                            tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(altbase)):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase))] = [' '] * (len(refbase) - len(altbase))
                                        elif len(refbase) < len(altbase):  # situation TTA     TTAAACTTCTATACTA;T       TATA;
                                            if len(refbase) == 1:
                                                tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1]] = len(altbase)
                                            else:
                                                tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1]):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase) - 1)] = list(altbase[0:(len(refbase) - 1)])
                                                tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase) - 1] = altbase[(len(refbase) - 1):]
                                        
                                        else:  # len(refbase)==len(altbase)==1
                                            tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1]] = altbase
                                            if snp[:5] in linetoCDSMap:
                                                linetoCDSMap[snp[:5]][1] = linetoCDSMap[snp[:5]][1] + ";" + tscptID;linetoCDSMap[snp[:5]][2] = linetoCDSMap[snp[:5]][2] + ";geneID";linetoCDSMap[snp[:5]][3] = linetoCDSMap[snp[:5]][3] + ";+";linetoCDSMap[snp[:5]][4] = linetoCDSMap[snp[:5]][4] + ";" + str(cdsidx-3)
                                                print("mytest",linetoCDSMap[snp[:5]])
                                            else:
                                                linetoCDSMap[snp[:5]] = [snp[5:], tscptID, "geneID", "+", str(cdsidx-3)]
###########################################################################
                                elif snppos > elemEnd and snppos < geneGroup[gene_idx][cdsidx + 1][1]:
                                    Intron_idx = cdsidx - 3
                            if Intron_idx != -1:
                                print(*(list(snp) + [tscptID, "geneID", "+", Intron_idx]), sep="\t", file=intronVF)
####################           translate protein       ###################
                        gene=geneGroup[gene_idx]
                        mutat_amino_seq[tscptID] = []
                        ref_amino_seq[tscptID] = []
                        tscptSeqAllCds_mut_str = "".join(filter(lambda e:e.strip() != "", tscptSeqAllCds_mut[tscptID]))
#                             tscptSeqAllCds_mut[tscptID] = list(tscptSeqAllCds_mut_str)
                        for i in range(cds_frame[tscptID][sorted(cds_frame[tscptID].keys())[0]][0], len(tscptSeqAllCds_mut_str), 3):  # produce mutat_amino_seq
                            codon_m = tscptSeqAllCds_mut_str[i:i + 3].lower()
                            try:
                                mutat_amino_seq[tscptID].append(CodonTable[codon_m])
                            except KeyError:
                                mutat_amino_seq[tscptID].append('X')
                        for i in range(cds_frame[tscptID][sorted(cds_frame[tscptID].keys())[0]][0], len(tscptSeqAllCds[tscptID]), 3):  # produce linetoCDSMap
                            codon = "".join(tscptSeqAllCds[tscptID][i:i + 3]).lower()
                            codon_m = "".join(tscptSeqAllCds_mut[tscptID][i:i + 3]).lower()
                            if codon != codon_m:
                                try:
                                    snppos_cds, ref_base_cds, alt_base_cds = Util.getSNPrecInCDS(i, len(tscptSeqAllCds[tscptID]), codon, codon_m, cds_frame[tscptID], gene)
                                    if (currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)])==10:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] + ";" + CodonTable[codon];linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] + ";" + CodonTable[codon_m]
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)] += [codon, CodonTable[codon], codon_m, CodonTable[codon_m]]
                                    else:
                                        print(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                except KeyError:
                                    snppos_cds, ref_base_cds, alt_base_cds = Util.getSNPrecInCDS(i, len(tscptSeqAllCds[tscptID]), codon, codon_m, cds_frame[tscptID], gene)
                                    if (currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)])==10:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] + ";X";linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] + ";X"
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)] += [codon, "X", codon_m, "X"]
                                    else:
                                        print(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                                                        
                            try:
                                ref_amino_seq[tscptID].append(CodonTable["".join(tscptSeqAllCds[tscptID][i:i + 3]).lower()])
                            except KeyError:
                                ref_amino_seq[tscptID].append("X")
                                    
                    else:  # strand == '-'
                        snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]) + " and snp_pos<" + str(geneGroup[gene_idx][3]) + " order by snp_pos")
                        for snp in snps:
                            snppos = snp[1]
                            refbase = snp[3]
                            altbase = snp[4]
                            cdsidx = 3
                            Intron_idx = -1
                            if re.search(r'[^a-zA-Z]', altbase) != None:  # contain ',' ie. multiple alle
                                continue  # go to the next snp
                            for feature, elemStart, elemEnd, frame in geneGroup[gene_idx][4:]:
                                
                                cdsidx += 1
                                if snppos <= elemEnd and snppos >= elemStart:
                                    if feature == 'CDS' or feature == 'stop_codon':
                                        
                                        if snppos + len(refbase) - 1 > elemEnd or snppos + len(altbase) - 1 > elemEnd:
                                            print(snp, tscptID, "indelatEdgeofCDS")
                                            break  # go to the next snp
################################################################
                                        if len(refbase) > len(altbase):  # situation TAA     TA;ACG     A
                                            tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1]):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(altbase))] = list(altbase)
                                            tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(altbase)):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase))] = [' '] * (len(refbase) - len(altbase))
                                        elif len(refbase) < len(altbase):  # situation TTA     TTAAACTTCTATACTA;T       TATA;
                                            if len(refbase) == 1:
                                                tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1]] = len(altbase)
                                            else:
                                                tscptSeqAllCds_mut[tscptID][(snppos - elemStart + cds_frame[tscptID][cdsidx][1]):(snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase) - 1)] = list(altbase[0:(len(refbase) - 1)])
                                                tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1] + len(refbase) - 1] = altbase[(len(refbase) - 1):]
                                        
                                        else:  # len(refbase)==len(altbase)==1
                                            tscptSeqAllCds_mut[tscptID][snppos - elemStart + cds_frame[tscptID][cdsidx][1]] = altbase
                                            print("1",snp[:5],linetoCDSMap,file=open("debug",'a'))
                                            if snp[:5] in linetoCDSMap:
                                                linetoCDSMap[snp[:5]][1] = linetoCDSMap[snp[:5]][1] + ";" + tscptID;linetoCDSMap[snp[:5]][2] = linetoCDSMap[snp[:5]][2] + ";geneID";linetoCDSMap[snp[:5]][3] = linetoCDSMap[snp[:5]][3] + ";-";linetoCDSMap[snp[:5]][4] = linetoCDSMap[snp[:5]][4] + ";" + str(len(cds_frame[tscptID]) - (cdsidx-4))
                                            else:
                                                linetoCDSMap[snp[:5]] = [snp[5:],tscptID, "geneID", "-", str(len(cds_frame[tscptID]) - (cdsidx-4))]   
########################################################
                                elif snppos > elemEnd and snppos < geneGroup[gene_idx][cdsidx + 1][1]:
                                    Intron_idx = len(cds_frame[tscptID]) - (cdsidx-3)
                            if Intron_idx != -1:
                                print(*(list(snp) + [tscptID, "geneID", "-", Intron_idx]), sep="\t", file=intronVF)
# translate protein
                        gene=geneGroup[gene_idx]
                        mutat_amino_seq[tscptID] = []
                        ref_amino_seq[tscptID] = []
                        tscptSeqAllCds_mut_str = "".join(filter(lambda e:e.strip() != "", tscptSeqAllCds_mut[tscptID]))
                        tscptSeqAllCds_mut_str = Util.complementary(tscptSeqAllCds_mut_str)
                        tscptSeqAllCds_mut_str = tscptSeqAllCds_mut_str[::-1]
#                             tscptSeqAllCds_mut[tscptID] = list(tscptSeqAllCds_mut_str)
#############################      produce linetoCDSMap  ############
                        tscptSeqAllCds_Revr_Cmplm = Util.complementary(tscptSeqAllCds[tscptID])
                        tscptSeqAllCds_Revr_Cmplm.reverse()
                        tscptSeqAllCds_mut_Revr_Cmplm = Util.complementary(tscptSeqAllCds_mut[tscptID])
                        tscptSeqAllCds_mut_Revr_Cmplm.reverse()
                        for bases_idx in range(len(tscptSeqAllCds_mut_Revr_Cmplm)):  # reverse every element of the list,ie . reverse every str of the list
                            tscptSeqAllCds_mut_Revr_Cmplm[bases_idx] = tscptSeqAllCds_mut_Revr_Cmplm[bases_idx][::-1]
                        if len(tscptSeqAllCds_Revr_Cmplm) != len(tscptSeqAllCds_mut_Revr_Cmplm):
                            print("length should be equal")
                            exit(-1)
                        for i in range(cds_frame[tscptID][sorted(cds_frame[tscptID].keys())[-1]][0], len(tscptSeqAllCds_Revr_Cmplm), 3):
                            codon = "".join(tscptSeqAllCds_Revr_Cmplm[i:i + 3]).lower()
                            codon_m = "".join(tscptSeqAllCds_mut_Revr_Cmplm[i:i + 3]).lower()
                            
                            if codon != codon_m:
                                try:
                                    snppos_cds, ref_base_cds, alt_base_cds = Util.getSNPrecInCDS(i, len(tscptSeqAllCds[tscptID]), codon, codon_m, cds_frame[tscptID], gene)
                                    if (currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        print((currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds),"in",linetoCDSMap,file=open("debug",'a'))
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)])==10:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] + ";" + CodonTable[codon];linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] + ";" + CodonTable[codon_m]
                                            print("mytest",linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)],file=open("debug",'a'))
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)] += [codon, CodonTable[codon], codon_m, CodonTable[codon_m]]
                                    else:
                                        print((currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds),"not in",linetoCDSMap,file=open("debug",'a'))
                                        print(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                except KeyError:
                                    print("except KEYERROR")
                                    if (currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)])==10:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][6] + ";X";linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)][8] + ";X"
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds)] += [codon, "X", codon_m, "X"]
                                    else:
                                        print(currentchrID,snppos_cds,".",ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",linetoCDSMap,tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                            try:
                                mutat_amino_seq[tscptID].append(CodonTable[codon_m])
                                ref_amino_seq[tscptID].append(CodonTable["".join(tscptSeqAllCds_Revr_Cmplm[i:i + 3]).lower()])
                            except KeyError:
                                ref_amino_seq[tscptID].append("X")
                            
                    print(">transcript:" + tscptID, file=testmutcds)
                    print(">" + tscptID, file=testrefaa)
                    print(">" + tscptID, file=mutaa)
                    k = 0
                    cdsstrline = "".join(tscptSeqAllCds_mut[tscptID][k:k + 60])
                    while len(cdsstrline) == 60:
                        print(cdsstrline, file=testmutcds);k += 60
                        cdsstrline = "".join(tscptSeqAllCds_mut[tscptID][k:k + 60])
                    else:
                        print(cdsstrline, file=testmutcds)
                    k = 0
                    aastrline = "".join(mutat_amino_seq[tscptID][k:k + 60])
                    while len(aastrline)==60:
                        print(aastrline, end="\n", file=mutaa);k += 60
                        aastrline = "".join(mutat_amino_seq[tscptID][k:k + 60])
                    else:
                        print(aastrline, file=mutaa)    
                    k = 0
                    aastrline = "".join(ref_amino_seq[tscptID][k:k + 60])
                    while len(aastrline) == 60:
                        print(aastrline, end="\n", file=testrefaa);k += 60
                        aastrline = "".join(ref_amino_seq[tscptID][k:k + 60])
                    else:
                        print(aastrline, file=testrefaa)                                    
                                    
                    for snpInCDS in linetoCDSMap:
                        print(*(list(snpInCDS) +list(linetoCDSMap[snpInCDS][0])+ linetoCDSMap[snpInCDS][1:]), sep="\t", file=cdsVF)
    intergenicVF.close()
    cdsVF.close()
    intronVF.close()
    utrVF.close()
    genegrouptest.close()
    testrefaa.close()
    testmutcds.close()
    mutaa.close()
    print("finish")
