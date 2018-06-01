'''
Created on 2014-12-13

@author: liurui
'''
import copy
from optparse import OptionParser
import pickle
import re
from NGS.BasicUtil import VCFutil
from src.NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm

# testfile


parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa",
                  help="reference.fa")
parser.add_option("-g", "--gtffile", dest="gtffile", help="gtffile")
parser.add_option("-v", "--variantstable", dest="variantstable", help="variants")
parser.add_option("-V", "--vcffilename", dest="vcffilename", help="variants")
parser.add_option("-b", "--bedfiles", dest="bedfiles", action="append", default=[], help="bedfiles")
parser.add_option("-o", "--outputpath", dest="outputpath", help="default infile1_infile2")
parser.add_option("-c", "--chromtablename", dest="chromtablename")
parser.add_option("-m", "--minlength", dest="minlength")
parser.add_option("-5", "--TSSregion", dest="TSSregion", default="0", help="")
parser.add_option("-3", "--utr3_region", dest="utr3_region", default="0")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
refFastaName = options.reffa
reffa=options.reffa
reffahandler=open(reffa,"r")


minlength = options.minlength
chromtable = options.chromtablename#Util.pekingduckchromtable#
dbchromtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
if options.vcffilename!=None:
    variantstablename = re.search(r"[^/]*$",options.vcffilename).group(0)
    vcfobj=VCFutil.VCF_Data(options.vcffilename)
    titlelist=vcfobj.VcfIndexMap["title"]
    titlelist=titlelist[0:2]+titlelist[3:5]+titlelist[7:]
    pos_idx=0;ref_idx=1;alt_idx=2
else:
    variantstablename=options.variantstable.strip()
    dbvariantstools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.vcfdbname)
    pos_idx=1;ref_idx=3;alt_idx=4
    titlelist = [a[0].strip() for a in dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + variantstablename + "'")]

gtfMap,utrMap,allgeneSetMap = Util.getGtfMap(options.gtffile)
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

intergenicVF = open(outputpath + variantstablename+".intergenic", 'w')
cdsVF = open(outputpath + variantstablename+".cds", 'w')
intronVF = open(outputpath + variantstablename+".intron", 'w')
utrVF = open(outputpath + variantstablename+".utr", 'w')
print(*(titlelist + ["trscptID", "geneID", "strand", "cdsidx", "refcodon", "refaa", "altcodon", "altaa"]), sep="\t", file=cdsVF)
print(*(titlelist + ["trscptID", "geneID", "strand", "intronidx"]), sep="\t", file=intronVF)
print(*(titlelist + ["trscptID", "geneID", "strand", "5'/3'"]), sep="\t", file=utrVF)
print(*(titlelist+["tag"]), sep="\t", file=intergenicVF)

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
    try:
        refidxByChr = pickle.load(open(reffa+ ".myfasteridx", 'rb'))
#         refidxByChr = pickle.load(open(reffastaidxName, 'rb'))
    except IOError:
        Util.generateFasterRefIndex(reffa, reffa+ ".myfasteridx")
        refidxByChr = pickle.load(open(reffa+ ".myfasteridx", 'rb'))
#         Util.generateIndexByChrom(refFastaName, reffastaidxName)
#         refidxByChr = pickle.load(open(reffastaidxName, 'rb'))
        
    totalChroms = dbchromtools.operateDB("select", "select count(*) from " + chromtable + " where chrlength>=" + minlength)[0][0]
    for i in range(0, totalChroms, 20):
        currentsql = sql + " order by " + primaryID + " limit " + str(i) + ",20"
#         testcurrentsql=sql + " and chrID='KB742382.1' order by " + primaryID + " limit " + str(i) + ",20"
        result = dbchromtools.operateDB("select", currentsql)
        for row in result:
            seektuple = ()
            currentchrID = row[0]
            currentchrLen = int(row[1])
            if currentchrID not in gtfMap:
                if options.vcffilename !=None:
                    snps=vcfobj.getVcfListByChrom(currentchrID)
                elif options.variantstable!=None:
                    snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(0) + " and snp_pos<" + str(currentchrLen) + " order by snp_pos")
                for snp in snps:
                    print(currentchrID,*(snp[pos_idx:]+tuple(["nogene"])), sep="\t", file=intergenicVF)
                print("no gene in the chrom:", currentchrID)
                continue
            ####################### 
            
            
            ########################### collect snp in  upstream and downstream of the genegroup
            if gtfMap[currentchrID][0][1]=="+":
                firstutrstartpos=gtfMap[currentchrID][0][2]-TSSregionlen
            else:
                firstutrstartpos=gtfMap[currentchrID][0][2]-utr3_region
            if currentchrID in utrMap:
                if gtfMap[currentchrID][0][0] in utrMap[currentchrID]:
                    if utrMap[currentchrID][gtfMap[currentchrID][0][0]][0][2]<gtfMap[currentchrID][0][2]:
                        firstutrstartpos=utrMap[currentchrID][gtfMap[currentchrID][0][0]][0][1]
            if options.vcffilename !=None:
                snps=vcfobj.getVcfListByChrom(currentchrID,1,firstutrstartpos)
            elif options.variantstable!=None:
                
                snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(0) + " and snp_pos<" + str(firstutrstartpos) + " order by snp_pos")
            for snp in snps:
                print(currentchrID,*(snp[pos_idx:]+tuple(["begin"])), sep="\t", file=intergenicVF)
            GeneGrouplist = Util.getGeneGrouplist(gtfMap[currentchrID])
            

            
            for gene in GeneGrouplist[-1][1:]:
                if gene[3]==GeneGrouplist[-1][0]:
                    transcript_id=gene[0]
                    if gene[1]=="+":
                        lastutr=gene[3]+utr3_region
                    else:
                        lastutr=gene[3]+TSSregionlen
                    if currentchrID in utrMap:
                        if transcript_id in utrMap[currentchrID]:
                            if utrMap[currentchrID][transcript_id][-1][2]>gene[3]:
                                lastutr=utrMap[currentchrID][transcript_id][-1][2]                    
                    break
            if options.vcffilename!=None:
                
                snps=vcfobj.getVcfListByChrom(currentchrID, lastutr, int(currentchrLen))
            else:
                snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(lastutr) + " and snp_pos<" + str(currentchrLen) + " order by snp_pos")

            for snp in snps:
                print(currentchrID,*(snp[pos_idx:]+tuple(["end"])), sep="\t", file=intergenicVF)              
#  ###################################################           print("collect intergenic rear part of chr :",currentchrID,GeneGrouplist[-1][0],str(currentchrLen))
            lasttscptID_endpos=None
            for geneGroup in GeneGrouplist:
                grouplist_idx=GeneGrouplist.index(geneGroup)
#                 if grouplist_idx>=1:
#                     intergennicsnps=dbvariantstools.operateDB("select","select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(lasttscptID_endpos) + " and snp_pos<" + str(geneGroup[1][2]) + " order by snp_pos");
#                     for snp in intergennicsnps:
#                         print(*snp,sep='\t',file=intergenicVF)
                RefSeqMap = Util.getRefSeqBypos_faster(reffahandler, refidxByChr, currentchrID, geneGroup[1][2], geneGroup[0], currentchrLen, seektuple)
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
                    
                linetoCDSMap = {};linetoIntronMap = {}
                if grouplist_idx>=1 and grouplist_idx<len(GeneGrouplist):

                    if options.vcffilename!=None:
                        intergennicsnps=vcfobj.getVcfListByChrom(currentchrID,GeneGrouplist[grouplist_idx-1][0]+utr3_region,geneGroup[1][2]-TSSregionlen)                      
                    else:
                        intergennicsnps=dbvariantstools.operateDB("select","select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(GeneGrouplist[grouplist_idx-1][0]) + " and snp_pos<" + str(geneGroup[1][2]) + " order by snp_pos");
                    for snp in intergennicsnps:
                        print(currentchrID,*(snp[pos_idx:]+tuple(["before"+tscptID])),sep='\t',file=intergenicVF)     
                for gene_idx in range(1, len(geneGroup)):
                    tscptID = geneGroup[gene_idx][0]
                    currenttscptID_endpos=geneGroup[gene_idx][3]+utr3_region
                    currenttscptID_startpos=geneGroup[gene_idx][2]-TSSregionlen
                    if geneGroup[gene_idx][1] == '+':
                        
#                       collect  both  5UTR and 3UTR
                        if currentchrID in utrMap:
                            if tscptID in utrMap[currentchrID]:
                                i=0;exist_3utr=False;exist_5utr=False
                                for utr_type,utr_s,utr_e in utrMap[currentchrID][tscptID]:
                                    if i>=1:
                                        if options.vcffilename!=None:
                                            
                                            snps=vcfobj.getVcfListByChrom(currentchrID, utrMap[currentchrID][i-1][2], utr_s)
                                        else:
                                            snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(utrMap[currentchrID][i-1][2]) + " and snp_pos<" + str(utr_s) + " order by snp_pos")
#                                         for snp in snps:
#                                             print(currentchrID,*(snp[pos_idx:]+tuple(["betweenutr"+tscptID])), sep="\t", file=intergenicVF)
                                    elif i+1==len(utrMap[currentchrID][tscptID]) and False:
                                        print("unfinished")
                                        
                                    if utr_e>geneGroup[gene_idx][3]:
                                        utrMap[currentchrID][i]=["3UTR",utr_s,utr_e]
                                        currenttscptID_endpos=utr_e
                                        exist_3utr=True
                                    if utr_s<geneGroup[gene_idx][2]:
                                        utrMap[currentchrID][i]=["5UTR",utr_s,utr_e]
                                        if not exist_5utr:
                                            currenttscptID_startpos=utr_s
                                        exist_5utr=True
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID, utr_s,utr_e)
                                    else:
                                        snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(utr_s) + " and snp_pos<" + str(utr_e) + " order by snp_pos")
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","+",utrMap[currentchrID][i][0],sep="\t",file=utrVF)               
                                    i+=1
                                if not exist_3utr:
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID, geneGroup[gene_idx][3],geneGroup[gene_idx][3]+utr3_region)
                                    else:
                                        
                                        snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][3]) + " and snp_pos<" + str(geneGroup[gene_idx][3]+utr3_region) + " order by snp_pos")
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","+","3UTR",sep="\t",file=utrVF)
                                if not exist_5utr:
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID, geneGroup[gene_idx][2]-TSSregionlen,geneGroup[gene_idx][2])
                                    else:
                                        snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]-TSSregionlen) + " and snp_pos<" + str(geneGroup[gene_idx][2]) + " order by snp_pos")
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","+","5UTR",sep="\t",file=utrVF)
                        else:#no utr annotation in this gene
                            if options.vcffilename!=None:
                                
                                snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][3],geneGroup[gene_idx][3]+utr3_region)
                            else:
                                snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][3]) + " and snp_pos<" + str(geneGroup[gene_idx][3]+utr3_region) + " order by snp_pos")
                            for snp in snps:
                                print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                print(tscptID,"geneID","+","3UTR",sep="\t",file=utrVF)
                            if options.vcffilename!=None:
                                
                                snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][2]-TSSregionlen,geneGroup[gene_idx][2])
                            else:
                                snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]-TSSregionlen) + " and snp_pos<" + str(geneGroup[gene_idx][2]) + " order by snp_pos")
                            
                            for snp in snps:
                                print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                print(tscptID,"geneID","+","5UTR",sep="\t",file=utrVF)
            ######################################## UTR collection finish
                        if options.vcffilename!=None:
                            
                            snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][2],geneGroup[gene_idx][3])
                        else:
                            
                            snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]) + " and snp_pos<" + str(geneGroup[gene_idx][3]) + " order by snp_pos")  
                        for snp in snps:
                            snppos = snp[pos_idx]
                            refbase = snp[ref_idx]
                            altbase = snp[alt_idx]
                            cdsidx = 3
                            Intron_idx = -1
                            if re.search(r'[^a-zA-Z]', altbase) != None or len(altbase)>1 or len(refbase)>1:  # contain ',' ie. multiple alle
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
                                            if tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1])) in linetoCDSMap:
                                                linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1] + ";" + tscptID;linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2] + ";geneID";linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3] + ";+";linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4] + ";" + str(cdsidx-3)
                                            else:
                                                linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))] = [snp[alt_idx+1:], tscptID, "geneID", "+", str(cdsidx-3)]
###########################################################################
                                elif snppos > elemEnd and snppos < geneGroup[gene_idx][cdsidx + 1][1]:
                                    Intron_idx = cdsidx - 3
                            if Intron_idx != -1:
                                if tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1])) in linetoIntronMap:
                                    linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1]+";"+tscptID;linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2]+";geneID";linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3]+";+";linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4]+";"+str(Intron_idx)
                                else:
                                    linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))]=[snp[alt_idx+1:],tscptID, "geneID", "+", str(Intron_idx)]
#                                 print(*(list(snp) + [tscptID, "geneID", "+", Intron_idx]), sep="\t", file=intronVF)
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
                                    if (currentchrID,snppos_cds,ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)])==9:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] + ";" + CodonTable[codon];linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] + ";" + CodonTable[codon_m]
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)] += [codon, CodonTable[codon], codon_m, CodonTable[codon_m]]
                                    else:
                                        print(currentchrID,snppos_cds,ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                except KeyError:
#                                     snppos_cds, ref_base_cds, alt_base_cds = Util.getSNPrecInCDS(i, len(tscptSeqAllCds[tscptID]), codon, codon_m, cds_frame[tscptID], gene)
                                    if (currentchrID,snppos_cds,ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)])==9:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] + ";X";linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] + ";X"
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)] += [codon, "X", codon_m, "X"]
                                    else:
                                        print(currentchrID,snppos_cds,ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                                                        
                            try:
                                ref_amino_seq[tscptID].append(CodonTable["".join(tscptSeqAllCds[tscptID][i:i + 3]).lower()])
                            except KeyError:
                                ref_amino_seq[tscptID].append("X")
                                    
                    else:  # strand == '-'
                        ###############################
                        if currentchrID in utrMap:
                            if tscptID in utrMap[currentchrID]:
                                i=0;exist_3utr=False;exist_5utr=False
                                for utr_type,utr_s,utr_e in utrMap[currentchrID][tscptID]:
                                    if i>=1:
                                        if options.vcffilename!=None:
                                            
                                            snps=vcfobj.getVcfListByChrom(currentchrID,utrMap[currentchrID][i-1][2],utr_s)
                                        else:
                                            snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(utrMap[currentchrID][i-1][2]) + " and snp_pos<" + str(utr_s) + " order by snp_pos")
#                                         for snp in snps:
#                                             print(currentchrID,*(snp[pos_idx:]+tuple(["betweenutr"+tscptID])), sep="\t", file=intergenicVF)
                                    elif i+1==len(utrMap[currentchrID][tscptID]) and False:
                                        if options.vcffilename!=None:
                                            
                                            snps=vcfobj.getVcfListByChrom(currentchrID,utrMap[currentchrID][i-1][2],geneGroup[gene_idx][2])
                                        else:
                                            snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(utrMap[currentchrID][i-1][2]) + " and snp_pos<" + str(geneGroup[gene_idx][2]) + " order by snp_pos")
                                        print("unfinished")
                                        

                                    if utr_e>geneGroup[gene_idx][3]:
                                        utrMap[currentchrID][i]=["5UTR",utr_s,utr_e]
                                        currenttscptID_endpos=utr_e
                                        exist_5utr=True
                                    if utr_s<geneGroup[gene_idx][2]:
                                        utrMap[currentchrID][i]=["3UTR",utr_s,utr_e]
                                        if not exist_3utr:
                                            currenttscptID_startpos=utr_s
                                        exist_3utr=True
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID,utr_s,utr_e)
                                    else:
                                        snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(utr_s) + " and snp_pos<" + str(utr_e) + " order by snp_pos")
                                    
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","-",utrMap[currentchrID][i][0],sep="\t",file=utrVF)               
                                    i+=1
                                if not exist_3utr:
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][2]-utr3_region,geneGroup[gene_idx][2])
                                    else:
                                        snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]-utr3_region) + " and snp_pos<" + str(geneGroup[gene_idx][2]) + " order by snp_pos")
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","-","3UTR",sep="\t",file=utrVF)
                                if not exist_5utr:
                                    if options.vcffilename!=None:
                                        
                                        snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][3],geneGroup[gene_idx][3]+TSSregionlen)
                                    else:
                                        snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][3]) + " and snp_pos<" + str(geneGroup[gene_idx][3]+TSSregionlen) + " order by snp_pos")
                                    for snp in snps:
                                        print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                        print(tscptID,"geneID","-","5UTR",sep="\t",file=utrVF)
                        else:
                            if options.vcffilename!=None:
                                
                                snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][3]-utr3_region,geneGroup[gene_idx][3])
                            else:
                                snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][3]-utr3_region) + " and snp_pos<" + str(geneGroup[gene_idx][3]) + " order by snp_pos")
                            for snp in snps:
                                print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                print(tscptID,"geneID","-","3UTR",sep="\t",file=utrVF)
                            if options.vcffilename!=None:
                                snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][2],geneGroup[gene_idx][2]+TSSregionlen)
                            else:
                                snps=dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]) + " and snp_pos<" + str(geneGroup[gene_idx][2]+TSSregionlen) + " order by snp_pos")
                            for snp in snps:
                                print(currentchrID,*snp[pos_idx:],sep="\t",end="\t",file=utrVF)
                                print(tscptID,"geneID","-","5UTR",sep="\t",file=utrVF)
                        ################################# utr region collection finished #############################
                        if options.vcffilename!=None:
                            snps=vcfobj.getVcfListByChrom(currentchrID,geneGroup[gene_idx][2],geneGroup[gene_idx][3])
                        else:
                            snps = dbvariantstools.operateDB("select", "select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>" + str(geneGroup[gene_idx][2]) + " and snp_pos<" + str(geneGroup[gene_idx][3]) + " order by snp_pos")
                        for snp in snps:
                            snppos = snp[pos_idx]
                            refbase = snp[ref_idx]
                            altbase = snp[alt_idx]
                            cdsidx = 3
                            Intron_idx = -1
                            if re.search(r'[^a-zA-Z]', altbase) != None or len(altbase)>1 or len(refbase)>1:  # contain ',' ie. multiple alle
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
                                            if tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1])) in linetoCDSMap:
                                                linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1] + ";" + tscptID;linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2] + ";geneID";linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3] + ";-";linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4] = linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4] + ";" + str(len(cds_frame[tscptID]) - (cdsidx-4))
                                            else:
                                                linetoCDSMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))] = [snp[alt_idx+1:],tscptID, "geneID", "-", str(len(cds_frame[tscptID]) - (cdsidx-4))]   
########################################################
                                elif snppos > elemEnd and snppos < geneGroup[gene_idx][cdsidx + 1][1]:
                                    Intron_idx = len(cds_frame[tscptID]) - (cdsidx-3)
                            if Intron_idx != -1:
                                if tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1])) in linetoIntronMap:
                                    linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][1]+";"+tscptID;linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][2]+";geneID";linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][3]+";-";linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4]=linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))][4]+";"+str(Intron_idx)
                                else:
                                    linetoIntronMap[tuple([currentchrID]+list(snp[pos_idx:pos_idx+1]+snp[ref_idx:alt_idx+1]))]=[snp[alt_idx+1:],tscptID, "geneID", "-", str(Intron_idx)]
#                                 print(*(list(snp) + [tscptID, "geneID", "-", Intron_idx]), sep="\t", file=intronVF)
# translate protein
                        gene=geneGroup[gene_idx]
                        mutat_amino_seq[tscptID] = []
                        ref_amino_seq[tscptID] = []
                        tscptSeqAllCds_mut[tscptID].reverse()
#                             tscptSeqAllCds_mut[tscptID] = list(tscptSeqAllCds_mut_str)
#############################      produce linetoCDSMap  ############
                        tscptSeqAllCds_Revr_Cmplm = Util.complementary(tscptSeqAllCds[tscptID])
                        tscptSeqAllCds_Revr_Cmplm.reverse()
                        tscptSeqAllCds_mut_Revr_Cmplm = Util.complementary(tscptSeqAllCds_mut[tscptID])
                        tscptSeqAllCds_mut[tscptID]=tscptSeqAllCds_mut_Revr_Cmplm
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
#                                     print((currentchrID,snppos_cds,ref_base_cds, alt_base_cds))
                                    if (currentchrID,snppos_cds,ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)])==9:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] + ";" + CodonTable[codon];linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] + ";" + CodonTable[codon_m]
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)] += [codon, CodonTable[codon], codon_m, CodonTable[codon_m]]
                                    else:
                                        print(currentchrID,snppos_cds,ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                                except KeyError:
                                    print("except KEYERROR")
                                    if (currentchrID,snppos_cds,ref_base_cds, alt_base_cds) in linetoCDSMap:
                                        if len(linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)])==9:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][5] + ";" + codon;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][6] + ";X";linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][7] + ";" + codon_m;linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] = linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)][8] + ";X"
                                        else:
                                            linetoCDSMap[(currentchrID,snppos_cds,ref_base_cds, alt_base_cds)] += [codon, "X", codon_m, "X"]
                                    else:
                                        print(currentchrID,snppos_cds,ref_base_cds, alt_base_cds,i,"should in the linetoCDSMap:",linetoCDSMap,tscptID,codon, codon_m,cds_frame[tscptID],"\n",linetoCDSMap,file=open("wrong.txt",'a'))
                            try:
                                mutat_amino_seq[tscptID].append(CodonTable[codon_m])
                                ref_amino_seq[tscptID].append(CodonTable["".join(tscptSeqAllCds_Revr_Cmplm[i:i + 3]).lower()])
                            except KeyError:
                                ref_amino_seq[tscptID].append("X")
                            
                    print(">"+tscptID+" transcript:" + tscptID, file=testmutcds)
                    print(">" + tscptID, file=testrefaa)
                    print(">transcript:" + tscptID, file=mutaa)
                    k = 0
                    cdsstrline = "".join(tscptSeqAllCds_mut[tscptID][k:k + 60]).upper()
                    while len(cdsstrline) == 60:
                        print(cdsstrline, file=testmutcds);k += 60
                        cdsstrline = "".join(tscptSeqAllCds_mut[tscptID][k:k + 60]).upper()
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
                    if gene_idx>1:
    #                     print("intergenicVF","select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(lasttscptID_endpos) + " and snp_pos<" + str(currenttscptID_startpos) + " order by snp_pos")
                        if options.vcffilename!=None:
                            intergennicsnps=vcfobj.getVcfListByChrom(currentchrID,lasttscptID_endpos,currenttscptID_startpos)
                        else:
                            intergennicsnps=dbvariantstools.operateDB("select","select * from " + variantstablename + " where chrID='" + currentchrID + "' and snp_pos>=" + str(lasttscptID_endpos) + " and snp_pos<" + str(currenttscptID_startpos) + " order by snp_pos");
                        for snp in intergennicsnps:
                            print(currentchrID,*(snp[pos_idx:]+tuple(["before"+tscptID])),sep='\t',file=intergenicVF)                
                    lasttscptID_endpos=currenttscptID_endpos
                for snpInIntron in sorted(linetoIntronMap.keys()):
                    print(*(list(snpInIntron)+list(linetoIntronMap[snpInIntron][0])+linetoIntronMap[snpInIntron][1:]),sep="\t",file=intronVF)           
                for snpInCDS in sorted(linetoCDSMap.keys()):
                    print(*(list(snpInCDS) +list(linetoCDSMap[snpInCDS][0])+ linetoCDSMap[snpInCDS][1:]), sep="\t", file=cdsVF)
                
    intergenicVF.close()
    cdsVF.close()
    intronVF.close()
    utrVF.close()
    testrefaa.close()
    testmutcds.close()
    mutaa.close()
    print("finish")
