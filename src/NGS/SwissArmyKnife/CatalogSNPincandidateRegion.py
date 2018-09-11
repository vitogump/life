'''
Created on 2015-11-13

@author: liurui
'''
import copy
from optparse import OptionParser
import os,config
import pickle
import re

from Bio import SeqIO, pairwise2

from src.NGS.BasicUtil import geneUtil, Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa", help="reference.fa")
parser.add_option("-t", "--topleveltable", dest="topleveltable", help="reference.fa")
parser.add_option("-g", "--gtffile", dest="gtffile",nargs=3, help="gtffile upextend downextend")
parser.add_option("-c", "--catalogVariantfold", dest="catalogVariantfold", help="catalogVariantfold")
parser.add_option("-d", "--domesticsnptable", dest="domesticsnptable",action="append", help="variants")
parser.add_option("-w", "--wildsnptable", dest="wildsnptable",action="append", help="variants")
parser.add_option("-R", "--candidateRegion", dest="candidateRegion", action="append", default=[], help="bedfiles")
parser.add_option("-o", "--outputpath", dest="outputpath", help="default infile1_infile2")
parser.add_option("-f","--filterfreq",dest="filterfreq",default="0",help="filterfreq ")
parser.add_option("-E", "--elementfold", dest="elementfold",action="append",nargs=2,help="fold targetseqnamesubstr")


(options, args) = parser.parse_args()
dbvariantstools = dbm.DBTools(config.ip, config.username, config.password, config.vcfdbname)
pathtoblastn="/pub/tool/blast_set/ncbi-blast-2.2.29+/bin/blastn "
reffa=options.reffa
reffahander=open(reffa,"r")
dbchromtools = dbm.DBTools(config.ip, config.username, config.password, config.genomeinfodbname)
try:
    refidxByChr = pickle.load(open(reffa+ ".myfasteridx", 'rb'))
except IOError:
    Util.generateFasterRefIndex(reffa, reffa+ ".myfasteridx")
    refidxByChr = pickle.load(open(reffa+ ".myfasteridx", 'rb'))
if __name__ == '__main__':
    #prepare elementfoldbed,gtfMap,utrMap
    chromSet=set()
    for bedlikeselectedgenefile in options.candidateRegion:
        f=open(bedlikeselectedgenefile,'r');f.readline()#title
        for regionline in f:
            linelist=re.split(r"\s+",regionline.strip())
            chromSet.add(linelist[0].strip())
        f.close()
        
        
        
    elementfoldbed={}
    for elementfold,targetseqnamesubstr in options.elementfold[:]:
        if elementfold.endswith("/") or elementfold.endswith("\\"):
            elementfold=elementfold[:-1]
        allseqtobed=geneUtil.make_getElemBed(elementfold, targetseqnamesubstr, pathtoblastn, reffa)
        elementfoldbed=dict(elementfoldbed,**allseqtobed)#merge into elementfoldbed
    gtfMap,utrMap,allgeneSetMap = geneUtil.getGtfMap(options.gtffile[0])
    upextend=int(options.gtffile[1]);downextend=int(options.gtffile[2])
    if options.catalogVariantfold.endswith("/") or options.catalogVariantfold.endswith("\\"):
        catalogVariantfold=options.catalogVariantfold[:-1]
    else:
        catalogVariantfold=options.catalogVariantfold
    w_cds_map_list=[];d_cds_map_list=[]
    w_intron_map_list=[];d_intron_map_list=[]
    w_utr_map_list=[];d_utr_map_list=[]
    for wild_table in options.wildsnptable[:]:
        w_cdf_file=open(catalogVariantfold+"/"+wild_table+".cds","r");w_cds_map={};titlelist=re.split(r"\s+",w_cdf_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for cdsline in w_cdf_file:
            cdslinelist=re.split(r"\s+",cdsline.strip())
            if cdslinelist[0].strip() not in chromSet:
                continue
            if cdslinelist[0] in w_cds_map:
                w_cds_map[cdslinelist[0]].append((int(cdslinelist[1]),cdslinelist[3],cdslinelist[4],";".join(cdslinelist[tp_idx:]),float(cdslinelist[af_idx]),cdslinelist[dp_an_idx]))
            else:
                w_cds_map[cdslinelist[0]]=[(int(cdslinelist[1]),cdslinelist[3],cdslinelist[4],";".join(cdslinelist[tp_idx:]),float(cdslinelist[af_idx]),cdslinelist[dp_an_idx])]
        w_cdf_file.close();w_cds_map_list.append(copy.deepcopy(w_cds_map))
         
        w_intron_file=open(catalogVariantfold+"/"+wild_table+".intron","r");w_intron_map={} ;titlelist=re.split(r"\s+",w_intron_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for intronline in w_intron_file:
            intronlinelist=re.split(r"\s+",intronline.strip())
            if intronlinelist[0].strip() not in chromSet:
                continue
            if intronlinelist[0] in w_intron_map:
                w_intron_map[intronlinelist[0]].append((int(intronlinelist[1]),intronlinelist[3],intronlinelist[4],";".join(intronlinelist[tp_idx:]),float(intronlinelist[af_idx]),intronlinelist[dp_an_idx]))
            else:
                w_intron_map[intronlinelist[0]]=[(int(intronlinelist[1]),intronlinelist[3],intronlinelist[4],";".join(intronlinelist[tp_idx:]),float(intronlinelist[af_idx]),intronlinelist[dp_an_idx])]
        w_intron_file.close();w_intron_map_list.append(copy.deepcopy(w_intron_map))
        w_utr_file=open(catalogVariantfold+"/"+wild_table+".utr","r");w_utr_map={};titlelist=re.split(r"\s+",w_utr_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for utrline in w_utr_file:
            intergeneticlinelist=re.split(r"\s+",utrline)
            if intergeneticlinelist[0].strip() not in chromSet:
                continue
            
            if intergeneticlinelist[0] in w_utr_map:
                w_utr_map[intergeneticlinelist[0]].append((int(intergeneticlinelist[1]),intergeneticlinelist[3],intergeneticlinelist[4],";".join(intergeneticlinelist[tp_idx:]),float(intergeneticlinelist[af_idx]),intergeneticlinelist[dp_an_idx]))
            else:
                w_utr_map[intergeneticlinelist[0]]=[(int(intergeneticlinelist[1]),intergeneticlinelist[3],intergeneticlinelist[4],";".join(intergeneticlinelist[tp_idx:]),float(intergeneticlinelist[af_idx]),intergeneticlinelist[dp_an_idx])]
        w_utr_file.close();w_utr_map_list.append(copy.deepcopy(w_utr_map))
    for dom_table in options.domesticsnptable[:]:
        d_cdf_file=open(catalogVariantfold+"/"+dom_table+".cds","r");d_cds_map={};titlelist=re.split(r"\s+",d_cdf_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for cdsline in d_cdf_file:
            cdslinelist=re.split(r"\s+",cdsline)
            if cdslinelist[0].strip() not in chromSet:
                continue
            if cdslinelist[0] in d_cds_map:
                d_cds_map[cdslinelist[0]].append((int(cdslinelist[1]),cdslinelist[3],cdslinelist[4],";".join(cdslinelist[tp_idx:]),float(cdslinelist[af_idx]),cdslinelist[dp_an_idx]))
            else:
                d_cds_map[cdslinelist[0]]=[(int(cdslinelist[1]),cdslinelist[3],cdslinelist[4],";".join(cdslinelist[tp_idx:]),float(cdslinelist[af_idx]),cdslinelist[dp_an_idx])]
        d_cdf_file.close();d_cds_map_list.append(copy.deepcopy(d_cds_map))
        d_intron_file=open(catalogVariantfold+"/"+dom_table+".intron","r");d_intron_map={};titlelist=re.split(r"\s+",d_intron_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for intronline in d_intron_file:
            intronlinelist=re.split(r"\s+",intronline.strip())
            if intronlinelist[0].strip() not in chromSet:
                continue
            if intronlinelist[0] in d_intron_map:
                d_intron_map[intronlinelist[0]].append((int(intronlinelist[1]),intronlinelist[3],intronlinelist[4],";".join(intronlinelist[tp_idx:]),float(intronlinelist[af_idx]),intronlinelist[dp_an_idx]))
            else:
                d_intron_map[intronlinelist[0]]=[(int(intronlinelist[1]),intronlinelist[3],intronlinelist[4],";".join(intronlinelist[tp_idx:]),float(intronlinelist[af_idx]),intronlinelist[dp_an_idx])]
        d_intron_file.close();d_intron_map_list.append(copy.deepcopy(d_intron_map))
        d_utr_file=open(catalogVariantfold+"/"+dom_table+".utr","r");d_utr_map={};titlelist=re.split(r"\s+",d_utr_file.readline().strip());af_idx=titlelist.index("AF");tp_idx=titlelist.index("trscptID")
        if "DP" in titlelist:
            dp_an_idx=titlelist.index("DP")
        elif "AN" in titlelist:
            dp_an_idx=titlelist.index("AN")
        for utrline in d_utr_file:
            intergeneticlinelist=re.split(r"\s+",utrline)            
            if intergeneticlinelist[0].strip() not in chromSet:
                continue

            if intergeneticlinelist[0] in d_utr_map:
                d_utr_map[intergeneticlinelist[0]].append((int(intergeneticlinelist[1]),intergeneticlinelist[3],intergeneticlinelist[4],";".join(intergeneticlinelist[tp_idx:]),float(intergeneticlinelist[af_idx]),intergeneticlinelist[dp_an_idx]))
            else:
                d_utr_map[intergeneticlinelist[0]]=[(int(intergeneticlinelist[1]),intergeneticlinelist[3],intergeneticlinelist[4],";".join(intergeneticlinelist[tp_idx:]),float(intergeneticlinelist[af_idx]),intergeneticlinelist[dp_an_idx])]
        d_utr_file.close();d_utr_map_list.append(copy.deepcopy(d_utr_map))
    ###########order##################
#     for utr_map in d_utr_map_list:
        
    ###sql statement
    sqlselectstatementpart_count_left="select chrID,snp_pos,"+options.wildsnptable[0] +".ref_base,"+options.wildsnptable[0] +".alt_base"# from "+options.wildsnptable[0];
    sqlselectstatementpart_count_right="select chrID,snp_pos,"+options.wildsnptable[0] +".ref_base,"+options.wildsnptable[0] +".alt_base"# from "+options.wildsnptable[0]
    for wild_domtable in options.wildsnptable[1:]+options.domesticsnptable[:]:
        sqlselectstatementpart_count_left+=(","+wild_domtable +".ref_base,"+wild_domtable +".alt_base")
        sqlselectstatementpart_count_right+=(","+wild_domtable +".ref_base,"+wild_domtable +".alt_base")
    sqlselectstatementpart_count_left+=(","+options.wildsnptable[0] +".AF");sqlselectstatementpart_count_right+=(","+options.wildsnptable[0] +".AF")
    for wild_domtable in options.wildsnptable[1:]+options.domesticsnptable[:]:
        sqlselectstatementpart_count_left+=(","+wild_domtable+ ".AF")
        sqlselectstatementpart_count_right+=(","+wild_domtable+".AF")
    sqlselectstatementpart_count_left+=" from "+options.topleveltable+" left join "+options.wildsnptable[0]+" using(chrID,snp_pos) ";sqlselectstatementpart_count_right+=" from "+options.topleveltable+" left join "+options.wildsnptable[0]+" using(chrID,snp_pos) "
    for wild_domtable in options.wildsnptable[1:]+options.domesticsnptable[:]:
        sqlselectstatementpart_count_left+=(" left join "+wild_domtable+" using(chrID,snp_pos) ")
        sqlselectstatementpart_count_right+=(" right join "+wild_domtable+" using(chrID,snp_pos) ")

    
    ###process candidate region with gene
#     outputpath=options.outputpath
#     if outputpath.endswith("/") or outputpath.endswith("\\"):
#         outputpath=outputpath[:-1]
    outfile=open(options.outputpath+".conserved","w")
    utroutfile=open(options.outputpath+".utr","w")
    print("delta_AF\tw_af\td_af\tsnp_pos\tref_base\talt_base\tpID",file=utroutfile)
    cdsoutfile=open(options.outputpath+".cds","w")
    print("delta_AF\tw_af\td_af\tsnp_pos\tref_base\talt_base",file=cdsoutfile)
    intronoutfile=open(options.outputpath+".intron","w")
    print("delta_AF\tw_af\td_af\tsnp_pos\tref_base\talt_base",file=intronoutfile)
    intergenic=open(options.outputpath+".intergenic","w")
    print("delta_AF\tw_af\td_af\tsnp_pos\tref_base\talt_base",file=intergenic)
#     testf=open("f.txt",'w')
#     testff=open("ff.txt",'w')
    print("delta_AF\tsnppos\tref_base\talt_base\tw_af\tsnp_idx_in_seq\t",file=outfile)
#     for speciese in sorted(muscleout_seqmap.keys()):
#         print(speciese,end="\t",file=outfile) 
    print(utrMap)
    for bedlikeselectedgenefile in options.candidateRegion:
        candidateRegionbed={}
        f=open(bedlikeselectedgenefile,'r')
        f.readline()#title
        for regionline in f:
            linelist=re.split(r"\s+",regionline.strip())
            chrom=linelist[0].strip();regionstart=int(linelist[1]);regionend=int(linelist[2])
#for every region,that is every line   [(pos,deltaAF,REF,ALT),{}]
                    #overlapwithgene
            if len(linelist)>=9 and (linelist[8].find("8")==linelist[8].find("7")):#only return -1 the equal are
                trscpts=re.split(r",",linelist[7])
#                 print("overlapwithgene",regionline,file=testf)
                print("region:",regionline.strip(),file=utroutfile)
                print("region:",regionline.strip(),file=cdsoutfile)
                print("region:",regionline.strip(),file=intronoutfile)
                
                for tpID in trscpts:
                    allutrsnps_fromfile=[];allutrsnps_fromsql=[];allcdssnps=[];allintronsnps=[]
                    tpIDidx=-1
                    IngtfMap=False#tpID not in grfMap,that means no cds
                    if chrom in gtfMap :
                        for  tplist  in gtfMap[chrom]:
                            tpIDidx+=1
                            if tplist[0]==tpID:
                                vcflist_A_chrom_container=[]

                                overlap_start=max(tplist[2],regionstart)
                                overlap_end=min(tplist[3],regionend)
                                for snpmap in w_cds_map_list+d_cds_map_list:
                                    tempmap={}
                                    if chrom not in snpmap:
                                        tempmap[chrom]=[]
                                    else:
                                        tempmap[chrom]=snpmap[chrom]
                                    vcflist_A_chrom_container.append(copy.deepcopy(tempmap))

                                MultipleVcfMap_cds=Util.alinmultPopSnpPos(vcflist_A_chrom_container,"o")
                                a=len(allcdssnps)
                                allcdssnps+=geneUtil.collectSNP_locatInRegion(MultipleVcfMap_cds,chrom,overlap_start,overlap_end)
                                b=len(allcdssnps)
                                if a==b:
                                    print(regionline,"no cds snp" ,tplist[2],overlap_start,overlap_end)
                                vcflist_A_chrom_container_intron=[]
                                for snpmap in w_intron_map_list+d_intron_map_list:
                                    tempmap={}
                                    if chrom not in snpmap:
                                        tempmap[chrom]=[]
                                    else:
                                        tempmap[chrom]=snpmap[chrom]
                                    vcflist_A_chrom_container_intron.append(copy.deepcopy(tempmap))
                                 
                                MultipleVcfMap_intron=Util.alinmultPopSnpPos(vcflist_A_chrom_container_intron,"o")
 
                                allintronsnps+=geneUtil.collectSNP_locatInRegion(MultipleVcfMap_intron,chrom,overlap_start,overlap_end)
                                IngtfMap=True
                                break
                        else:
                            print("search all gtf,no found",tpID)
                            IngtfMap=False
                         
                    elif not IngtfMap and (chrom in allgeneSetMap) and (tpID in allgeneSetMap[chrom]):
                        print("exon,unknown cds gene,collect snp from mysql ")
                    
                    
                    
                    ####################
                    sql_file=""
#                     if chrom in utrMap and (tpID  in utrMap[chrom]):
#                         vcflist_A_chrom_container=[]
# #                         if chrom not in (d_utr_map_list+w_utr_map_list):
# #                             print(chrom,tpID,"give up")
# #                             continue
#                         for snpmap in w_utr_map_list+d_utr_map_list:
#                             tempmap={}
#                             if chrom not in snpmap:
#                                 tempmap[chrom]=[]
#                             else:
#                                 tempmap[chrom]=snpmap[chrom]
#                             vcflist_A_chrom_container.append(copy.deepcopy(tempmap))
#                         MultipleVcfMap=Util.alinmultPopSnpPos(vcflist_A_chrom_container, "o")
#                         if utrMap[chrom][tpID][0][1]<(allgeneSetMap[chrom][tpID][1]+allgeneSetMap[chrom][tpID][2])/2:
#                             sql_file="u"
#                         if utrMap[chrom][tpID][-1][2]>(allgeneSetMap[chrom][tpID][1]+allgeneSetMap[chrom][tpID][2])/2:
#                             sql_file+="d"
#                         for a,utrstart,utrend in utrMap[chrom][tpID]:
#                             if (utrend>regionstart and utrend<regionend) or (utrstart>regionstart and utrstart<regionend):
#                                 overlap_start=max(regionstart,utrstart)
#                                 overlap_end=min(regionend,utrend)
#                                 print("chrom",chrom,"overlap_start",overlap_start,"overlap_end", overlap_end)
#                                 allutrsnps_fromfile+=geneUtil.collectSNP_locatInRegion(MultipleVcfMap, chrom, overlap_start, overlap_end)
#                         sql_file="file"
                    if chrom in gtfMap and (tpID == gtfMap[chrom][tpIDidx][0]):
                        vcflist_A_chrom_container=[]
#                         if chrom not in (w_utr_map_list+d_utr_map_list):
#                             print(chrom,tpID,"give up")
#                             continue

                        if sql_file.find("u")==-1:
                            overlap_start=max(regionstart,(gtfMap[chrom][tpIDidx][2]-upextend))
                            overlap_end=min(regionend,gtfMap[chrom][tpIDidx][2])
                            print(tpID,sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
                            allutrsnps_fromsql+=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
#                             if len(temp)>1:
#                                 sql_file+="u"
#                                 allutrsnps_fromfile+=temp
                        if sql_file.find("d")==-1:
                            overlap_start=max(regionstart,(gtfMap[chrom][tpIDidx][3]))
                            overlap_end=min(regionend,(gtfMap[chrom][tpIDidx][3]+downextend))
                            print(tpID,sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
                            allutrsnps_fromsql+=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
#                             temp=geneUtil.collectSNP_locatInRegion(MultipleVcfMap, chrom, overlap_start, overlap_end)
#                             if len(temp)>1:
#                                 allutrsnps_fromfile+=temp
#                                 sql_file+="d" 
                        sql_file+="sql"   
                    elif chrom in allgeneSetMap and (tpID in allgeneSetMap[chrom]):

                        if sql_file.find("u")==-1 : 
                            overlap_start=max(regionstart,(allgeneSetMap[chrom][tpID][1]-upextend))
                            overlap_end=min(regionend,(allgeneSetMap[chrom][tpID][1]))                            
                            print(sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)+" union "+sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))#+" union "+sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)
                            allutrsnps_fromsql+=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
 
                        if sql_file.find("d")==-1:
                            overlap_start=max(regionstart,(allgeneSetMap[chrom][tpID][2]))
                            overlap_end=min(regionend,(allgeneSetMap[chrom][tpID][2]+downextend))                            
                            print(sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)+" union "+sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))#+" union "+sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)
                            allutrsnps_fromsql+=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
                        sql_file+="sql"
                    else:
                        pass
#                         sql_file="notfound"

#                     allcdssnps=list(set(allcdssnps));allintronsnps=list(set(allintronsnps));allutrsnps=list(set(allutrsnps))
                    snplist_posdeltaAFrefalt=[]
                    for snp in allcdssnps:
                        w_af=0;d_af=0;delta_af=0
                        for w_idx in range(len(options.wildsnptable)):
                            if snp[3+w_idx]==None:
                                w_af+=0
                            else:
 
                                w_af+=snp[3+w_idx][1]
                        for d_idx in range(len(options.domesticsnptable)):
                            if snp[3+len(options.wildsnptable)+d_idx]==None:
                                d_af+=0
                            else:
                                d_af+=snp[3+len(options.wildsnptable)+d_idx][1]
                        w_af=w_af/len(options.wildsnptable)
                        d_af=d_af/len(options.domesticsnptable)
                        delta_af=abs(w_af-d_af)
#                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                        if delta_af>=float(options.filterfreq):                         
                            snplist_posdeltaAFrefalt.append((delta_af,w_af,d_af,copy.deepcopy(snp)))
                    snplist_posdeltaAFrefalt.sort(key=lambda  listRec:listRec[0],reverse=True)
                    for delta_af,w_af,d_af,snp in snplist_posdeltaAFrefalt:
                        print('%.5f'%delta_af,w_af,d_af,*snp,sep="\t",file=cdsoutfile)
                    snplist_posdeltaAFrefalt=[]
                    if sql_file.find("u")!=-1 or sql_file.find("d")!=-1:
                        for snp in allutrsnps_fromfile:
                            print(snp)
                            w_af=0;d_af=0;delta_af=0
                            for w_idx in range(len(options.wildsnptable)):
                                if snp[3+w_idx]==None:
                                    w_af+=0
                                else:
    
                                    w_af+=snp[3+w_idx][1]
                            for d_idx in range(len(options.domesticsnptable)):
                                if snp[3+len(options.wildsnptable)+d_idx]==None:
                                    d_af+=0
                                else:
                                    d_af+=snp[3+len(options.wildsnptable)+d_idx][1]
                            w_af=w_af/len(options.wildsnptable)
                            d_af=d_af/len(options.domesticsnptable)
                            delta_af=abs(w_af-d_af)
    #                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                            if delta_af>float(options.filterfreq):                            
                                snplist_posdeltaAFrefalt.append((delta_af,w_af,d_af,copy.deepcopy(snp)))
                        snplist_posdeltaAFrefalt.sort(key=lambda  listRec:listRec[0],reverse=True)
                        for delta_af,w_af,d_af,snp in snplist_posdeltaAFrefalt:
                            print('%.5f'%delta_af,w_af,d_af,*snp,sep="\t",file=utroutfile)
                    if sql_file.find("sql")!=-1:
                        for snp in allutrsnps_fromsql:
                            w_af=0;d_af=0;delta_af=0
                            for w_idx in range(len(options.wildsnptable)):
                                if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx]==None:
                                    w_af+=0
                                else:
                                    ref_base=snp[2+w_idx*2]
                                    w_base=snp[3+w_idx*2]
                                    w_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx])
                            for d_idx in range(len(options.domesticsnptable)):
                                if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx]==None:
                                    d_af+=0
                                else:
                                    ref_base=snp[2+d_idx*2+len(options.wildsnptable)*2]
                                    w_base=snp[3+d_idx*2+len(options.wildsnptable)*2]
                                    d_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx])
                            w_af=w_af/len(options.wildsnptable)
                            d_af=d_af/len(options.domesticsnptable)
                            delta_af=abs(w_af-d_af)
#                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                            if delta_af>float(options.filterfreq):
                                snplist_posdeltaAFrefalt.append((int(snp[1]),delta_af,ref_base,w_base,w_af,d_af,tpID))
                        snplist_posdeltaAFrefalt.sort(key=lambda  listRec:listRec[1],reverse=True)
                        for snp_pos,delta_af,ref_base,w_base,w_af,d_af,tpID in snplist_posdeltaAFrefalt:
                            print('%.5f'%delta_af,w_af,d_af,snp_pos,ref_base,w_base,tpID,sep="\t",file=utroutfile)
                    else:
                        pass
#                         print(sql_file,file=utroutfile)

                
                    snplist_posdeltaAFrefalt=[]
                    for snp in allintronsnps:
                        w_af=0;d_af=0;delta_af=0
                        for w_idx in range(len(options.wildsnptable)):
                            if snp[3+w_idx]==None:
                                w_af+=0
                            else:

                                w_af+=snp[3+w_idx][1]
                        for d_idx in range(len(options.domesticsnptable)):
                            if snp[3+len(options.wildsnptable)+d_idx]==None:
                                d_af+=0
                            else:
                                d_af+=snp[3+len(options.wildsnptable)+d_idx][1]
                        w_af=w_af/len(options.wildsnptable)
                        d_af=d_af/len(options.domesticsnptable)
                        delta_af=abs(w_af-d_af)
#                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                        if delta_af>float(options.filterfreq):
                            snplist_posdeltaAFrefalt.append((delta_af,w_af,d_af,copy.deepcopy(snp)))
                    snplist_posdeltaAFrefalt.sort(key=lambda  listRec:listRec[0],reverse=True)
                    for delta_af,w_af,d_af,snp in snplist_posdeltaAFrefalt:
                        print('%.5f'%delta_af,w_af,d_af,*snp,sep="\t",file=intronoutfile)
#                     start=gtfMap[linelist[0]]
#                     min(int(linelist[1]),)
#                     candidateRegionbed[(linelist[0],linelist[1],linelist[2],linelist[3],linelist[4])]="unfinished"
            #search element        
#             else:
            else:

                allsweepRegionsnps_fromsql=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+linelist[1] + " and snp_pos< "+linelist[2])
                snplist_posdeltaAFrefalt=[]
                for snp in allsweepRegionsnps_fromsql:
                    w_af=0;d_af=0;delta_af=0
                    for w_idx in range(len(options.wildsnptable)):
                        if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx]==None:
                            w_af+=0
                        else:
                            ref_base=snp[2+w_idx*2]
                            w_base=snp[3+w_idx*2]
                            w_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx])
                    for d_idx in range(len(options.domesticsnptable)):
                        if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx]==None:
                            d_af+=0
                        else:
                            ref_base=snp[2+d_idx*2+len(options.wildsnptable)*2]
                            w_base=snp[3+d_idx*2+len(options.wildsnptable)*2]
                            d_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx])
                    w_af=w_af/len(options.wildsnptable)
                    d_af=d_af/len(options.domesticsnptable)
                    delta_af=abs(w_af-d_af)
#                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                    if delta_af>float(options.filterfreq):
                        snplist_posdeltaAFrefalt.append((int(snp[1]),delta_af,ref_base,w_base,w_af,d_af))
                snplist_posdeltaAFrefalt.sort(key=lambda  listRec:listRec[1],reverse=True)
                print("region:",regionline.strip(),file=intergenic)
                for snp_pos,delta_af,ref_base,w_base,w_af,d_af in snplist_posdeltaAFrefalt:
                    print(delta_af,snp_pos,ref_base,w_base,w_af,d_af,sep="\t",file=intergenic)
            if chrom in allseqtobed:
                for sstartpos,sendpos,fafilename,qs,qe,revcom,total_bases,gap_open in allseqtobed[chrom]:
                    if (sendpos>regionstart and sendpos<regionend) or (sstartpos>regionstart and sstartpos<regionend):
                        overlap_start=max(regionstart,sstartpos)
                        overlap_end=min(regionend,sendpos)
#                             sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)
#                             sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end)
#                         muscleout_seqgenerator=SeqIO.parse(elementfold+"/"+fafilename,"fasta")
#                         muscleout_seqmap={}
#                         for seq_rec in muscleout_seqgenerator:
#                             muscleout_seqmap[seq_rec.id]=seq_rec.seq
                        print(sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))#+" union "+sqlselectstatementpart_count_right+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))                      
                        allsnps=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" where chrID='"+chrom+"' and snp_pos>"+str(overlap_start) + " and snp_pos< "+str(overlap_end))
                        snplist_posdeltaAFrefalt=[]
                        for snp in allsnps:
                            w_af=0;d_af=0;delta_af=0
                            for w_idx in range(len(options.wildsnptable)):
                                if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx]==None:
                                    w_af+=0
                                else:
                                    ref_base=snp[2+w_idx*2]
                                    w_base=snp[3+w_idx*2]
                                    w_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+w_idx])
                            for d_idx in range(len(options.domesticsnptable)):
                                if snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx]==None:
                                    d_af+=0
                                else:
                                    ref_base=snp[2+d_idx*2+len(options.wildsnptable)*2]
                                    w_base=snp[3+d_idx*2+len(options.wildsnptable)*2]
                                    d_af+=float(snp[2+len(options.wildsnptable+options.domesticsnptable)*2+len(options.wildsnptable)+d_idx])
                            w_af=w_af/len(options.wildsnptable)
                            d_af=d_af/len(options.domesticsnptable)
                            delta_af=abs(w_af-d_af)
#                                 print(snp,"======",delta_af,w_af,d_af,ref_base,w_base)
                            if delta_af>=float(options.filterfreq):
                                snplist_posdeltaAFrefalt.append((int(snp[1]),delta_af,ref_base,w_base,w_af,d_af))
#                             targetspecies_baseidx=[];
#                             targetspecies_mutseq=[]
#                             slah_count_before_q=0
#                             seqbeforeqs=0                           
#                             if revcom=="revcom":
#                                 
#                                 cur_ref_pos=sendpos;qe_r=total_bases-qe+1;cur_ref_pos_p=cur_ref_pos+qe_r-1
#                                 print(fafilename,"qe_r",qe_r,"cur_ref_pos_p",cur_ref_pos_p,cur_ref_pos)
#                                 targetspecies_seq=list(Util.complementary(list(muscleout_seqmap[targetseqnamesubstr])))
#                                 targetspecies_seq.reverse()
#                                 currentChromNOlen=int(dbchromtools.operateDB("select","select chrlength from pekingduckchrominfo where chrID='"+chrom+"'")[0][0])
#                                 refseq=Util.getRefSeqBypos_faster(reffahander, refidxByChr, chrom, sstartpos, sendpos, currentChromNOlen)
#                                 alignments=pairwise2.align.globalxx("".join(refseq[chrom][1:]),"".join(targetspecies_seq).replace("-", "N"))
#                                 print(alignments[0])
#                                 alignments[0]sssss
#                                 while seqbeforeqs!=qe_r-1:
#                                     if  targetspecies_seq[seqbeforeqs+slah_count_before_q]!="-":
#                                         seqbeforeqs+=1
#                                     else:
#                                         slah_count_before_q+=1
#                                 for base in reversed(targetspecies_seq[0:qe_r+slah_count_before_q-1]):#process - and bases before qs
#                                     if re.search(r"[ATCGNnatcg]",base)!=None:
#                                         targetspecies_mutseq.append(base)
#                                         targetspecies_baseidx.append(cur_ref_pos_p);cur_ref_pos_p-=1
#                                     else:
#                                         targetspecies_mutseq.append("-")
#                                         targetspecies_baseidx.append(0)
#                             #seq collected above may be wrong,but it's not used else,so it doesn't matter
#                                 snplist_posdeltaAFrefalt.sort(key=lambda tp:tp[0]);snp_idx=0
#                                 mapsnppos_valueidx={}
#                                 num_slash=0
#                                 for base in targetspecies_seq[qe_r+slah_count_before_q-1:]:
#                                     print(cur_ref_pos,fafilename,base,file=testff)
#                                     if re.search(r"[ATCGNnatcg]",base)!=None:
#                                         if snplist_posdeltaAFrefalt[snp_idx][0]==cur_ref_pos:#find snp postion
#                                             if len(snplist_posdeltaAFrefalt[snp_idx][3])==1 and len(snplist_posdeltaAFrefalt[snp_idx][2])>1:
#                                                 num_slash=len(snplist_posdeltaAFrefalt[snp_idx][2])-1
#                                             elif len(snplist_posdeltaAFrefalt[snp_idx][3])>1 and len(snplist_posdeltaAFrefalt[snp_idx][2])>1:
#                                                 print("what's wrong with the indel ref base length")
#                                             mapsnppos_valueidx[cur_ref_pos]=len(targetspecies_mutseq)
#                                             print(cur_ref_pos,len(targetspecies_mutseq))
#                                             targetspecies_mutseq.append(snplist_posdeltaAFrefalt[snp_idx][3])
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                             snp_idx+=1  
#                                         if num_slash!=0:
#                                             num_slash-=1
#                                             targetspecies_mutseq.append("-")
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                         elif snplist_posdeltaAFrefalt[snp_idx][0]!=cur_ref_pos:
#                                             targetspecies_mutseq.append(base)# not snp position
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                         cur_ref_pos-=1
#                                                                            
#                             elif revcom=="forward":#forward
#                                 cur_ref_pos=sstartpos;cur_ref_pos_p=cur_ref_pos-qs-1
#                                 targetspecies_seq=list(muscleout_seqmap[targetseqnamesubstr])
#                                 while seqbeforeqs!=qs-1:#collection "-" before qs (or qe)                  
#                                     if targetspecies_seq[seqbeforeqs+slah_count_before_q]!="-":
#                                         seqbeforeqs+=1
#                                     else:
#                                         slah_count_before_q+=1         
#                                 for base in reversed(targetspecies_seq[0:qs+slah_count_before_q-1]):#process - and bases before qs
#                                     if re.search(r"[ATCGNnatcg]",base)!=None:
#                                         targetspecies_mutseq.append(base)
#                                         targetspecies_baseidx.append(cur_ref_pos_p);cur_ref_pos_p+=1
#                                     else:
#                                         targetspecies_mutseq.append("-")
#                                         targetspecies_baseidx.append(0)
#                             #seq collected above may be wrong,but it's not used else,so it doesn't matter
#                                 snplist_posdeltaAFrefalt.sort(key=lambda tp:tp[0]);snp_idx=0
#                                 mapsnppos_valueidx={}
#                                 num_slash=0
#                                 for base in targetspecies_seq[qs+slah_count_before_q-1:]:
#                                     print(cur_ref_pos,fafilename,base,file=testff)
#                                     if re.search(r"[ATCGNnatcg]",base)!=None:
#                                         if snplist_posdeltaAFrefalt[snp_idx][0]==cur_ref_pos:#find snp postion
#                                             if len(snplist_posdeltaAFrefalt[snp_idx][3])==1 and len(snplist_posdeltaAFrefalt[snp_idx][2])>1:
#                                                 num_slash=len(snplist_posdeltaAFrefalt[snp_idx][2])-1
#                                             elif len(snplist_posdeltaAFrefalt[snp_idx][3])>1 and len(snplist_posdeltaAFrefalt[snp_idx][2])>1:
#                                                 print("what's wrong with the indel ref base length")
#                                             mapsnppos_valueidx[cur_ref_pos]=len(targetspecies_mutseq)
#                                             print(cur_ref_pos,len(targetspecies_mutseq))
#                                             targetspecies_mutseq.append(snplist_posdeltaAFrefalt[snp_idx][3])
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                             snp_idx+=1
#                                         if num_slash!=0:
#                                             num_slash-=1
#                                             targetspecies_mutseq.append("-")
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                         elif snplist_posdeltaAFrefalt[snp_idx][0]!=cur_ref_pos:
#                                             targetspecies_mutseq.append(base)# not snp position
#                                             targetspecies_baseidx.append(cur_ref_pos)
#                                         cur_ref_pos+=1
#                                     else:
#                                         targetspecies_mutseq.append("-")
#                                         targetspecies_baseidx.append(0)
                        snplist_posdeltaAFrefalt.sort(key=lambda tp:tp[1],reverse=True)
                        print("region:",regionline.strip(),fafilename,file=outfile)
#                             print(targetspecies_seq,len(targetspecies_baseidx))
#                             print(targetspecies_mutseq,len(targetspecies_baseidx))
#                             print(targetspecies_baseidx,len(targetspecies_baseidx))
#                             print(snplist_posdeltaAFrefalt,mapsnppos_valueidx,len(snplist_posdeltaAFrefalt),len(mapsnppos_valueidx))
                        for snp_pos,delta_af,ref_base,w_base,w_af,d_af in snplist_posdeltaAFrefalt:
                            print('%.5f'%delta_af,snp_pos,ref_base,w_base,w_af,d_af,sep='\t',file=outfile)
#                                 for speciese in sorted(muscleout_seqmap.keys()):
#                                     tempseqlist=list(muscleout_seqmap[speciese])
#                                     print(tempseqlist[mapsnppos_valueidx[snp_pos]],end="\t",file=outfile)
#                                 print("",file=outfile)
                    else:
                        print("conserved region no overlap with gene ")
        f.close()
        
    intergenic.close()
    outfile.close()
    utroutfile.close()
    cdsoutfile.close()
    intronoutfile.close()

#     testf.close()
#     testff.close()
    os.system("awk '$1!~/region/&& $1 !~/notfound/{print $0}' "+options.outputpath+".utr>"+options.outputpath+".utr_")
    os.system("awk '$1!~/region/&& $1 !~/notfound/{print $0}' "+options.outputpath+".conserved>"+options.outputpath+".conserved_")
    os.system("awk '$1!~/region/&& $1 !~/notfound/{print $0}' "+options.outputpath+".intergenic>"+options.outputpath+".intergenic_")
    os.system("""awk '$1!~/region/{OFS="\\t";print $1,$2,$3,$4,$5,$6}' """+options.outputpath+".cds>"+options.outputpath+".cds_")
    os.system("""awk '$1!~/region/{OFS="\\t";print $1,$2,$3,$4,$5,$6}' """+options.outputpath+".intron>"+options.outputpath+".intron_")
    print(options.outputpath+".conserved<-read.delim('"+options.outputpath+".conserved_',header=T)")
    print(options.outputpath+".utr<-read.delim('"+options.outputpath+".utr_',header=T)")
    print(options.outputpath+".intergenic<-read.delim('"+options.outputpath+".intergenic_',header=T)")
    print(options.outputpath+".cds<-read.delim('"+options.outputpath+".cds_',header=T)")
    print(options.outputpath+".intron<-read.delim('"+options.outputpath+".intron_',header=T)")
    print("plot(density("+options.outputpath+".conserved$delta_AF),main="+"'"+options.outputpath+" conserved',xlim=c(0,1))")
    print("plot(density("+options.outputpath+".utr$delta_AF),main="+"'"+options.outputpath+" utr',xlim=c(0,1))")
    print("plot(density("+options.outputpath+".intron$delta_AF),main="+"'"+options.outputpath+" intron',xlim=c(0,1))")
    print("plot(density("+options.outputpath+".cds$delta_AF),main="+"'"+options.outputpath+" cds',xlim=c(0,1))")
    print("plot(density("+options.outputpath+".intergenic$delta_AF),main="+"'"+options.outputpath+" intergenic',xlim=c(0,1))")
    