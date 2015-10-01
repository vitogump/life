import re, copy,math,numpy
'''
Created on 2013-7-2

@author: rui
'''
class Caculator():
    def process(self, T):
        pass
    def getResult(self):
        pass
class Caculate_SNPsPerBIN(Caculator):
    def __init__(self, winwidth, considerINDEL="no", MethodToSeq="pool"):
        self.considerINDEL = considerINDEL.lower()
        self.winwidth = winwidth
        self.COUNTED = 0
        self.MethodToSeq = MethodToSeq
    def process(self, T, seqerrorrate=0.01):
        if self.considerINDEL == "no" and (len(T[1]) != 1 or len(T[2]) != 1):
            return
        if self.considerINDEL == "just" and (len(T[1]) == 1 and len(T[2]) == 1):
            return
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T[3])
        refdep = 0;altalleledep = 0
        if dp4 != None:  # vcf from samtools 
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
        else:
            if self.MethodToSeq == "pool":
                AD_idx = (re.split(":", T[4])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[5]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])
                        altalleledep += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")

                if refdep <= seqerrorrate * (refdep + altalleledep):  # not fixed
                    return
            elif self.MethodToSeq == "indvd":
                AN = int(re.search(r"AN=(\d+);", T[3]).group(1))
                AF = int(float(re.search(r"AF=([\d\.]+);", T[3]).group(1)))
                if AF == 1:
                    return
#                 refdep=AN-AC
#                 altalleledep=AC                
#        if refdep+altalleledep<10:
#            return
        self.COUNTED += 1
    def getResult(self):
        snpsinthiswin = self.COUNTED
        snpsdensity = 1000 * snpsinthiswin / self.winwidth
        self.COUNTED = 0
        return snpsinthiswin, snpsdensity
class Caculate_phastConsValue(Caculator):
    def __init__(self):
        super().__init__()
        self.conservationvalue = 0
        self.totalPostionsAwin = 0
    def process(self, T, NumOfPositions):
        self.conservationvalue += T[2] * NumOfPositions
        self.totalPostionsAwin += NumOfPositions
    def getResult(self):
        winvalue = "NA"
        if self.totalPostionsAwin == 0:
            self.conservationvalue = 0
            print("getResult")
            return "NA"
        else:
            winvalue = self.conservationvalue
            print(self.conservationvalue, self.totalPostionsAwin)
            winvalue = self.conservationvalue / self.totalPostionsAwin
            self.conservationvalue = 0
            self.totalPostionsAwin = 0
            return winvalue

class Caculate_Dstatistics(Caculator):
    def __init__(self, considerFixed=False):
        super().__init__()
        self.ABBA = 0
        self.BABA = 0
        self.numerator_fixed = 0
        self.denominator_fixed = 0
        self.numerator_snp = 0
        self.denominator_snp = 0
        self.considerFixed = considerFixed
        self.COUNTEDforSNP_notonlyfixed = 0
    def process(self, T, seqerrorrate=0.01):
        """T:(pos,"a,b","c,d","e,f",A_base_idx)     1 - A_base_idx= B_base_idx ie. T[4] is the A idx . 1-T[4] is the B idx
        """
        p1A = int(re.split(r",", T[1])[T[4]]);p1B = int(re.split(r",", T[1])[1 - T[4]])
        p2A = int(re.split(r",", T[2])[T[4]]);p2B = int(re.split(r",", T[2])[1 - T[4]])
        p3A = int(re.split(r",", T[3])[T[4]]);p3B = int(re.split(r",", T[3])[1 - T[4]])
        if p3A == 0 and p3B != 0:  # p3 fixed as B
            if p2A == 0 and p2B != 0:  # p2 fixed as B
                if p1B == 0 and p1A != 0:  # p1 fixed as A
                    self.ABBA += 1
                    print(T, "abba")
                    self.numerator_fixed += 1
                    self.denominator_fixed += 1
            elif p2B == 0 and p2A != 0:  # p2 fixed as A
                if p1A == 0 and p1B != 0:  # p1 fixed as B
                    self.BABA += 1
                    print(T, "baba")
                    self.numerator_fixed += -1
                    self.denominator_fixed += 1
        if (not self.considerFixed) and(p1A == 0 or p1B == 0 or p2A == 0 or p2B == 0 or p3A == 0 or p3B == 0):
            return
        try:
            self.numerator_snp += p3B / (p3B + p3A) * ((p1A / (p1A + p1B)) * (p2B / (p2A + p2B)) - (p1B / (p1A + p1B)) * (p2A / (p2A + p2B)))
            self.denominator_snp += p3B / (p3B + p3A) * ((p1A / (p1A + p1B)) * (p2B / (p2A + p2B)) + (p1B / (p1A + p1B)) * (p2A / (p2A + p2B)))
            self.COUNTEDforSNP_notonlyfixed += 1
        except ZeroDivisionError:
            print(self.denominator_snp, self.numerator_snp, T)
    def getResult(self):
        ABBAcount = self.ABBA
        BABAcount = self.BABA
        try:
            D_fixed = self.numerator_fixed / self.denominator_fixed
            D_snp = self.numerator_snp / self.denominator_snp
        except ZeroDivisionError:
            D_fixed = 'NA'
            D_snp = 'NA'
        self.numerator_fixed = 0;self.denominator_fixed = 0;self.ABBA = 0;self.BABA = 0
        self.numerator_snp = 0;self.denominator_snp = 0;noofsnps = copy.deepcopy(self.COUNTEDforSNP_notonlyfixed);self.COUNTEDforSNP_notonlyfixed = 0
        return ABBAcount, BABAcount, D_fixed, D_snp, noofsnps
        
        

class Caculate_Hp(Caculator):
    def __init__(self, SeqMethodlist=["pool"], minsnps=10):
        super().__init__()
        self.minsnps = minsnps
        self.COUNTED = [0] * len(SeqMethodlist)
        self.CNMI = [0] * len(SeqMethodlist)
        self.CNMA = [0] * len(SeqMethodlist)
        self.sum_mean_2pq = 0
        self.SeqMethodlist = SeqMethodlist     
    def process(self, T, seqerrorrate=0.005, mode=1):
        if len(T[1]) != len(T[2]) or len(T[2])!=1 or len(T[2])!=1:
            return
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            MethofToSeq = self.SeqMethodlist[MethodToSeq_idx]
            if T[3 + MethodToSeq_idx] == None:
                continue
            if MethofToSeq == "pool":
                refdep = 0;altalleledep = 0
                AD_idx = (re.split(":", T[3 + MethodToSeq_idx][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[3 + MethodToSeq_idx][2]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])
                        altalleledep += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")
            elif MethofToSeq == "indvd":
                AF = float(re.search(r"AF=([\d\.]+);", T[3 + MethodToSeq_idx][0]).group(1))
                AN = int(re.search(r"AN=(\d+);", T[3 + MethodToSeq_idx][0]).group(1))
                AC = int(re.search(r"AC=(\d+);", T[3 + MethodToSeq_idx][0]).group(1))
                refdep = AN - AC
                altalleledep = AC
            if refdep <= seqerrorrate * (refdep + altalleledep):  # skip fixed as altallele ,ie refdep == 0
                continue
            if refdep + altalleledep < 10:
                continue
            self.COUNTED[MethodToSeq_idx] += 1
            if refdep < altalleledep:
                self.CNMI[MethodToSeq_idx] += refdep
                self.CNMA[MethodToSeq_idx] += altalleledep
            else:
                self.CNMA[MethodToSeq_idx] += refdep
                self.CNMI[MethodToSeq_idx] += altalleledep
    def getResult(self):
        HETEROZY = ['NA'] * len(self.SeqMethodlist)
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            try:
                HETEROZY[MethodToSeq_idx] = self.CNMA[MethodToSeq_idx] * self.CNMI[MethodToSeq_idx] * 2 / ((self.CNMA[MethodToSeq_idx] + self.CNMI[MethodToSeq_idx]) ** 2)
            except ZeroDivisionError:
                # print("the Heterozigosity value of currentwindow is dividsion by zero,so set it to be NA")
                HETEROZY[MethodToSeq_idx] = 'NA'
        het_count = 0;het_sum = 0;pop_idx=0
        for pop_idx in range(len(HETEROZY)) :
            if HETEROZY[pop_idx] != 'NA' and self.COUNTED[pop_idx]>=self.minsnps:
                het_count += 1
                het_sum += HETEROZY[pop_idx]
#                 print(HETEROZY[pop_idx],end="\t")
#         print()
        noofsnpcount = min(self.COUNTED)
        if het_count == 0 :
            HETEROZY_toreturn = 'NA'
        else:
            HETEROZY_toreturn = het_sum / het_count
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        return noofsnpcount, HETEROZY_toreturn
class Caculate_depth_judge(Caculator):
    def __init__(self, total_samples, winsize, mindepth, speciesorder=[], sampleidxlisttocount={}):
        self.mindepth = int(mindepth)
        self.total_samples = total_samples
        self.winsize = winsize
        if speciesorder == [] and (not sampleidxlisttocount):
            self.COVERED_COUNT = [0] * total_samples
            self.AVERAGE_DEPTH = [0] * total_samples
            self.speciesorder = None;
            self.sampleidxlisttocount = None
        elif len(speciesorder) != 0 and len(sampleidxlisttocount.keys()) != 0:
            self.COVERED_COUNT = [0] * len(speciesorder)
            self.AVERAGE_DEPTH = [0] * len(speciesorder)
            self.speciesorder = speciesorder
            self.sampleidxlisttocount = sampleidxlisttocount
    def process(self, T, seqerrorrate=0.01):
        """
        T=(pos,sample1dp,sample2dp,,,,,,)
        """
#         print(T,"\n",self.AVERAGE_DEPTH)
        if self.speciesorder == [] and (not self.sampleidxlisttocount):
            for sampleidx in range(1, len(T)):
                self.AVERAGE_DEPTH[sampleidx - 1] += int(T[sampleidx])
                if int(T[sampleidx]) >= self.mindepth:
                    self.COVERED_COUNT[sampleidx - 1] += 1
        elif len(self.speciesorder) != 0 and len(self.sampleidxlisttocount.keys()) != 0:
            for species in self.speciesorder:
                totaldepth = 0
                for sampleidx in self.sampleidxlisttocount[species]:
                    totaldepth += int(T[sampleidx])
                self.AVERAGE_DEPTH[self.speciesorder.index(species)] += totaldepth
                if totaldepth >= self.mindepth:
                    self.COVERED_COUNT[self.speciesorder.index(species)] += 1
                
            
    def getResult(self):
        """pecentage of cover,average depth
        """
        countlist = copy.deepcopy(self.COVERED_COUNT);average = copy.deepcopy(self.AVERAGE_DEPTH)
        del self.AVERAGE_DEPTH[:]
        del self.COVERED_COUNT[:]
        if self.speciesorder == [] and (not self.sampleidxlisttocount):
            self.COVERED_COUNT = [0] * self.total_samples
            self.AVERAGE_DEPTH = [0] * self.total_samples
        elif len(self.speciesorder) != 0 and len(self.sampleidxlisttocount.keys()) != 0:
            self.COVERED_COUNT = [0] * len(self.speciesorder)
            self.AVERAGE_DEPTH = [0] * len(self.speciesorder)
        return "empty", ([a / self.winsize for a in countlist], [a / self.winsize for a in average])
class Caculate_S_ObsExp_difference(Caculator):
    def __init__(self,mindepthtojudefixed,N_of_targetpop,N_of_refpop,dbvariantstoolstojudgeancestral,topleveltablejudgeancestralname):
        super().__init__()
        self.dbvariantstoolstojudgeancestral=dbvariantstoolstojudgeancestral
        self.topleveltablejudgeancestralname=topleveltablejudgeancestralname
        self.MethodToSeqpoplist=[]
        self.mindepthtojudefixed=20
        self.N_of_targetpop=N_of_targetpop
        self.N_of_refpop=N_of_refpop
        self.depthobjlist=[]
        self.species_idx_list=[]
        self.currentchrID=None
        self.COUNT=0
        self.obsseq=[]
        self.CEXP=0
        self.CfixedDerived=0
        self.freq_xaxisKEY_yaxisVALUERelation=None
    def process(self,T):
        """T=[pos,ref,alt,pop1,pop2,.....,popn]"""
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        snp=self.dbvariantstoolstojudgeancestral.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+self.currentchrID+"' and snp_pos='"+str(T[0])+"'")
        if not snp:
#             print(self.currentchrID,T,"snp not find,skip")
            return
        else:
            A_base_idx=100
            fanyadepthlist=re.split(r",",snp[0][9])
            if len(fanyadepthlist)==2 and int(fanyadepthlist[1]) >=self.mindepthtojudefixed and fanyadepthlist[0].strip()=="0":
                A_base_idx=1
            elif len(fanyadepthlist)==2 and int(fanyadepthlist[0])>=self.mindepthtojudefixed and fanyadepthlist[1].strip()=="0":
                A_base_idx=0
            else:
                print("skip snp",snp[0][1],snp[0][7],snp[0][9],snp[0][11],snp[0][13])
                return
#             depthlist1=re.split(r",",snp[0][7])
#             depthlist2=re.split(r",",snp[0][9])
#             if len(depthlist1)==2 and len(depthlist2)==2 and (int(depthlist1[0]) + int(depthlist1[1])>=self.mindepthtojudefixed or int(depthlist2[0]) + int(depthlist2[1])>=self.mindepthtojudefixed) and ((depthlist1[0].strip()=="0" and depthlist2[0].strip()=="0") or (depthlist1[1].strip()=="0" and depthlist2[1].strip()=="0") ):
#                 if depthlist1[0].strip()=="0" and depthlist2[0].strip()=="0":
#                     A_base_idx=1
#                 elif depthlist1[1].strip()=="0" and depthlist2[1].strip()=="0":
#                     A_base_idx=0
#                 else:
#                     print(snp,"never get here!")
#             elif (len(depthlist1)==2 and  snp[0][9] == "no covered" and int(depthlist1[0]) + int(depthlist1[1])>=self.mindepthtojudefixed and (depthlist1[0].strip()=="0" or depthlist1[1].strip()=="0" ))   or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=self.mindepthtojudefixed and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
#                 if (snp[0][9] == "no covered" and depthlist1[0].strip()=="0") or (snp[0][7]=="no covered" and depthlist2[0].strip()=="0"):
#                     A_base_idx=1
#                 elif (snp[0][9] == "no covered" and depthlist1[1].strip()=="0") or (snp[0][7]=="no covered" and depthlist2[1].strip()=="0"):
#                     A_base_idx=0
#                 else:
#                     print(snp,"never get here!")
#             else:
# #                 print(snp,"skip snp")
#                 return
        ancestrallcontext=snp[0][5].strip()[0].upper()+snp[0][3+A_base_idx].strip().upper()+snp[0][5].strip()[2].upper()
        if "CG" in ancestrallcontext or "GC" in ancestrallcontext:
#             print("skip CG site",ancestrallcontext)
            return
        ##########x-axis
        countedAF=0;target_DAF_sum=0
        for tpopidx in range(3,self.N_of_targetpop+3):
            if T[tpopidx]==None:
                if self.depthobjlist==[]:
                    print("skip this pos",T)
                    continue
                else:
                    depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for idx in self.species_idx_list[tpopidx-3][:]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
                        continue
            else:
                if self.MethodToSeqpoplist[tpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.]+);", T[tpopidx][0]).group(1))
                    AN = float(re.search(r"AN=([\d]+);", T[tpopidx][0]).group(1))
                    if AN<5:
                        continue
                elif self.MethodToSeqpoplist[tpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[tpopidx][1])).index("AD")
                    for sample in T[tpopidx][2]:
                        if len(re.split(":", sample)) == 1:  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
            if A_base_idx==0:
                DAF=1-AF
            elif A_base_idx==1:
                DAF=AF
            target_DAF_sum+=DAF;countedAF+=1
        if target_DAF_sum==0 or countedAF==0:
#             print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",T,snp)
            return
        target_DAF=target_DAF_sum/countedAF
        #########y-axis
        countedAF=0;rer_DAF_sum=0
        for rpopidx in range(3+self.N_of_targetpop,self.N_of_refpop+self.N_of_targetpop+3):
            if T[rpopidx]==None:
                if self.depthobjlist==[]:
                    print("skip this snp",T)
                    continue
                else:
                    depth_linelist=self.depthobjlist[rpopidx-3-self.N_of_targetpop].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for idx in self.species_idx_list[rpopidx-3-self.N_of_targetpop][:]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
                        continue
            else:
                if self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="indvd":
                    AF=float(re.search(r"AF=([\d\.]+);", T[rpopidx][0]).group(1))
                elif self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[rpopidx][1])).index("AD")
                    for sample in T[rpopidx][2]:
                        if len(re.split(":",sample))==1:
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if refdep==altalleledep and altalleledep==0:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                rer_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0:
#             print("skip this snp,because it  no covered in this pos in ref pops",T,snp)
            return
        for a,b in sorted(self.freq_xaxisKEY_yaxisVALUERelation.keys()):
            if target_DAF>a and target_DAF<=b:
                self.CEXP+=self.freq_xaxisKEY_yaxisVALUERelation[(a,b)]
                break
        self.obsseq.append(rer_DAF_sum/countedAF)
        self.COUNT+=1
        if rer_DAF_sum/countedAF==1:
            self.CfixedDerived+=1
    def getResult(self):
        S1="NA"
        S2="NA"
        try:
            S1=math.log(numpy.sum(self.obsseq)/self.CEXP)
            S2=(numpy.sum(self.obsseq)-self.CEXP)/numpy.std(self.obsseq,ddof=1)
        except:
            S1="NA"
            S2="NA"
        noofsnp=self.COUNT
        self.COUNT=0
        self.CEXP=0
        self.obsseq=[]
        self.CfixedDerived=0
        if S1=="NA" and S2=="NA":
            return noofsnp,"NA"
        return noofsnp,[S1,S2]
        
                                     
                            
class Caculate_Fst(Caculator):
    def __init__(self, MethodToSeqpop1="pool", MethodToSeqpop2="indvd", minsnps=10):
        super().__init__()
        self.minsnps = minsnps
        self.considerfixdiffinfst=False
        self.CNk = 0
        self.CDk = 0
        self.COUNTED = [0,0]#fst used snp,fixed difference snp
#         self.considerFixed = considerFixed
        self.MethodToSeqpop1 = MethodToSeqpop1
        self.MethodToSeqpop2 = MethodToSeqpop2
#         self.depthforcurrentchrom=None
        self.depthobjmap=None
        self.species_idx_map=None
        self.currentchrID=None
        self.pop1_indvdsormediandepth=None#=6#when pop1 is none at a pos,and no depth information
        self.pop2_indvdsormediandepth=None#=6
    def process(self, T, seqerrorrate=0.01):
        """T=[pos,ref,alt,pop1,pop2]"""
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        refdep_1 = 0;refdep_2 = 0
        altalleledep_1 = 0;altalleledep_2 = 0
#        T=(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist))
        pop1 = T[3]
        pop2 = T[4]
#         dp4_1 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop1[0])
#         if dp4_1 != None:  # vcf from samtools
#         if None !=None:  # vcf from samtools
#             pass
#             refdep_1 = int(dp4_1.group(1)) + int(dp4_1.group(2))
#             altalleledep_1 = int(dp4_1.group(3)) + int(dp4_1.group(4))
#             dp4_2 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop2[0])
#             refdep_2 = int(dp4_2.group(1)) + int(dp4_2.group(2))
#             altalleledep_2 = int(dp4_2.group(3)) + int(dp4_2.group(4))
#         else:  # vcf from gatk
        if self.MethodToSeqpop1 == "pool":              
            if pop1==None:
                if self.depthobjmap==None:
                    refdep_1=self.pop1_indvdsormediandepth
                    altalleledep_1=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop1_ref"].getdepthByPos_optimized(self.currentchrID,T[0])#  re.split(r"\t",self.depthforcurrentchrom["vcfpop1_ref"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop1_ref"]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>=self.pop1_indvdsormediandepth:
                        refdep_1=self.pop1_indvdsormediandepth
                        altalleledep_1=0
                    else:
                        return
            else:
                AD_idx_1 = (re.split(":", pop1[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in pop1[2][:]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx_1])
                    try:
                        refdep_1 += int(AD_depth[0])
                        altalleledep_1 += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")
        elif self.MethodToSeqpop1 == "indvd":
            if pop1==None:
                if self.depthobjmap==None:
                    refdep_1=self.pop1_indvdsormediandepth
                    altalleledep_1=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop1_ref"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop1_ref"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop1_ref"]:
                        if int(depth_linelist[idx])>4:
                            sum_depth+=1
#                         sum_depth+=int(depth_linelist[idx])
#                     if sum_depth>=self.pop1_indvdsormediandepth:
                    refdep_1=sum_depth*2
                    altalleledep_1=0
#                     else:
#                         return
            else:
                AN = int(re.search(r"AN=(\d+);", pop1[0]).group(1))
                AC = int(re.search(r"AC=(\d+);", pop1[0]).group(1))
                refdep_1 = AN - AC
                altalleledep_1 = AC
        if self.MethodToSeqpop2 == "pool":
            
            if pop2==None:
                if self.depthobjmap==None:
                    refdep_2=self.pop2_indvdsormediandepth
                    altalleledep_2=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop2"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop2"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop2"]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>=self.pop2_indvdsormediandepth:
                        refdep_2=self.pop2_indvdsormediandepth
                        altalleledep_2=0
                    else:
                        return
            else:
                AD_idx_2 = (re.split(":", pop2[1])).index("AD")
                for sample in pop2[2][:]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx_2])
                    try:
                        refdep_2 += int(AD_depth[0])
                        altalleledep_2 += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")
        elif self.MethodToSeqpop2 == "indvd":
            if pop2==None:
                if self.depthobjmap==None:
                    refdep_2=self.pop2_indvdsormediandepth
                    altalleledep_2=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop2"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop2"][int(T[0])-1])

                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop2"]:
                        if int(depth_linelist[idx])>4:
                            sum_depth+=1
#                         sum_depth+=int(depth_linelist[idx])
#                     if sum_depth>=self.pop2_indvdsormediandepth:
                    refdep_2=sum_depth*2
                    altalleledep_2=0
#                     else:
#                         return
            else:
                AN = int(re.search(r"AN=(\d+);", pop2[0]).group(1))
                AC = int(re.search(r"AC=(\d+);", pop2[0]).group(1))
                refdep_2 = AN - AC
                altalleledep_2 = AC

              
                
#         if  (refdep_1 <= seqerrorrate * (refdep_1 + altalleledep_1) or refdep_2 <= seqerrorrate * (refdep_2 + altalleledep_2)):
#             return  # NOTICT HERE
        if refdep_1==0 and refdep_2==0:#skip both fixed as alt
            return
        if ((refdep_1 + altalleledep_1 - 1) * (refdep_1 + altalleledep_1))==0 or  ((refdep_2 + altalleledep_2 - 1) * (refdep_2 + altalleledep_2))==0:
            return
        if (refdep_1==0 and altalleledep_2==0 and altalleledep_1>=self.pop1_indvdsormediandepth and refdep_2>=self.pop2_indvdsormediandepth) or (altalleledep_1==0 and refdep_2==0 and refdep_1>=self.pop1_indvdsormediandepth and altalleledep_2>=self.pop2_indvdsormediandepth):#fixed difference
            self.COUNTED[1]+=1
            if self.considerfixdiffinfst:
                print(T,"fixdifferent in Fst")
                pass
            else:
                print(T,"fixdiffernet not in Fst")
                return
        self.COUNTED[0] += 1
        h_1 = refdep_1 * altalleledep_1 / ((refdep_1 + altalleledep_1 - 1) * (refdep_1 + altalleledep_1))
        h_2 = refdep_2 * altalleledep_2 / ((refdep_2 + altalleledep_2 - 1) * (refdep_2 + altalleledep_2))
        Nk = ((refdep_1 / (refdep_1 + altalleledep_1) - refdep_2 / (refdep_2 + altalleledep_2)) ** 2 - h_1 / (refdep_1 + altalleledep_1) - h_2 / (refdep_2 + altalleledep_2))
        self.CNk += Nk
        self.CDk += (Nk + h_1 + h_2)
#         print("self.CNk",self.CNk,"NK",Nk,"self.CDk",self.CDk)
    def getResult(self):
        Fst = 'NA'
        try:
            Fst = self.CNk / self.CDk
        except ZeroDivisionError:
            # print("the Fst value of currentwindow is dividsion by zero,so set it to be NA")
            Fst = 'NA'
#         if self.COUNTED<=self.minsnps:
#             Fst='NA'
        self.CDk = 0
        self.CNk = 0
        noofsnp = copy.copy(self.COUNTED[0])
        nooffixdifference=copy.copy(self.COUNTED[1])
        self.COUNTED = [0,0]
        if noofsnp<self.minsnps:
            Fst="NA"
        return [noofsnp,nooffixdifference], Fst
