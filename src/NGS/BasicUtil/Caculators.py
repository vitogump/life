import re,copy
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
    def __init__(self,winwidth,considerINDEL="no"):
        self.considerINDEL=considerINDEL.lower()
        self.winwidth=winwidth
        self.COUNTED=0
    def process(self, T,seqerrorrate=0.01):
        if self.considerINDEL=="no" and (len(T[1])!=1 or len(T[2])!=1):
            return
        if self.considerINDEL=="just" and (len(T[1])==1 and len(T[2])==1):
            return
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T[3])
        refdep=0;altalleledep=0
        if dp4!=None:#vcf from samtools 
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
        else:
#             if len(T[1])!=1 or len(T[2])!=1:
#                 return
            AD_idx=(re.split(":",T[4])).index("AD")#gatk GT:AD:DP:GQ:PL
            for sample in T[5]:
                if len(re.split(":",sample))==1:# ./.
                    continue
                AD_depth=re.split(",",re.split(":",sample)[AD_idx])
                try :
                    refdep+=int(AD_depth[0])
                    altalleledep+=int(AD_depth[1])
                except ValueError:
                    print(sample,end="|")

        if refdep<=seqerrorrate*(refdep+altalleledep):#not fixed
            return
#        if refdep+altalleledep<10:
#            return
        self.COUNTED+=1
    def getResult(self):
        snpsinthiswin=self.COUNTED
        snpsdensity=snpsinthiswin/self.winwidth
        self.COUNTED=0
        return snpsinthiswin,snpsdensity
class Caculate_phastConsValue(Caculator):
    def __init__(self):
        super().__init__()
        self.conservationvalue=0
        self.totalPostionsAwin=0
    def process(self,T,NumOfPositions):
        self.conservationvalue+=T[2]*NumOfPositions
        self.totalPostionsAwin+=NumOfPositions
    def getResult(self):
        winvalue="NA"
        if self.totalPostionsAwin==0:
            self.conservationvalue=0
            print("getResult")
            return "NA"
        else:
            winvalue=self.conservationvalue
            print(self.conservationvalue,self.totalPostionsAwin)
            winvalue=self.conservationvalue/self.totalPostionsAwin
            self.conservationvalue=0
            self.totalPostionsAwin=0
            return winvalue

class Caculate_Dstatistics(Caculator):
    def __init__(self,considerFixed=False):
        super().__init__()
        self.ABBA=0
        self.BABA=0
        self.numerator_fixed=0
        self.denominator_fixed=0
        self.numerator_snp=0
        self.denominator_snp=0
        self.considerFixed=considerFixed
        self.COUNTEDforSNP_notonlyfixed=0
    def process(self,T,seqerrorrate=0.01):
        """T:(pos,"a,b","c,d","e,f",A_base_idx)     1 - A_base_idx= B_base_idx ie. T[4] is the A idx . 1-T[4] is the B idx
        """
        p1A=int(re.split(r",",T[1])[T[4]]);p1B=int(re.split(r",",T[1])[1-T[4]])
        p2A=int(re.split(r",",T[2])[T[4]]);p2B=int(re.split(r",",T[2])[1-T[4]])
        p3A=int(re.split(r",",T[3])[T[4]]);p3B=int(re.split(r",",T[3])[1-T[4]])
        if p3A==0 and p3B!=0:#p3 fixed as B
            if p2A==0 and p2B!=0:#p2 fixed as B
                if p1B==0 and p1A!=0:#p1 fixed as A
                    self.ABBA+=1
                    print(T,"abba")
                    self.numerator_fixed+=1
                    self.denominator_fixed+=1
            elif p2B==0 and p2A!=0:#p2 fixed as A
                if p1A==0 and p1B!=0:#p1 fixed as B
                    self.BABA+=1
                    print(T,"baba")
                    self.numerator_fixed+=-1
                    self.denominator_fixed+=1
        if (not self.considerFixed) and( p1A==0 or p1B==0 or p2A==0 or p2B==0 or p3A==0 or p3B==0):
            return
        try:
            self.numerator_snp+=p3B/(p3B+p3A) * ((p1A/(p1A+p1B))*(p2B/(p2A+p2B)) - (p1B/(p1A+p1B))*(p2A/(p2A+p2B)))
            self.denominator_snp+=p3B/(p3B+p3A) * ((p1A/(p1A+p1B))*(p2B/(p2A+p2B)) + (p1B/(p1A+p1B))*(p2A/(p2A+p2B)))
            self.COUNTEDforSNP_notonlyfixed+=1
        except ZeroDivisionError:
            print(self.denominator_snp,self.numerator_snp,T)
    def getResult(self):
        ABBAcount=self.ABBA
        BABAcount=self.BABA
        try:
            D_fixed=self.numerator_fixed/self.denominator_fixed
            D_snp=self.numerator_snp/self.denominator_snp
        except ZeroDivisionError:
            D_fixed='NA'
            D_snp='NA'
        self.numerator_fixed=0;self.denominator_fixed=0;self.ABBA=0;self.BABA=0
        self.numerator_snp=0;self.denominator_snp=0;noofsnps=copy.deepcopy(self.COUNTEDforSNP_notonlyfixed);self.COUNTEDforSNP_notonlyfixed=0
        return ABBAcount,BABAcount,D_fixed,D_snp,noofsnps
        
        
        
class Caculate_Hp(Caculator):
    def __init__(self,minsnps=3,considerFixed=False):
        super().__init__()
        self.minsnps=minsnps
        self.COUNTED=0
        self.CNMI = 0
        self.CNMA = 0
        self.considerFixed=considerFixed
    def process(self, T,seqerrorrate=0.01):
        
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T[3])
        refdep=0;altalleledep=0
        if dp4!=None:#vcf from samtools 
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
        else:
            AD_idx=(re.split(":",T[4])).index("AD")#gatk GT:AD:DP:GQ:PL
            for sample in T[5]:
                if len(re.split(":",sample))==1:# ./.
                    continue
                AD_depth=re.split(",",re.split(":",sample)[AD_idx])
                try :
                    refdep+=int(AD_depth[0])
                    altalleledep+=int(AD_depth[1])
                except ValueError:
                    print(sample,end="|")
                            
            

        if (not self.considerFixed) and refdep<=seqerrorrate*(refdep+altalleledep):#not fixed
            return
        if refdep+altalleledep<10:
            return
        self.COUNTED+=1
        if refdep < altalleledep:
            self.CNMI += refdep
            self.CNMA += altalleledep
        else:
            self.CNMA += refdep
            self.CNMI += altalleledep
    def getResult(self):
        HETEROZY = 'NA'
        try:
            HETEROZY = self.CNMA * self.CNMI * 2 / ((self.CNMA + self.CNMI) ** 2)
        except ZeroDivisionError:
            #print("the Heterozigosity value of currentwindow is dividsion by zero,so set it to be NA")
            HETEROZY = 'NA'
        
        noofsnpcount=self.COUNTED
        print(noofsnpcount)
#         if self.COUNTED<=self.minsnps:
#             HETEROZY= 'NA'
        self.CNMA = 0
        self.CNMI = 0
        self.COUNTED=0
        print(noofsnpcount)
        return noofsnpcount,HETEROZY
class Caculate_depth_judge(Caculator):
    def __init__(self,sampleNo,winsize,mindepth):
        self.mindepth=int(mindepth)
        self.sampleNo=sampleNo
        self.winsize=winsize
        self.COVERED_COUNT=[0]*sampleNo
        self.AVERAGE_DEPTH=[0]*sampleNo
    def process(self,T,seqerrorrate=0.01):
        """
        T=(pos,sample1dp,sample2dp,,,,,,)
        """
#         print(T,"\n",self.AVERAGE_DEPTH)
        for sampleNo in range(1,len(T)):
            self.AVERAGE_DEPTH[sampleNo-1]+=int(T[sampleNo])
            if int(T[sampleNo])>=self.mindepth:
                self.COVERED_COUNT[sampleNo-1]+=1
            
    def getResult(self):
        """pecentage of cover,average depth
        """
        countlist=copy.deepcopy(self.COVERED_COUNT);average =  copy.deepcopy(self.AVERAGE_DEPTH)
        del self.AVERAGE_DEPTH[:]
        del self.COVERED_COUNT[:]
        self.COVERED_COUNT=[0]*self.sampleNo
        self.AVERAGE_DEPTH=[0]*self.sampleNo
        return "empty",([a/self.winsize for a in countlist],[a/self.winsize for a in average])
class Caculate_Fst(Caculator):
    def __init__(self,minsnps=3,considerFixed=False):
        super().__init__()
        self.minsnps=minsnps
        self.CNk = 0
        self.CDk = 0
        self.COUNTED=0
        self.considerFixed=considerFixed
    def process(self, T,seqerrorrate=0.01):
        
        refdep_1=0;refdep_2=0
        altalleledep_1=0;altalleledep_2=0
#        T1 = (T[0], T[1], T[2], T[3])
#        T2 = (T[4], T[5], T[6], T[7])
        pop1=T[3]
        pop2=T[4]
        dp4_1 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop1[0])
        if dp4_1!=None:#vcf from samtools
            refdep_1 = int(dp4_1.group(1)) + int(dp4_1.group(2))
            altalleledep_1 = int(dp4_1.group(3)) + int(dp4_1.group(4))
            dp4_2 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop2[0])
            refdep_2 = int(dp4_2.group(1)) + int(dp4_2.group(2))
            altalleledep_2 = int(dp4_2.group(3)) + int(dp4_2.group(4))
        else:#vcf from gatk
            AD_idx_1=(re.split(":",pop1[1])).index("AD")#gatk GT:AD:DP:GQ:PL
            AD_idx_2=(re.split(":",pop2[1])).index("AD")
            for sample in pop1[2][:]:
                if len(re.split(":",sample))==1:# ./.
                    continue
                AD_depth=re.split(",",re.split(":",sample)[AD_idx_1])
                try:
                    refdep_1+=int(AD_depth[0])
                    altalleledep_1+=int(AD_depth[1])
                except ValueError:
                    print(sample,end="|")

            for sample in pop2[2][:]:
                if len(re.split(":",sample))==1:# ./.
                    continue
                AD_depth=re.split(",",re.split(":",sample)[AD_idx_2])
                try:
                    refdep_2+=int(AD_depth[0])
                    altalleledep_2+=int(AD_depth[1])
                except ValueError:
                    print(sample,end="|")
                             
        if (not self.considerFixed) and ( refdep_1<=seqerrorrate*(refdep_1+altalleledep_1) or refdep_2<=seqerrorrate*(refdep_2+altalleledep_2)):
            return  #NOTICT HERE
        self.COUNTED+=1
        h_1 = refdep_1 * altalleledep_1 / ((refdep_1 + altalleledep_1 - 1) * (refdep_1 + altalleledep_1))
        h_2 = refdep_2 * altalleledep_2 / ((refdep_2 + altalleledep_2 - 1) * (refdep_2 + altalleledep_2))
        Nk = ((refdep_1 / (refdep_1 + altalleledep_1) - refdep_2 / (refdep_2 + altalleledep_2)) ** 2 - h_1 / (refdep_1 + altalleledep_1) - h_2 / (refdep_2 + altalleledep_2))
        self.CNk += Nk
        self.CDk += (Nk + h_1 + h_2)
    def getResult(self):
        Fst = 'NA'
        try:
            Fst = self.CNk / self.CDk
        except ZeroDivisionError:
            #print("the Fst value of currentwindow is dividsion by zero,so set it to be NA")
            Fst = 'NA'
#         if self.COUNTED<=self.minsnps:
#             Fst='NA'
        self.CDk = 0
        self.CNk = 0
        noofsnp=self.COUNTED
        self.COUNTED=0
        return noofsnp,Fst