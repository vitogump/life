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
    def __init__(self,considerINDEL="no"):
        self.considerINDEL=considerINDEL.lower()
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

        if refdep<=seqerrorrate*(refdep+altalleledep):
            return
#        if refdep+altalleledep<10:
#            return
        self.COUNTED+=1
    def getResult(self):
        snpsinthiswin=self.COUNTED
        self.COUNTED=0
        return snpsinthiswin
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
        
class Caculate_Hp(Caculator):
    def __init__(self,minsnps=3):
        super().__init__()
        self.minsnps=minsnps
        self.COUNTED=0
        self.CNMI = 0
        self.CNMA = 0
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
                            
            

        if refdep<=seqerrorrate*(refdep+altalleledep):
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
        if self.COUNTED<=self.minsnps:
            HETEROZY= 'NA'
        self.CNMA = 0
        self.CNMI = 0
        self.COUNTED=0
        return HETEROZY
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
        return ([a/self.winsize for a in countlist],[a/self.winsize for a in average])
class Caculate_Fst(Caculator):
    def __init__(self,minsnps=3):
        super().__init__()
        self.minsnps=minsnps
        self.CNk = 0
        self.CDk = 0
        self.COUNTED=0
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
                             
        if refdep_1<=seqerrorrate*(refdep_1+altalleledep_1) or refdep_2<=seqerrorrate*(refdep_2+altalleledep_2):
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
        if self.COUNTED<=self.minsnps:
            Fst='NA'
        self.CDk = 0
        self.CNk = 0
        self.COUNTED=0
        return Fst