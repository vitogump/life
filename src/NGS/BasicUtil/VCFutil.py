# -*- coding: UTF-8 -*-
'''
Created on 2013-6-30

@author: rui
'''


import random
import re, pickle, copy


class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfIndexMap = {}
        self.chromOrder=[]
        self.NumOfRecbychromOrder=[]
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        self.chromOrder=self.VcfIndexMap["chromOrder"]
        self.NumOfRecbychromOrder=self.VcfIndexMap["NumOfRecbychromOrder"]
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')
        vcfChromIndex = {}
        chromOrder=[]
        NumOfRecbychromOrder=[]
        line = vcffile.readline()
        vcfChromIndex["header"]=[line]
        while re.search(r'^##', line) != None:
            line = vcffile.readline()
            vcfChromIndex["header"].append(line)
        
        if re.search(r'^#CHROM', line) != None:
            vcfChromIndex["title"] = re.split(r'\s+', line.strip())
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)        
        currentChrom = "temptodele"
        lastPosition = vcffile.tell()
        lastChromend_currentChromstartPostion=lastPosition
#         vcfChromIndex[currentChrom]=(lastPosition,0)
        print(line)
        line = vcffile.readline()
        print(line)
        i=0
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != linelist[0]:
                chromOrder.append(currentChrom)
                NumOfRecbychromOrder.append(i);i=0
                vcfChromIndex[currentChrom]=(lastChromend_currentChromstartPostion,lastPosition)
                lastChromend_currentChromstartPostion=lastPosition
                currentChrom = linelist[0]
                
#                 vcfChromIndex[currentChrom] = (lastPosition,0)
            
            lastPosition = vcffile.tell()

            line = vcffile.readline()
            i+=1
        else:
            chromOrder.append(currentChrom)
            NumOfRecbychromOrder.append(i-1)
            vcfChromIndex[currentChrom]=(lastChromend_currentChromstartPostion,lastPosition)

        vcfChromIndex.pop("temptodele")
        i=chromOrder.index("temptodele")
        if i!=0:
            print("wrong indexVCF")
            exit(-1)
        b=chromOrder.pop(i)
        a=NumOfRecbychromOrder.pop(i)
        print(i,a,b)
        vcfChromIndex["chromOrder"]=chromOrder
        
        vcfChromIndex["NumOfRecbychromOrder"]=NumOfRecbychromOrder
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
#     def extractVcfRecByChroms(self,vcfFileName,chromlist,replacechromlist,outfile):
#         vcfFile = open(vcfFileName, 'r')
#         if len(chromlist)!=len(replacechromlist):
#             print("the length of the chromlist and replacechromlist should be the same")
#             return
#         for chrom in chromlist:
#             vcfFile.seek(self.VcfIndexMap[chrom][0])
    @staticmethod        
    def Vcf2Ped(vcfFileName,outputfileprefix,software,VcfIndexMap=None,withheader=False):
        vcffile = open(vcfFileName, "r")
        mapfile = open(outputfileprefix+".map", "w")
        pedfile = open(outputfileprefix+".ped", "w")
        positionlist=[]
        pedmap={}
        if withheader:
            line = vcffile.readline()
            while re.search(r"^##", line) != None:
                line = vcffile.readline()
            title=re.split(r"\s+",line.strip())
            total_individ= len(title) -9
            print(title,len(title),total_individ)
            for outName in title[len(title)-total_individ:]:
                pedmap[outName]=[]
        else:
            title=VcfIndexMap["title"]
            total_individ=len(VcfIndexMap["title"])-9
            print(VcfIndexMap["title"],len(VcfIndexMap["title"]),total_individ)
            for outName in VcfIndexMap["title"][len(VcfIndexMap["title"])-total_individ:]:
                pedmap[outName]=[]

        currentChromSome=None
        print("exclude those sites which ref is not N ,INDEL ,or multiple alleles")
        for line in vcffile:
            linelist = re.split(r"\s+",line)
            if linelist[3].strip().upper()=='N' or len(linelist[3].strip()) > 1 or len(linelist[4].strip())>1:#when ref is N ,or INDEL ,or multiple allels 
                continue
            
            positionlist.append((linelist[0].replace("scaffold",""),linelist[0]+"_"+linelist[1],0,linelist[1]))
            if software.upper()=="GATK":
                GT_idx=(re.split(":",linelist[8])).index("GT")#gatk GT:AD:DP:GQ:PL
                PL_idx=(re.split(":",linelist[8])).index("PL")
                for i in range(total_individ):
                    sample=linelist[i+9]
                    if len(re.split(":",sample))==1 or re.split(":",sample)[GT_idx]=="./." or  len(re.split(r",",re.split(":",sample)[PL_idx]))!=3:# ./.
                        pl="0,0,0"
                    else:
                        pl=re.split(":",sample)[PL_idx]
                        genotype = re.split(":",sample)[GT_idx]                
                    if pl!="0,0,0":
                        a1=int(re.search(r"(\d)/(\d)",genotype).group(1))
                        a2=int(re.search(r"(\d)/(\d)",genotype).group(2))
                        alle1=linelist[a1+3].strip()
                        alle2=linelist[a2+3].strip()
                        pedmap[title[9+i]]+=[alle1,alle2]
                    else:
                        pedmap[title[9+i]]+=['0','0']
            
            elif software.upper()=="SAMTOOLS":
                pass
                for i in range(total_individ):
                    Sample=linelist[i+9]
                    genotype = re.search(r"([^:]+):([^:]+):([^:]+)",Sample.strip()).group(1)
                    pl = re.search(r"([^:]+):([^:]+):([^:]+)",Sample.strip()).group(2)
                    if pl!="0,0,0":
                        a1=int(re.search(r"(\d)/(\d)",genotype).group(1))
                        a2=int(re.search(r"(\d)/(\d)",genotype).group(2))
                        alle1=linelist[a1+3].strip()
                        alle2=linelist[a2+3].strip()
                        pedmap[title[9+i]]+=[alle1,alle2]
                    else:
                        pedmap[title[9+i]]+=['0','0']
                    
        for elem in positionlist:
            print(elem[0],elem[1],elem[2],elem[3],sep='\t',file=mapfile)
        i=1
        for name in pedmap.keys():
            print(i,name,"0","0","1","1","\t".join(pedmap[name]),sep='\t',file=pedfile)
            i+=1       
        mapfile.close()
        pedfile.close()   
        vcffile.close()
    def getVcfListByChrom(self, vcfFileName, chrom,dilute=1,posUniq=True,considerINDEL=False):
        """
            return a list that contain all vcf record of a chrom
        """
        print(chrom)
        VcfList_A_Chrom = []
        if chrom not in self.chromOrder:
            return []
            print(chrom + "didn't find in " + vcfFileName)
        i=self.chromOrder.index(chrom.strip())
        if dilute!=1:
            VcfRecRandomSelectIdxlist=random.sample([j for j in range(self.NumOfRecbychromOrder[i])],int(dilute*self.NumOfRecbychromOrder[i]))
            VcfRecRandomSelectIdxlist.sort()
        elif self.NumOfRecbychromOrder[i]<1000:
            dilute=1
        vcfFile = open(vcfFileName, 'r')
        print("getVcfListByChrom", self.VcfIndexMap[chrom], chrom,int(dilute*self.NumOfRecbychromOrder[i]),self.NumOfRecbychromOrder[i])            
        vcfFile.seek(self.VcfIndexMap[chrom][0])
        line = vcfFile.readline().strip()
        i = 1

        while line and (re.split(r'\s+', line))[0] == chrom:
            if dilute!=1 and len(VcfRecRandomSelectIdxlist)==0:
                break
            elif dilute!=1 and i != VcfRecRandomSelectIdxlist[0]:
                line = vcfFile.readline();i+=1
                continue
            elif dilute!=1 and i == VcfRecRandomSelectIdxlist[0]:
                VcfRecRandomSelectIdxlist.pop(0)
                
                    
            linelist = re.split(r'\s+', line.strip())
            samples=linelist[9:len(linelist)]
            chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            if considerINDEL and len(REF)>1 and len(ALT)>1:
                continue
            INFO = linelist[7]
            FORMAT=linelist[8]
            line = vcfFile.readline();i+=1
            if posUniq and VcfList_A_Chrom and pos==VcfList_A_Chrom[0]:
                print("VCFutil unique the vcf pos",line,VcfList_A_Chrom[-1])
                continue
            VcfList_A_Chrom.append((pos, REF, ALT, INFO,FORMAT,samples))
        vcfFile.close()
        print(chrom,len(VcfList_A_Chrom))    
        return copy.deepcopy(VcfList_A_Chrom)
        

    def getVcfMap(self, vcfFileName):
        """
        this func is from bio\test\posAroundGene\func.py ,and did some improvement,that is add  INFO = collist[7],and add INFO into
        read the vcffile into a map which keys are chrom,values are a list of tuple
        {chrNo:[(pos,REF,ALT,INFO),(pos,REF,ALT,INFO),,,,,],chrNo:[],,,,,,},the order of the tuples in the list,is according pos,
        you we can search a record by  binary chop search
        no matter self.VcfMap_AllChrom has or has not value,the value will be clean
        """
        vcfMap = {}
        vcfFile = open(vcfFileName, 'r')
        
        line = vcfFile.readline()
        while re.search(r'^##', line) != None:
    #        print(line)
            line = vcfFile.readline()
        if re.search(r'^#', line) != None:
            lineslist = vcfFile.readlines()
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'\n" + line)
            exit(-1)
    #    print("pass")
        currentLine = 0
        totalRecs = len(lineslist)
        while currentLine != totalRecs:
    #        print(currentLine)
            collist = re.split(r'\s+', lineslist[currentLine])
            samples=collist[9:len(collist)]
            chrom = collist[0].strip()
            pos = int(collist[1].strip())
            REF = collist[3].strip()
            ALT = collist[4].strip()
            INFO = collist[7]
            FORMAT = collist[8]
            if chrom in vcfMap:
                vcfMap[chrom].append((pos, REF, ALT, INFO,FORMAT,samples))
            else:
                vcfMap[chrom] = [(pos, REF, ALT, INFO,FORMAT,samples)]
            currentLine += 1
        vcfFile.close()
        self.VcfMap_AllChrom = vcfMap
#         for line in self.VcfMap["scaffold8"]:
#             print(line,file =open("vcfMapdata.txt",'a'))
