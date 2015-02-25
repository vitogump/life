# -*- coding: UTF-8 -*-
'''
Created on 2013-6-30

@author: rui
'''


import re, pickle, copy





     
class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfIndexMap = {}
        self.chromOrder=[]
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        self.chromOrder=self.VcfIndexMap["chromOrder"]
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')
        vcfChromIndex = {}
        chromOrder=[]
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
        vcfChromIndex[currentChrom]=(lastPosition,0)
        print(line)
        line = vcffile.readline()
        print(line)
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != linelist[0]:
                vcfChromIndex[currentChrom][1]=lastPosition
                currentChrom = linelist[0]
                chromOrder.append(currentChrom)
                vcfChromIndex[currentChrom] = (lastPosition,0)
            lastPosition = vcffile.tell()
            
    
            line = vcffile.readline()
        vcfChromIndex.pop("temptodele")
        vcfChromIndex["chromOrder"]=chromOrder
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
#     def extractVcfRecByChroms(self,vcfFileName,chromlist,replacechromlist,outfile):
#         vcfFile = open(vcfFileName, 'r')
#         if len(chromlist)!=len(replacechromlist):
#             print("the length of the chromlist and replacechromlist should be the same")
#             return
#         for chrom in chromlist:
#             vcfFile.seek(self.VcfIndexMap[chrom][0])
            
    def Vcf2Ped(self,vcfFileName,outputfileprefix,software,withheader=False):
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
            total_individ=len(self.VcfIndexMap["title"])-9
            print(self.VcfIndexMap["title"],len(self.VcfIndexMap["title"]),total_individ)
            for outName in self.VcfIndexMap["title"][len(self.VcfIndexMap["title"])-total_individ:]:
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
    def getVcfListByChrom(self, vcfFileName, chrom,posUniq=True,considerINDEL=False):
        """
            return a list that contain all vcf record of a chrom
        """
        VcfList_A_Chrom = []
        vcfFile = open(vcfFileName, 'r')
        try:
            print("getVcfListByChrom", self.VcfIndexMap[chrom], chrom)            
            vcfFile.seek(self.VcfIndexMap[chrom][0])
            line = vcfFile.readline().strip()
        except KeyError:
            print(chrom + "didn't find in " + vcfFileName)
            return []
        while line and (re.split(r'\s+', line))[0] == chrom:
            linelist = re.split(r'\s+', line)
            samples=linelist[9:len(linelist)]
            chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            if considerINDEL and len(REF)>1 and len(ALT)>1:
                continue
            INFO = linelist[7]
            FORMAT=linelist[8]
            line = vcfFile.readline()
            if posUniq and VcfList_A_Chrom and pos==VcfList_A_Chrom[0]:
                print("VCFutil unique the vcf pos",line,VcfList_A_Chrom[-1])
                continue
            VcfList_A_Chrom.append((pos, REF, ALT, INFO,FORMAT,samples))
            
        return copy.deepcopy(VcfList_A_Chrom)
        vcfFile.close()

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
