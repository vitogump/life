# -*- coding: UTF-8 -*-
'''
Created on 2017年11月6日

@author: liurui
'''
import copy,pickle,random,re

"""
第二 试述这段代码的功能
"""
def alinmultPopSnpPos(vcfMaplist,jointmode="i"):
    """input:
    two or more map fomart like this [chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],{chrNo:[]},,,,,,]
    output:
    one map like this {chrNo:[(pos,REF,ALT,(INFO,FORMAT,sample,...),(INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
                                            from pop1                        from pop2
    """
    multipleVcfMap={}
    if len(vcfMaplist)==1 or jointmode=="o" or jointmode=="l":

        for currentChrom in vcfMaplist[0].keys():
            multipleVcfMap[currentChrom]=[]
            for SNPrec in vcfMaplist[0][currentChrom]:
                posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt0",'a'))
                RefInPop1 = SNPrec[1]
                AltInPop1 = SNPrec[2]
                multipleVcfMap[currentChrom].append([posInPop1,RefInPop1,AltInPop1,SNPrec[3:]])
        if len(vcfMaplist)==1:
            return copy.deepcopy(multipleVcfMap)
        vcfMap_obj_idx=0
        for vcfMap in vcfMaplist[1:]:
            vcfMap_obj_idx+=1
            for currentChrom in vcfMap:
                for SNPrec in vcfMap[currentChrom]:
                    posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt"+str(vcfMap_obj_idx),'a'))
                    RefInPop1 = SNPrec[1]
                    AltInPop1 = SNPrec[2]
                    low=0;high=len(multipleVcfMap[currentChrom])-1
                    while low<=high:
                        mid = (low + high)>>1
                        if multipleVcfMap[currentChrom][mid][0]<posInPop1:
                            low=mid+1
                        elif multipleVcfMap[currentChrom][mid][0]>posInPop1:
                            high=mid-1
                        else:
                            if AltInPop1 == multipleVcfMap[currentChrom][mid][2]:#same alt alle
                                fillNoneNum=vcfMap_obj_idx-(len(multipleVcfMap[currentChrom][mid])-3)
                                for i in range(fillNoneNum):
                                    multipleVcfMap[currentChrom][mid].append(None)
                                multipleVcfMap[currentChrom][mid].append(SNPrec[3:])
                            break
                    else:
                        if jointmode=="o":
                            insertelem=[posInPop1,RefInPop1,AltInPop1]
                            for i in range(0,vcfMap_obj_idx):
                                insertelem.append(None)
                            insertelem.append(SNPrec[3:])
                            multipleVcfMap[currentChrom].insert(low,insertelem)
    #list(multipleVcfMap.keys())[0]==currentChrom
    #when a pos only exist in the former several pops,but not exist in the rear several pops,the loop block under are neccessary
        for REC_idx in range(0,len(multipleVcfMap[list(multipleVcfMap.keys())[0]])):
            #
            for i in range(len(vcfMaplist)+3-len(multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx])):
                multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx].append(None)

        return copy.deepcopy(multipleVcfMap)
    
    for currentChrom in vcfMaplist[0].keys():
#             self.FstMapByChrom[currentChrom] = []
        multipleVcfMap[currentChrom]=[]
        for SNPrec in vcfMaplist[0][currentChrom]:
            posInPop1 = SNPrec[0]
            RefInPop1 = SNPrec[1]
            AltInPop1 = SNPrec[2]
            elementToAppend=[posInPop1,RefInPop1,AltInPop1,SNPrec[3:]]
            if len(vcfMaplist)==1:
                multipleVcfMap[currentChrom].append(elementToAppend)
                continue
            for vcfMap_obj_idx in range(1,len(vcfMaplist[:])):
                vcfMap_obj=vcfMaplist[vcfMap_obj_idx]
                if currentChrom not in vcfMap_obj or len(vcfMap_obj[currentChrom])==0:
                    print("alinmultPopSnpPos",currentChrom,"didn't find in vcfMap2")
                    break

                low = 0
                high = len(vcfMap_obj[currentChrom]) - 1
                
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
                    continue
#                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#                 print(dp4.group(0))
                
                while low <= high:
                    mid = (low + high)>>1
                    if vcfMap_obj[currentChrom][mid][0]<posInPop1:
                        low=mid+1
                    elif vcfMap_obj[currentChrom][mid][0]>posInPop1:
                        high=mid-1
                    else:
                        if AltInPop1 == vcfMap_obj[currentChrom][mid][2]:#same alt alle
                            if vcfMap_obj_idx!=len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                            elif vcfMap_obj_idx==len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                                multipleVcfMap[currentChrom].append(elementToAppend)
                        break
                else:
                    if jointmode=="i":
                        #ignore the rec
                        break

    return multipleVcfMap

"""
首先从这里读代码，描述getVcfListByChrom() 函数返回的数据结构，（给定vcf文件例子以及，主程序的代码执行之后getVcfListByChrom()返回什么样的数据结构，写出来）
是否有语法错误
"""
class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfIndexMap = {}
        self.chromOrder = []
        self.vcfFileName=vcffileName
        self.NumOfRecbychromOrder = []
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        self.chromOrder = self.VcfIndexMap["chromOrder"]
        self.NumOfRecbychromOrder = self.VcfIndexMap["NumOfRecbychromOrder"]
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')
        vcfChromIndex = {}
        chromOrder = []
        NumOfRecbychromOrder = []
        line = vcffile.readline()
        vcfChromIndex["header"] = [line]
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
        lastChromend_currentChromstartPostion = lastPosition
#         vcfChromIndex[currentChrom]=(lastPosition,0)
        print(line)
        line = vcffile.readline()
        print("first line:",line)
        i = 0
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != linelist[0]:
                chromOrder.append(currentChrom)
                NumOfRecbychromOrder.append(i);i = 0  # collect the  number of snp recs  of the last chrom
                vcfChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)
                lastChromend_currentChromstartPostion = lastPosition
                currentChrom = linelist[0]
                
#                 vcfChromIndex[currentChrom] = (lastPosition,0)
            
            lastPosition = vcffile.tell()

            line = vcffile.readline()
            i += 1
        else:
            chromOrder.append(currentChrom)
            NumOfRecbychromOrder.append(i - 1)  # collect the  number of snp recs  of the lastest chrom of all chroms
            vcfChromIndex[currentChrom] = (lastChromend_currentChromstartPostion, lastPosition)

        vcfChromIndex.pop("temptodele")
        i = chromOrder.index("temptodele")
        if i != 0:
            print("wrong indexVCF")
            exit(-1)
        b = chromOrder.pop(i)
        a = NumOfRecbychromOrder.pop(i)
#         print(i, a, b)
        vcfChromIndex["chromOrder"] = chromOrder
        
        vcfChromIndex["NumOfRecbychromOrder"] = NumOfRecbychromOrder
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
#     def extractVcfRecByChroms(self,vcfFileName,chromlist,replacechromlist,outfile):
#         vcfFile = open(vcfFileName, 'r')
#         if len(chromlist)!=len(replacechromlist):
#             print("the length of the chromlist and replacechromlist should be the same")
#             return
#         for chrom in chromlist:
#             vcfFile.seek(self.VcfIndexMap[chrom][0])

    def getVcfListByChrom(self, chrom,startpos=1,endpos=9999999999999999999999999999999999999999999999999999,  considerINDELandmultpleallele):
        """
            although dilute and dilutetodensity can exist at the same time,but it not make sense and may final produce a bug.
            return a list that contain all vcf record of a chrom
        """


        VcfList_A_Chrom = []
        if (chrom not in self.chromOrder) or startpos>=endpos:
            return []
            print(chrom + "didn't find in " + self.vcfFileName)
        i = self.chromOrder.index(chrom.strip())

        vcfFile = open(self.vcfFileName, 'r')
                    
        vcfFile.seek(self.VcfIndexMap[chrom][0])
        #find the first line
        filepos=vcfFile.tell()
        line=vcfFile.readline()
        
        while line:
#         for line in vcfFile:

            linelist = re.split(r'\s+', line.strip())
            samples = linelist[9:len(linelist)]
            c_chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            if chrom.strip()!=c_chrom:
                return VcfList_A_Chrom            
            if pos>=startpos : 
                break

            filepos=vcfFile.tell()
            line=vcfFile.readline()

            
        vcfFile.seek(filepos)
        linescontent=vcfFile.read(self.VcfIndexMap[chrom][1]-filepos)
        vcflineslist=re.split(r"\n",linescontent.strip())
#         line = vcfFile.readline().strip()
        recidx = 0
#         print(line)
        for line in vcflineslist:
#         while line and (re.split(r'\s+', line))[0] == chrom:

                
            # the code block blow will be excute only when dilute==1,or    dilute!=1 and recidx == VcfRecRandomSelectIdxlist[0]:
            linelist = re.split(r'\s+', line.strip())
            samples = linelist[9:len(linelist)]
            c_chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            if pos>endpos or chrom!=c_chrom:
                break
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            recidx += 1#line = vcfFile.readline();
            if not considerINDELandmultpleallele and (len(REF) > 1 or len(ALT)) > 1:
                continue
            INFO = linelist[7]
            FORMAT = linelist[8]

            VcfList_A_Chrom.append((pos, REF, ALT, INFO, FORMAT, samples))
        vcfFile.close()


        if VcfList_A_Chrom!=[]:
            print(VcfList_A_Chrom[-1][0])
        return VcfList_A_Chrom
        


if __name__ == '__main__':
    """
   3  从这段代码是否有语法错误，
   4 c_b_innerjoin2=alinmultPopSnpPos({b:bchr1},{c:cchr1},"i") 
    """
    1c=VCF_Data("campbell1.pool.withindel.vcf")
    1b=VCF_Data("beijing27.indvd.withindel.vcf")
    bchr1=1b.getVcfListByChrom("1")
    cchr1=1c.getVcfListByChrom("1")
    c_b_innerjoin2=alinmultPopSnpPos({b:bchr1},{c:cchr1},"i")
    """
    [132, 'T', 'C', ['AC=20;AF=0.370;AN=54;BaseQRankSum=1.322;DP=222;Dels=0.00;FS=4.463;HaplotypeScore=0.0358;InbreedingCoeff=0.1584;MLEAC=19;MLEAF=0.352;MQ=34.03;MQ0=13;MQRankSum=-9.958;QD=9.47;ReadPosRankSum=-0.632;SOR=0.254', 'GT:AD:DP:GQ:PL', ['0/0:9,0:9:21:0,21,280', '0/0:8,0:8:18:0,18,242', '1/1:0,9:9:21:222,21,0', '1/1:0,6:6:6:55,6,0', '0/0:7,0:7:21:0,21,276', '0/1:2,4:6:44:44,0,49', '0/1:3,5:8:92:92,0,101', '1/1:0,7:7:12:115,12,0', '0/1:2,4:6:26:106,0,26', '0/1:5,3:8:11:11,0,182', '0/0:9,0:9:27:0,27,335', '0/1:4,1:5:17:17,0,101', '0/1:4,5:9:38:38,0,137', '0/0:13,0:13:30:0,30,401', '0/0:5,0:5:15:0,15,191', '0/0:6,0:6:15:0,15,186', '0/0:7,0:7:21:0,21,261', '0/0:9,0:9:21:0,21,277', '0/0:12,0:12:30:0,30,395', '0/1:8,7:15:66:66,0,240', '0/0:6,0:6:18:0,18,224', '1/1:0,6:6:9:95,9,0', '1/1:0,5:5:9:97,9,0', '0/1:4,2:6:14:14,0,135', '0/0:12,0:12:30:0,30,404', '0/1:5,10:15:99:178,0,111', '0/1:2,6:8:63:107,0,63']], ['AC=1;AF=0.500;AN=2;BaseQRankSum=-1.034;DP=11;Dels=0.00;FS=0.000;HaplotypeScore=0.9947;MLEAC=1;MLEAF=0.500;MQ=25.69;MQ0=0;MQRankSum=-1.398;QD=7.89;ReadPosRankSum=-1.398;SOR=1.270', 'GT:AD:DP:GQ:PL', ['0/1:4,7:11:99:115,0,114']]]

    """
    c_b_innerjoin2["1"][0][4][2][0:5]=['indvd1', 'indvd2', 'indvd3', 'indvd4', 'indvd5']
    """
                请问执行完这一步，之后 如果后续程序需要继续使用bchr1 是否有bug
                
                如果 第266行代码使用 列表 即 VcfList_A_Chrom.append([pos, REF, ALT, INFO, FORMAT, samples]) 如果后续程序需要继续使用bchr1 是否有bug
    
                如果后续程序还将继续使用bchr1
    """