# -*- coding: UTF-8 -*-
'''
Created on 2017年11月6日
开卷测试，可以查文档和使用手册，测试运行代码。题目从226行开始，一共7道题。
答案写在文本文件或word中，第七题可以回家做，第二天上班前交
@author: liurui
'''
import copy,pickle,re

def func1(vcfMaplist,jointmode="b"):

    multipleVcfMap={}
    if len(vcfMaplist)==1 or jointmode=="a" :

        for currentChrom in vcfMaplist[0].keys():
            multipleVcfMap[currentChrom]=[]
            for SNPrec in vcfMaplist[0][currentChrom]:
                posInPop1 = SNPrec[0]
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
                    posInPop1 = SNPrec[0]
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
                        if jointmode=="a":
                            insertelem=[posInPop1,RefInPop1,AltInPop1]
                            for i in range(0,vcfMap_obj_idx):
                                insertelem.append(None)
                            insertelem.append(SNPrec[3:])
                            multipleVcfMap[currentChrom].insert(low,insertelem)
         #此处              for vcfMap in vcfMaplist[1:]: 循环结束      
        for REC_idx in range(0,len(multipleVcfMap[list(multipleVcfMap.keys())[0]])):
            for i in range(len(vcfMaplist)+3-len(multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx])):
                multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx].append(None)

        return copy.deepcopy(multipleVcfMap)
    # 此处if len(vcfMaplist)==1 or jointmode=="a" : 语句结束
    for currentChrom in vcfMaplist[0].keys():
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
                    print("func1",currentChrom,"didn't find in vcfMap2")
                    break
                low = 0
                high = len(vcfMap_obj[currentChrom]) - 1
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:# multiple allels
                    continue
                
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
                    if jointmode=="b":
                        break
    return multipleVcfMap

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
        vcfChromIndex["chromOrder"] = chromOrder
        
        vcfChromIndex["NumOfRecbychromOrder"] = NumOfRecbychromOrder
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()


    def getVcfListByChrom(self, chrom,startpos=1,endpos=9999999999999999999999999999999999,  considerINDEL_Mult):
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

        recidx = 0
        for line in vcflineslist:
            linelist = re.split(r'\s+', line.strip())
            samples = linelist[9:len(linelist)]
            c_chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            if pos>endpos or chrom!=c_chrom:
                break
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            recidx += 1#line = vcfFile.readline();
            if not considerINDEL_Mult and (len(REF) > 1 or len(ALT)) > 1:
                continue
            INFO = linelist[7]
            FORMAT = linelist[8]
            VcfList_A_Chrom.append((pos, REF, ALT, INFO, FORMAT, samples))
        vcfFile.close()

        if VcfList_A_Chrom!=[]:
            print(VcfList_A_Chrom[-1][0])
        return VcfList_A_Chrom

if __name__ == '__main__':#此句语法没有问题
    #主代码
    1c=VCF_Data("campbell1.pool.withindel.exam.vcf")
    1b=VCF_Data("beijing27.indvd.withindel.exam.vcf")
    bchr1=1b.getVcfListByChrom("KB742833.1")
    cchr1=1c.getVcfListByChrom("KB742833.1")
    c_b_innerjoin2=func1({"KB742833.1":bchr1},{"KB742833.1":cchr1},"b")
    c_b_innerjoin2["1"][0][4][2][0:5]=['indvd1', 'indvd2', 'indvd3', 'indvd4', 'indvd5']
    """
    1。   描述调用getVcfListByChrom() 之后，函数返回的数据的结构是什么样的，即bchr1的数据结构。（VCF文件的例子见    beijing27.indvd.withindel.exam.vcf 和campbell1.pool.withindel.exam.vcf  ）并检查类及方法的的定义及使用代码，以及主代码中是否有语法错误
    2. 为getVcfListByChrom()函数增加一个参数用来接收：将SNP密度随机稀释xxx（0到1之间的一个数）倍，请修改函数代码，使之能够完成位点数的稀释，并返回稀释后的数据。
    3. 试述func1函数的功能
    4。 已知 c_b_innerjoin2["1"][0] 的数据内容及结构为：[132, 'T', 'C', ['AC=20;AF=0.370;AN=54;BaseQRankSum=1.322;DP=222;Dels=0.00;FS=4.463;HaplotypeScore=0.0358;InbreedingCoeff=0.1584;MLEAC=19;MLEAF=0.352;MQ=34.03;MQ0=13;MQRankSum=-9.958;QD=9.47;ReadPosRankSum=-0.632;SOR=0.254', 'GT:AD:DP:GQ:PL', ['0/0:9,0:9:21:0,21,280', '0/0:8,0:8:18:0,18,242', '1/1:0,9:9:21:222,21,0', '1/1:0,6:6:6:55,6,0', '0/0:7,0:7:21:0,21,276', '0/1:2,4:6:44:44,0,49', '0/1:3,5:8:92:92,0,101', '1/1:0,7:7:12:115,12,0', '0/1:2,4:6:26:106,0,26', '0/1:5,3:8:11:11,0,182', '0/0:9,0:9:27:0,27,335', '0/1:4,1:5:17:17,0,101', '0/1:4,5:9:38:38,0,137', '0/0:13,0:13:30:0,30,401', '0/0:5,0:5:15:0,15,191', '0/0:6,0:6:15:0,15,186', '0/0:7,0:7:21:0,21,261', '0/0:9,0:9:21:0,21,277', '0/0:12,0:12:30:0,30,395', '0/1:8,7:15:66:66,0,240', '0/0:6,0:6:18:0,18,224', '1/1:0,6:6:9:95,9,0', '1/1:0,5:5:9:97,9,0', '0/1:4,2:6:14:14,0,135', '0/0:12,0:12:30:0,30,404', '0/1:5,10:15:99:178,0,111', '0/1:2,6:8:63:107,0,63']], ['AC=1;AF=0.500;AN=2;BaseQRankSum=-1.034;DP=11;Dels=0.00;FS=0.000;HaplotypeScore=0.9947;MLEAC=1;MLEAF=0.500;MQ=25.69;MQ0=0;MQRankSum=-1.398;QD=7.89;ReadPosRankSum=-1.398;SOR=1.270', 'GT:AD:DP:GQ:PL', ['0/1:4,7:11:99:115,0,114']]]

                      如果后续程序需要继续读取 bchr1 和cchr1 来对两个品种的vcf数据的 KB742833.1 号 scaffold的数据进行数据分析， 是否存在bug ?(#或者什么样的操作可能导致bug产生)
                
    5. 如果 第222行代码使用 列表 即 VcfList_A_Chrom.append([pos, REF, ALT, INFO, FORMAT, samples]) 同样执行完本段主程序的代码，后续程序需要继续使用bchr1 是否存在bug或者什么样的操作可能导致bug产生，如何解决
    
    """

    """
    6 下列说法正确的是 （可多选）：
    A
    >>> def func(a, b, c=3, d=4): print(a, b, c, d)
    ...
    >>> func(1, *(5,6))
            将输出：1 5 6 4
    B
    >>> cID="chr1"
    >>> snp=("sdf","sfd","1")
    >>> zFst=0.234123435
    >>> pos=234435
    >>> print(cID + "\t" + str(pos) + "\t"  + '%.12f'%(zFst),cID+"sdf",*(snp[1:]+tuple(["nogene"])), sep="\t")
            语法没有错误
    C
    >>> l=[1,2,3,4,5,6,7,8,9]
    >>> while l: # One way to code counter loops
    >>>    a=l.pop(0)
    >>>    print(a, end=' ')
    >>> else:
    >>>    print("finished")
    
            将输出
    0 1 2 3 4 5 6 7 8 9
    
    D
    >>> origin = [1, 2, [3, 4]]
    >>> cop1 = origin
    >>> origin[2][0] = "hey!" 
              此时 origin 的数据为
    [1, 2, ['hey!', 4]]
     cop1 的数据为
    [1, 2, [3, 4]]
    
    E
            执行下述代码
    l=[e for e in range(7)]
    for i in l:
        print(i,end=" ")
        a=l.pop(0)
               最终结果该程序将输出 1 2 3 4 5 6, 最终l==[],l为空列表   
     
    """
    """
    7. 对于一串整数序列，有正数有负数，写一程序找出序列中一段连续的子序列，使得该子序列数字之和最大。（IT人员增加思考附加题：能否设计出时间复杂度为O(n)的算法实现该功能，请写出代码）
                比如 2 -1 2 -2  8  -9  10  -3  -5  1  8  9  3  -8  4  6  -9
    """