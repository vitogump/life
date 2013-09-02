import re
'''
Created on 2013-7-2

@author: rui
'''
class Caculator():
    def process(self, T):
        pass
    def getResult(self):
        pass
class Caculate_Hp(Caculator):
    def __init__(self):
        super().__init__()
        self.COUNTED=0
        self.CNMI = 0
        self.CNMA = 0
    def process(self, T):
        self.COUNTED+=1
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T[3])
        refdep = int(dp4.group(1)) + int(dp4.group(2))
#         if refdep ==0:##########################################如果ref是0是否考虑该位点
#             return
        altalleledep = int(dp4.group(3)) + int(dp4.group(4))
        if refdep+altalleledep<10:
            return
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
        if self.COUNTED<=3:
            HETEROZY= 'NA'
        self.CNMA = 0
        self.CNMI = 0
        self.COUNTED=0
        return HETEROZY
class Caculate_Fst(Caculator):
    def __init__(self):
        super().__init__()
        self.CNk = 0
        self.CDk = 0
        self.COUNTED=0
    def process(self, T):
        self.COUNTED+=1
        T1 = (T[0], T[1], T[2], T[3])
        T2 = (T[4], T[5], T[6], T[7])
        dp4_1 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T1[3])
        refdep_1 = int(dp4_1.group(1)) + int(dp4_1.group(2))
        altalleledep_1 = int(dp4_1.group(3)) + int(dp4_1.group(4))
        dp4_2 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T2[3])
        refdep_2 = int(dp4_2.group(1)) + int(dp4_2.group(2))
        altalleledep_2 = int(dp4_2.group(3)) + int(dp4_2.group(4))
        if refdep_1==0 and refdep_2==0:
            return  #NOTICT HERE
        if T1[2].strip().upper() == T2[2].strip().upper():
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
        if self.COUNTED<=3:
            Fst='NA'
        self.CDk = 0
        self.CNk = 0
        self.COUNTED=0
        return Fst