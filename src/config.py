# -*- coding: UTF-8 -*-
CSRF_ENABLED=True
SECRET_KEY='you-will-never-guess'
'''
Created on 2013-9-2

@author: liurui
'''
import sys,inspect,os,configparser,random,string
import platform
if not hasattr(sys.modules[__name__], '__file__'):
    __file__ = inspect.getfile(inspect.currentframe())
    
currentpath=os.path.realpath(__file__)
if 'Windows' in platform.system():
    currentpath[:currentpath.find("life\\src")]+"life\\com\\config.properties"
    cfparser = configparser.ConfigParser()
    cfparser.read(currentpath[:currentpath.find("life\\src")]+"life\\com\\config.properties")
else:
    currentpath[:currentpath.find("life/src")]+"life/com/config.properties"
    cfparser = configparser.ConfigParser()
    cfparser.read(currentpath[:currentpath.find("life/src")]+"life/com/config.properties") 

ip=cfparser.get("mysqldatabase","ip")
scriptdir=cfparser.get("mysqldatabase","scriptdir")
print("load in config",currentpath,ip)#currentpath[:currentpath.find("life/src")]+"life/com/config.properties")
username=cfparser.get("mysqldatabase","username")
password=cfparser.get("mysqldatabase","password")
webdbname=cfparser.get("mysqldatabase","webdbname")
genomeinfodbname=cfparser.get("mysqldatabase","genomeinfodbname")
pekingduckchromtable=cfparser.get("mysqldatabase","pekingduckchromtable")
ghostdbname=cfparser.get("mysqldatabase","ghostdbname")
vcfdbname=cfparser.get("mysqldatabase","vcfdbname")
TranscriptGenetable=cfparser.get("mysqldatabase","TranscriptGenetable")
D2Bduckchromtable=cfparser.get("mysqldatabase","D2Bduckchromtable")
KB743256_1=cfparser.get("mysqldatabase","KB743256_1")
outgroupVCFBAMconfig_beijingref=cfparser.get("mysqldatabase","outgroupVCFBAMconfig_beijingref")
pathtoPython=cfparser.get("mysqldatabase", "pathtoPython")
beijingreffa=cfparser.get("mysqldatabase","beijingreffa")

def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])
class CDSMutation():
    def __init__(self):
        self.super()
        self.CodonTable = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
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
              'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}
def findMutation(self,vcfline):
    pass
#原来的SNPAnalysis程序是需要vcf文件排过序才行的