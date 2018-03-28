'''
Created on 2014-5-4

@author: liurui
'''


import datetime, configparser, os
import platform
import time

import markdown2
import mysql.connector
from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.sql.sqltypes import Date, DateTime

import web.entity as Entity


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
print(currentpath,ip)#currentpath[:currentpath.find("life/src")]+"life/com/config.properties")
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

db_config = {
    'host': ip,
    'user': username,
    'passwd': password,
    'db':'ninglabweb',
    'charset':'utf8'
}

engineweb = create_engine('mysql+mysqlconnector://%s:%s@%s/%s?charset=%s'%(db_config['user'],
                                                         db_config['passwd'],
                                                         db_config['host'],
                                                         db_config['db'],
                                                         db_config['charset']), echo=True,pool_recycle=3600)
enginevar = create_engine('mysql+mysqlconnector://%s:%s@%s/%s?charset=%s'%(db_config['user'],
                                                         db_config['passwd'],
                                                         db_config['host'],
                                                         vcfdbname,
                                                         db_config['charset']), echo=False,pool_recycle=3600)
ISOTIMEFORMAT = '%Y-%m-%d %X'
def getVarSession():
    Session = scoped_session(sessionmaker(autoflush=True,bind=enginevar))
    session = Session()
    return session
def getWebSession():
    Session = scoped_session(sessionmaker(autoflush=True,bind=engineweb))
    session = Session()
    return session
def addArticle(name,catalogue_id):
    print("aaaaaaaaaa")
    session = getWebSession()
    rc=Entity.Article(name,catalogue_id)
    session.add(rc)
    session.commit()

def addJobs2jobstate(scriptslist,foldername,state=0):
    
    session = getWebSession()
    for scriptname in scriptslist:
        sc=Entity.Jobstate(foldername=foldername,scriptname=scriptname,state=state)
        session.add(sc)
        session.commit()
def addJobs2jobs_recoder(scriptslist,foldername,logicalpurpose,state=0):
    
    session = getWebSession()
    for scriptname in scriptslist:
        sc=Entity.Jobs_recoder(foldername=foldername,scriptname=scriptname,logicalpurpose=logicalpurpose,state=state)
        session.add(sc)
        session.commit()
print("this line meaning this py module arc execute when import ")
# def addShell(scriptslist):
# session.add(jb)
# session.add(jb2)
# session.commit()


# 
# for i in l:
#     print(i.title,i.catalogue_id)
#c1=entity.Catalogue("mRNA/miRNA表达分析")
#c2=entity.Catalogue("自然选择和人工选择")
#c3=entity.Catalogue("基因印迹和表观遗传")
#c4=entity.Catalogue("基因定位(GWAS,Linkage, NGS等)")
#c5=entity.Catalogue("GBS相关")
#session.add(c1)
#session.add(c2)
#session.add(c3)
#session.add(c4)
#session.add(c5)
#session.add(rep)
#session.commit()


#article_table = Table('articles',metadata,Column('id',Integer,primary_key=True),Column('title',String(1000)))
#metadata.create_all(engine)
#
#metadata = MetaData(engine)
#users_table = Table('articles', metadata, autoload=True)
#print(users_table.columns)