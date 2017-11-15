'''
Created on 2014-5-4

@author: liurui
'''
# from datetime import datetime

import datetime,re
import time
import platform
from sqlalchemy import Column, ForeignKey
from sqlalchemy.engine import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.sql.schema import Sequence
from sqlalchemy.sql.sqltypes import Integer, String, Text, Date, DateTime

import datetime,configparser,os

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

engine = create_engine('mysql+mysqlconnector://%s:%s@%s/%s?charset=%s'%(db_config['user'],
                                                         db_config['passwd'],
                                                         db_config['host'],
                                                         db_config['db'],
                                                         db_config['charset']), echo=True)

Base = declarative_base()
ISOTIMEFORMAT = '%Y-%m-%d %X'
class Jobstate(Base):
    __tablename__="jobsstate"
    id = Column(Integer,Sequence("jobsstate_id_seq"),primary_key=True)
    
    scriptname=Column(String(1000))
    foldername=Column(String(1000))
    scriptcontent=Column(String(1000))
    outputinfo=Column(Text())
# startdate=Column(DateTime, default=time.strftime(ISOTIMEFORMAT, time.localtime()))
    startdate=Column(DateTime)
    finishdate=Column(DateTime)
    state=Column(Integer)
    def __init__(self,scriptname,foldername,state):
        self.scriptname=scriptname
        self.foldername=foldername
        self.state=state

    def __repr(self):
        return "<Jobstat('%s')>"%(self.scriptname)
class Jobs_recoder(Base):
    __tablename__="jobs_recoder"
    id = Column(Integer,Sequence("jobs_recoder_id_seq"),primary_key=True)
    logicalpurpose=Column(String(200))
    scriptname=Column(String(1000))
    foldername=Column(String(1000))
    inputdata=Column(String(1000))
    outputdata=Column(String(1000))
    scriptcontent=Column(String(1000))
    outputinfo=Column(Text())
# startdate=Column(DateTime, default=time.strftime(ISOTIMEFORMAT, time.localtime()))
    startdate=Column(DateTime)
    finishdate=Column(DateTime)
    state=Column(Integer)
    def __init__(self,scriptname,foldername,logicalpurpose,state):
        if foldername[-1]=="/":
            foldername=foldername[0:-1]
        scriptfile=open(foldername+"/"+scriptname,'r')
        line=scriptfile.readline()
        if re.search(r"^scriptinputdata=(.*)",line)!=None:
            self.inputdata=re.search(r"^scriptinputdata=(.*)",line).group(1)
        else:
            self.inputdata="unknow"
        line=scriptfile.readline()
        if re.search(r"^scriptoutputdata=(.*)",line)!=None:
            self.outputdata=re.search(r"^scriptoutputdata=(.*)",line).group(1)
        else:
            self.outputdata="unknow"
#         scriptcontent_all=scriptfile.read()
#         scriptfile.close()
#         if re.search(r"(.*(\n)*)cmdline=\s*(.*)",scriptcontent_all)!=None:
#             scriptcontent_cmdline=re.search(r"(.*(\n)*)cmdline=\s*(.*)",scriptcontent_all).group(3)
#             self.scriptcontent=scriptcontent_cmdline
        self.scriptname=scriptname
        self.foldername=foldername
        self.state=state
        self.logicalpurpose=logicalpurpose

    def __repr(self):
        return "<jobs_recoder('%s')>"%(self.scriptname)
class Catalogue(Base):
    __tablename__="catalogues"
    id = Column(Integer,Sequence("catalogue_id_seq"),primary_key=True)
    title = Column(String(1000))
    article=relationship("Article", backref="catalogues")
    def __init__(self,title):
        self.title=title
    def __repr(self):
        return "<Catalogue('%s')>"%(self.title)
class Article(Base):
    __tablename__ = "articles"
    id = Column(Integer, Sequence("article_id_seq"), primary_key=True)
    title = Column(String(1000))
    catalogue_id=Column(Integer,ForeignKey("catalogues.id"))
    catalogue=relationship("Catalogue",backref=backref('catalogues',order_by=id))
    replys=relationship("Reply", backref="articles")
    def __init__(self,title,catalogue_id):
        self.title=title
        self.catalogue_id=catalogue_id
    def __repr(self):
        return "<Article('%s')>"%(self.title)
class Reply(Base):
    __tablename__="replys"
    id=Column(Integer,primary_key=True)
    content=Column(Text)
    article_id=Column(Integer,ForeignKey("articles.id"))
    article=relationship("Article",backref=backref('articles', order_by=id))
#    replys_id=Column(Integer,ForeignKey("replys.id"))
#    replys=relationship("Reply",backref="reply",order_by=id)
    def __init__(self,content):
        self.content=content
    def __repr(self):
        return "<Reply('%s')>"%(self.content)
Base.metadata.create_all(engine)
