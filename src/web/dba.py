'''
Created on 2014-5-4

@author: liurui
'''

import datetime
import time

import markdown2
import mysql.connector
import src.web.Entity as Entity
from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.sql.sqltypes import Date, DateTime



db_config = {
    'host': '10.2.48.147',
    'user': 'root',
    'passwd': '1234567',
    'db':'ninglabweb',
    'charset':'utf8'
}

engine = create_engine('mysql+mysqlconnector://%s:%s@%s/%s?charset=%s'%(db_config['user'],
                                                         db_config['passwd'],
                                                         db_config['host'],
                                                         db_config['db'],
                                                         db_config['charset']), echo=True,pool_recycle=3600)

ISOTIMEFORMAT = '%Y-%m-%d %X'
def getSession():
    Session = scoped_session(sessionmaker(autoflush=True,bind=engine))
    session = Session()
    return session
def addArticle(name,catalogue_id):
    print("aaaaaaaaaa")
    session = getSession()
    rc=Entity.Article(name,catalogue_id)
    session.add(rc)
    session.commit()

def addJobs2jobstate(scriptslist,foldername,state=0):
    
    session = getSession()
    for scriptname in scriptslist:
        sc=Entity.Jobstate(foldername=foldername,scriptname=scriptname,state=state)
        session.add(sc)
        session.commit()
def addJobs2jobs_recoder(scriptslist,foldername,logicalpurpose,state=0):
    
    session = getSession()
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