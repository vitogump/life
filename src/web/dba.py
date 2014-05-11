'''
Created on 2014-5-4

@author: liurui
'''
from sqlalchemy import *
from sqlalchemy.orm import *
import mysql.connector
import web.entity as entity

db_config = {
    'host': '10.2.48.96',
    'user': 'root',
    'passwd': '1234567',
    'db':'ninglabweb',
    'charset':'utf8'
}

engine = create_engine('mysql+mysqlconnector://%s:%s@%s/%s?charset=%s'%(db_config['user'],
                                                         db_config['passwd'],
                                                         db_config['host'],
                                                         db_config['db'],
                                                         db_config['charset']), echo=True)
def addArticle(name,catalogue_id):
    print("aaaaaaaaaa")
    Session = sessionmaker(bind=engine)
    session = Session()
    rc=entity.Article(name,catalogue_id)
    session.add(rc)
    session.commit()
Session = sessionmaker(bind=engine)
session = Session()   
l = session.query(entity.Article).all()

for i in l:
    print(i.title,i.catalogue_id)
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