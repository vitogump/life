'''
Created on 2014-5-4

@author: liurui
'''
from sqlalchemy import Column, ForeignKey
from sqlalchemy.engine import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.sql.schema import Sequence
from sqlalchemy.sql.sqltypes import Integer, String, Text

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

Base = declarative_base()
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
