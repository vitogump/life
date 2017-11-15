# -*- coding: UTF-8 -*-
'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from flask.ext.wtf import Form
from wtforms import StringField,SubmitField
from wtforms.validators import Required

class ParaForm(Form):
    projectpath=StringField("请输入项目路径：",validators=[Required()])
    datadepth=StringField("请输入所在层级：",validators=[Required()])
#     sftWname=StringField("请输入所在层级：",validators=[Required()])
    submit = SubmitField('Submit')