# -*- coding: UTF-8 -*-
'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from flask_wtf import FlaskForm
from wtforms import StringField,SubmitField
from wtforms.validators import Required

class ParaForm(FlaskForm):
    projectpath=StringField("请输入项目路径：",validators=[Required()])
    datadepth=StringField("数据所在在层级：",validators=[Required()])
#     sftWname=StringField("请输入所在层级：",validators=[Required()])
    submit = SubmitField('提交')