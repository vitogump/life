# -*- coding: UTF-8 -*-
'''
Created on 2017年10月19日

@author: liurui
'''

from flask_wtf import FlaskForm
from flask_wtf.form import Form
from wtforms import FieldList, StringField, SubmitField, SelectField,HiddenField,TextAreaField
from wtforms.fields.core import FormField
from wtforms.validators import Required, DataRequired
import itertools
from web import web
from wtforms.utils import unset_value
from wtforms import validators as wtf_validators


class FileUploadForm(Form):
    pass

class StudentForm(Form):
    student_id = StringField('Student ID', validators = [DataRequired()])
    student_name = StringField('Student Name', validators = [DataRequired()])

class AddClassForm(Form):
    name = StringField('classname', validators=[DataRequired()])
    day = SelectField('classday', 
                      choices=[(1,"Monday"),(2,"Tuesday"),(3,"Wednesday"),(4,"Thursday"),(5,"Friday")],
                      coerce=int)

    students = FieldList(FormField(StudentForm), min_entries = 5) # show at least 5 blank fields by default

class UserForm(FlaskForm):
# must inherit from wtforms.Form, not flask-WTForms'
# see http://stackoverflow.com/questions/15649027/wtforms-csrf-flask-fieldlist
    first_name = StringField('First Name')
    last_name = StringField('Last Name')

    experience = SelectField('Experience', coerce=int)
class UsersForm(Form):
    users = FieldList(FormField(UserForm), min_entries=2)
_max_nb_entries = 9000
_max_len_per_entry = 360  
_delimiter = "#;_"
class FieldListFromString(FieldList):
    """
    The idea here is to have a FieldList but to store the data in a string format instead of a list
    """
    def process(self, formdata, data=unset_value):
        self.entries = []
        if data is unset_value or not data:
            try:
                data = self.default()
            except TypeError:
                data = self.default
                
        ## Modification from classic FieldList
        if data and 0<len(data):
            data = data.split(_delimiter)
        else:
            data = []
        #
            
        self.object_data = data

        if formdata:
            indices = sorted(set(self._extract_indices(self.name, formdata)))
            if self.max_entries:
                indices = indices[:self.max_entries]

            idata = iter(data)
            for index in indices:
                try:
                    obj_data = next(idata)
                except StopIteration:
                    obj_data = unset_value
                self._add_entry(formdata, obj_data, index=index)
        else:
            for obj_data in data:
                self._add_entry(formdata, obj_data)

        while len(self.entries) < self.min_entries:
            self._add_entry(formdata)


    def populate_obj(self, obj, name):
        values = getattr(obj, name, None)
        try:
            ivalues = iter(values)
        except TypeError:
            ivalues = iter([])

        candidates = itertools.chain(ivalues, itertools.repeat(None))
        _fake = type(str('_fake'), (object, ), {})
        output = []
        for field, data in zip(self.entries, candidates):
            fake_obj = _fake()
            fake_obj.data = data
            field.populate_obj(fake_obj, 'data')
            output.append(fake_obj.data)
        
        ## Modification from classic FieldList
        setattr(obj, name, _delimiter.join(output))
class HiddenFolderTagForm(Form):
    firstpart=HiddenField(id="firstPart")
    secondpart=HiddenField(id="secondPart")
class ParaForm(FlaskForm):
    projectpath=StringField("请输入项目根路径：",default="/USCIRL/project/RD/CTBCTG/liurui")
    datadepth=StringField("数据所在在层级：",default="2")
    batchofinputpath = TextAreaField('一批待处理数据路径',default="""/USCIRL/project/Fastq/lims/A00268/180614_A00268_L4_0092_BH5FY5DSXX/fastq/Project_RTCS/CS00GL18AW0002-1-56_H5FY5DSXX_L4_CS00GL18AW0002-1-56_20180614_P
/USCIRL/project/Fastq/lims/A00262/180621_A00262_L2_0122_BH5FW3DSXX/fastq/Project_RTCS/CS00GL18AW0003-1-55_H5FW3DSXX_L2_CS00GL18AW0003-1-55_20180621_P
""")
#     inputperfix=StringField("输入：")
    collectiondepth=StringField("数据收集层级：",default="2")
    outputpath=StringField("输出路径：            \t ",default="/USCIRL/project/RD/CTBCTG/liurui")  #waiting for changing to be addable with fixed field
    outputperfix=StringField("输出选项 及 后缀：",default="-o bam") 
    outputperfix2=StringField("输出选项 及 后缀：")
    mem=StringField("指定内存（mem）：")
    tagtoFolderlevel=StringField("tag目录层级:",default="1")
    filteredforderlevel=StringField("筛选目录层级:",render_kw={'disabled':'true'},default="1")
#     software=SelectField('选择软件/工具',choices=[("bowtie2","bowtie2"),("bismark","bismark"),("mv","mv"),("cp","cp"),("/Bioinfo/software/install/GTZ-2.0.0/gtz","gtz"),("python ~/software/kosaidtu-norgal-d61342edcdfd/norgal.py","python ~/software/kosaidtu-norgal-d61342edcdfd/norgal.py"),("python ~/life/src/usci/depthDistance.py","python ~/life/src/usci/depthDistance.py"),("selected","othercommand"),("rm","rm")])
#     linuxcommand=StringField("linux命令：",default="/USCIRL/project/RD/CTBCTG/liurui")

    commandParameters = TextAreaField('请输入命令参数（$$$$将被换为输入）：',default=" -p 8 -x /USCIMD/usr/liurui/databases/human/hg19_UCSC_pure --rg-id ID --rg-id PL --rg-id PU --rg-id LB --rg-id SM --rg 'PL:illumina' --rg 'PU:indvd' --rg 'LB:ninglab'  $$$$|samtools view -@ 8 -bS - ")
    messagecomment=TextAreaField('message/comment:',default="prepeocess of cancer early scaning; add tag ")
    ######### this solution is just a temp way , should use FieldList in the further, like filteredforders does#######
    tag1part1=HiddenField(label=None,id='foldertag1')#waiting for changing to be addable
    tag1part2=HiddenField(label=None,id='foldertag2')
    tag2part1=HiddenField(label=None,id='foldertag3')
    tag2part2=HiddenField(label=None,id='foldertag4')
    input1part1=HiddenField(label="输入：",id='pinput1')#waiting for change to be addable
    input1part2=HiddenField(label="输入：",id='pinput2')
    input2part1=HiddenField(label="输入：",id='pinput3')
    input2part2=HiddenField(label="输入：",id='pinput4')
    ##############this solution is just a temp way,#####
    filteredforders = FieldListFromString(StringField('筛选目录名称:',default='mallard',validators=[wtf_validators.Length(min=0, max=_max_len_per_entry)]),
                                  min_entries=1, max_entries=_max_nb_entries)
#     foldername=FieldList(StringField("org"),label='需过滤的文件夹名:',min_entries=2)
    NumOfThreads=StringField("分批运行所有sh命令，每次线程数：",default="1")
    
    submit = SubmitField('提交')