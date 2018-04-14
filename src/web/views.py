'''
Created on 2017年10月19日

@author: liurui
'''
from email.utils import unquote
import os, time, Service
import re

from flask import request, jsonify, send_from_directory, abort, render_template
from flask.helpers import url_for
from werkzeug.datastructures import MultiDict
from werkzeug.routing import BaseConverter
from werkzeug.utils import redirect

from src.web.forms import StudentForm
from web import web
from web.forms import ParaForm, UsersForm, AddClassForm, FileUploadForm


class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

@web.route('/mytest',methods=['GET','POST'])        
def testmyform():
    form=ParaForm()
    print("error",form.errors)
    if request.method=='POST'  :
        print("post",form.projectpath.data)
        if form.is_submitted():
            print("here")
            pass
#             return form.projectpath.data+"submitted"
        if form.validate_on_submit():
            print(request.form['projectpath'],"there\n",os.getcwd())
            return form.datadepth.data+"come on"

    return render_template('configsotware.html',form=form)
    

# @web.route('/login',methods=['GET','POST'])
# def login():
#     return render_template('hello.html')
# 
# @web.route('/hello',methods=['POST'])
# def hello():
#     ccc=request.form['textfield']
#     print(ccc,"there\n",os.getcwd())
#     return request.form['textfield']+"come on"

@web.route('/',methods=['GET','POST'])
def configsoftware():
    form=ParaForm()
    if form.validate_on_submit():
        print("configsoftware here",request.form)
        print("hidden value",form.tag1part1.data,form.tag2part1.data)
#         return form.tag1part1.data+form.tag1part2.data+form.tag2part1.data+form.tag2part2.data+form.datadepth.data+form.projectpath.data+"<br />"+form.outputpath.data+form.outputperfix.data+"<br />"+"<br />".join(form.filteredforders.data)
        ppl=re.split(r'[/\\]',form.projectpath.data.strip('/'));scriptdir="/home/liurui/pipeline/";print(ppl)
        for e in ppl:
            orde=0
            for c in e:
                orde+=ord(c)
            scriptdir+=chr(int(orde/len(e)))
        scriptdir+=time.strftime('%Y%m%d', time.localtime()).replace(":","")
        inputList=[];tagList=[]
        for inputpart1,inputpart2 in [(form.input1part1.data,form.input1part2.data),(form.input2part1.data,form.input2part2.data)]:
            if  inputpart2.strip()!="":
                inputList.append(inputpart1+" ${"+inputpart2+"}")
        for tagpart1,tagpart2 in [(form.tag1part1.data,form.tag1part2.data),(form.tag2part1.data,form.tag2part2.data)]:
            if tagpart1.strip()!="":
                tagList.append(tagpart1+"${tag}"+tagpart2)
        outputlist=[form.outputpath.data]+re.split(r"\s+",form.outputperfix.data)#,form.outputperfix2.data]
        tagtofolder=form.tagtoFolderlevel.data if form.tagtoFolderlevel.data.isdigit() else 0
        selectfolderlevel=form.filteredforderlevel.data if form.filteredforderlevel.data.isdigit() else "0"
        scriptsstorediruniq=Service.scriptproduce(form.datadepth.data,form.collectiondepth.data,scriptdir,form.projectpath.data,form.software.data,form.commandParameters.data,inputList,outputlist,lenOfdirtotag=tagtofolder,taglist=tagList,selecteddepth=selectfolderlevel,selecteddirs=list(form.filteredforders.data))
        try:
            t=int(form.NumOfThreads.data)
        except:
            print(form.NumOfThreads.data)
        tt=t if t>1 else 1
#         os.system("cd "+scriptsstorediruniq)
        os.system("chmod +x "+scriptsstorediruniq+"/*.sh")
        os.system("nohup "+Service.aaa.pathtoPython+" ../pipelinecontrol/JobTracker.py -d "+scriptsstorediruniq+" -t "+str(tt)+" -p purposeofthiscommand &")
        currentUstr=scriptsstorediruniq.replace(scriptdir,"")
#         Service.jobminitor(currentUstr)#should store in session
        return redirect("/jobmoinitor"+currentUstr)
#         Service.callsh_updateDB(scriptsstorediruniq,NumOfThread=tt,"purposeofthiscommand")

#         return form.tag1part1.data+form.tag1part2.data+form.tag2part1.data+form.tag2part2.data+scriptdir+form.datadepth.data+form.projectpath.data+"<br />"+form.outputpath.data+form.outputperfix.data+"<br />"+"<br />".join(form.filteredforders.data)
    else:
        print("didn't validate")
        return render_template('commandtemplate.html',form=form)
@web.route('/jobmoinitor/<ustr>', methods=['GET', 'POST'])
@web.route('/jobmoinitor', methods=['GET', 'POST'])
def jobmoinitor(ustr=None):
    html=Service.jobminitor(ustr)
    t="""
    <html>
    <head>
        <title>任务监控</title>
    </head>
    <body>
               <form name="myform" method="post" action="">
           <tr><td>筛选<select name="software" onchange="change(this.value)">
               <option value="1">我最近提交的一批任务</option>
               <option value="2">按日期筛选</option>
               <option value="3">安运行情况筛选</option>
               <option value="4">根据数据筛选</option>
           </select></td></tr><br />
                  </form>
    <br /><p>state:0 任务尚未启动     state:1 任务正在运行        state:2 任务已经完成       state:-1 任务运行失败</p><br />
    %s
        </body>
    </html>
    """
    return t%html
# normally student data is read in from a file uploaded, but for this demo we use dummy data
student_info=[("123","Bob Jones"),("234","Peter Johnson"),("345","Carly Everett"),
              ("456","Josephine Edgewood"),("567","Pat White"),("678","Jesse Black")]
@web.route('/ttt1', methods=['GET', 'POST'])
def addclass():
    fileform = FileUploadForm()
    classform = AddClassForm()

    # Check which 'submit' button was called to validate the correct form
    if 'addclass' in request.form and classform.validate_on_submit():
        # Add class to DB - not relevant for this example.
        return redirect(url_for('addclass'))

    if 'upload' in request.form and fileform.validate_on_submit():
        # get the data file from the post - not relevant for this example.
        # overwrite the classform by populating it with values read from file
        classform = PopulateFormFromFile()
        return render_template('dynamic2.html', classform=classform)

    return render_template('dynamic2.html', fileform=fileform, classform=classform)

def PopulateFormFromFile():
    classform = AddClassForm()

    # normally we would read the file passed in as an argument and pull data out, 
    # but let's just use the dummy data from the top of this file, and some hardcoded values
    classform.name.data = "Super Awesome Class"
    classform.day.data = 4 # Thursday

    # pop off any blank fields already in student info
    while len(classform.students) > 0:
        classform.students.pop_entry()

    for student_id, name in student_info:
    
        studentform = StudentForm()
        studentform.student_id = student_id     # not student_id.data
        studentform.student_name = name
    
        classform.students.append_entry(studentform)
#     for student_id, name in student_info:
        # either of these ways have the same end result.
        #
        # studentform = StudentForm()
        # studentform.student_id.data = student_id
        # studentform.student_name.data = name
        #
        # OR
#         student_data = MultiDict([('student_id',student_id), ('student_name',name)])
#         studentform = StudentForm(student_data)

#         classform.students.append_entry(studentform)

    return classform

@web.route('/ttt', methods=['GET', 'POST'])
def hello_world():
#     return render_template("hello.html")
    form = UsersForm(users=[{}, {}, {}])
    form.users[0].experience.choices=[(1, 'One'), (2, 'Two')]
    form.users[1].experience.choices = [(1, 'Uno'), (2, 'Du')]
    form.users[2].experience.choices = [(0, 'Zero')]
    return render_template("dynamic.html", form=form)
    
@web.route('/downloadfile/:urlpath#.+#')
def send_static(urlpath):
    print("send_static")
    filename=re.search(r'[^/]*$',unquote(urlpath)).group(0)
    path="../../"+re.search(r'^.*/',unquote(urlpath)).group(0)
    print(path,filename,unquote(urlpath))
#    print(urlpath,re.search(r'^.*/',urlpath).group(0),re.search(r'[^/]*$',urlpath).group(0))
    return render_template(path+filename)