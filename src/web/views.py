'''
Created on 2017年10月19日

@author: liurui
'''
from email.utils import unquote
import os, time, Service,re

from functools import wraps  # 装饰器(用于访问控制)
from flask import request, jsonify, send_from_directory, abort, render_template, Response
from flask.helpers import url_for
from src.web.forms import StudentForm
from web import web
from web.forms import *
from werkzeug.datastructures import MultiDict
from werkzeug.routing import BaseConverter
from werkzeug.utils import redirect
from email.policy import default


class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

# 前台登录装饰器(只能登录后才能访问会员中心)
# def user_login_req(f):
#     @wraps(f)
#     def decorated_function(*args, **kwargs):
#         if 'user' not in session:  # if session['user'] is None:
#             return redirect(url_for('login', next=request.url))
#         return f(*args, **kwargs)
# 
#     return decorated_function




@web.route('/',methods=['GET','POST'])
def upload_file():
    if request.method=='POST':
        file=request.files['file']
        if file:
            filename=file.filename.rsplit('.',1)[0]
            filetype=file.filename.rsplit('.',1)[1]
            filename=filename.replace(' ','')
            filename=filename.replace('.','-')
            filename=filename+'.'+filetype
            #解决命名冲突
#             records=Service.getrecords()
#             for record in records:
#                 if filename==record['filename']:
#                     flash('')
#                     return render_template("upload.html")
            file.save(os.path.join(web.config['UPLOAD_FOLDER'],filename))
#             date=time.strftime(,time.localtime())
#             record={"filename":filename,"date":date}
#             cursor.execute()
#             data_base.commit()
            return redirect(url_for('identifywork',filename=filename))
    return render_template('hlsbwelcome.html')
            
@web.route('/video_feed',defaults={'filename':'monkey.avi'})
@web.route('/video_feed/<filename>')
def video_feed(filename):
    path=web.config['UPLOAD_FOLDER']+"/"+filename
    print("i find "+path)
    return Response(Service.gen(path),
                    mimetype='multipart/x-mixed-replace; boundary=frame')            
@web.route('/demoexhibitionmp4')
def demo_exhibition():
    return render_template('video.html')
    return render_template('videomp4.html')
@web.route('/identifwork/<filename>',methods=['GET'])
def identifywork(filename):
    
    return render_template('hlsbwork.html',path=filename)
#     return """
#     <html>
#   <head>
#     <title>Video Streaming Demonstration</title>
#   </head>
#   <body>
#     <h1>Video Streaming Demonstration</h1>
#     <img src="{{ url_for('video_feed') }}">
#   </body>
# </html>
#     """
#     return ''' <!DOCTYPE html> 
#     <html> <body> 
#     <video width="700" height="500" controls="controls"> 
#     <source src="static/video/my_movie.mp4" type="video/mp4" /> 
#     </video> 
#     </body> 
#     </html> 
#     '''


@web.route('/hom',methods=['GET','POST'])
def guide():
    return redirect(url_for('tutorial'))
    return render_template("hlsbwelcome.html")

@web.route('/login',methods=['GET','POST'])
def login():
    form=LoginForm()
    if form.validate_on_submit():
        data=form.data
    return render_template("login.html",form=form)
@web.route('/mytest',methods=['GET','POST'])        
def testmyform():
    form=ParaForm()
    print("error",form.errors)
    if request.method=='POST'  :
        print("post",form.data)
        return "what"

        if form.is_submitted():
            print("here")
            pass
#             return form.projectpath.data+"submitted"
        if form.validate_on_submit():
            #print(request.form['projectpath'],"there\n",os.getcwd())
            return "come on"

    return render_template('entryelem.html',form=form)
@web.route("/analysisdata2",methods=['GET','POST'])
def processData():
    return render_template('index.html')

@web.route('/analysisdata',methods=['GET','POST'])
def configsoftware():
    data=request.get_json(silent=True)
    print(request.args.to_dict())
    print(request.json,data,request.form)
    form=ParaForm()
    if form.validate_on_submit():
        print("configsoftware here",type(request.form),request.form)
        print("hidden value",form.tag1part1.data,form.tag2part1.data,"software:",request.form['software1'],type(request.form['software1']))
        scriptdir=Service.scriptdir
        tagtofolder=int(form.tagtoFolderlevel.data) if form.tagtoFolderlevel.data else 0
        if "Checkbox1" in request.form:
            selectfolderlevel=form.filteredforderlevel.data
            print("request.form Checkbox1 seletecfolderlevel",selectfolderlevel)
        else:
            selectfolderlevel="0"
#         return form.tag1part1.data+form.tag1part2.data+form.tag2part1.data+form.tag2part2.data+form.datadepth.data+form.projectpath.data+"<br />"+form.outputpath.data+form.outputperfix.data+"<br />"+"<br />".join(form.filteredforders.data)
        if form.projectpath.data.strip('/').strip()!="":
            ppl=re.split(r'[\\]',form.projectpath.data.strip('/'))
            print("type form.projectpath.data",type(form.projectpath.data),"ppl",ppl)
            for e in ppl:
                orde=0
                for c in e:
                    orde+=ord(c)
                print(orde/len(e))
                scriptdir+="/"
#             scriptdir+=chr(int(orde/len(e)))+"/"
        scriptdir=scriptdir.rstrip()+"/"+time.strftime('%Y%m%d', time.localtime()).replace(":","")
        inputList=[];tagList=[];outputlist=[form.outputpath.data.rstrip(os.sep)]
        for inputpart1,inputpart2 in [(form.input1part1.data,form.input1part2.data),(form.input2part1.data,form.input2part2.data)]:# this code need to be modified !
            if  inputpart2.strip()!="":
                if inputpart1.endswith("="):
                    inputList.append((inputpart1.strip(),inputpart2))
                else:
                    inputList.append((inputpart1+" ",inputpart2))
        for tagpart1,tagpart2 in [(form.tag1part1.data,form.tag1part2.data),(form.tag2part1.data,form.tag2part2.data)]:
            if tagpart1.strip()!="" and form.tagtoFolderlevel.data.strip()!="":
                tagList.append(tagpart1+"${tag}"+tagpart2)
        for outputpart1 in list(form.outOptSuffixs.data):
            if outputpart1.strip()!="":
                outputlist.append(re.split(r"\s+",outputpart1))#,form.outputperfix2.data]
        batchofInList=re.split(r"\s+",form.batchofinputpath.data.strip());print(type(form.batchofinputpath.data),batchofInList)
        softwareconfig= request.form['MyLinuxCommand'] if request.form['software1'] =="othercommand" else request.form['software1']
        msg=re.sub(r"[\s\W]+","_",form.messagecomment.data.strip())
        scriptsstorediruniq=Service.scriptproduce(form.datadepth.data,form.collectiondepth.data,(scriptdir,msg),form.projectpath.data.rstrip(os.sep),softwareconfig,form.commandParameters.data,inputList,outputlist,batchofInList,lenOfdirtotag=tagtofolder,taglist=tagList,selecteddepth=selectfolderlevel,selecteddirs=list(form.filteredforders.data))
        try:
            t=int(form.NumOfThreads.data)
        except:
            print(form.NumOfThreads.data)
        tt=t if t>1 else 1
#         os.system("cd "+scriptsstorediruniq)
        os.system("chmod +x "+scriptsstorediruniq+"/*.sh")
        qsubMorNot=" " if not form.mem.data.strip() else (" -m "+form.mem.data)
        currentUstr=scriptsstorediruniq.replace(scriptdir,"").strip("/").replace(os.sep,"_")
        if softwareconfig.strip()=="rm" or softwareconfig.strip()=="mv": return ('mv rm  command not execute,please manually execute; redirect("/jobmoinitor/"+currentUstr)'+currentUstr)
        os.system("nohup "+Service.aaa.pathtoPython+" ../pipelinecontrol/JobTracker.py -d "+scriptsstorediruniq+" -t "+str(tt)+qsubMorNot+" -p "+ msg+" &")
        
        print("currentUstr",currentUstr,"return url:","/jobmoinitor/"+currentUstr)
#         Service.jobminitor(currentUstr)#should store in session

        return redirect("/jobmoinitor/"+currentUstr)
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