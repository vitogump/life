#!python
#cython: language_level=3, boundscheck=False
'''
Created on 2018年9月10日

@author: Dr.liu
'''
# 命令行下：（在项目目录下打开命令行或者shell，该命令只能编译一个文件，编译之后会发现出现三个文件，yourmod.c、yourmod.html、yourmod-win_amd64.pyd，此时将c、html和原py文件删除，将pyd文件命名更改为yourmod就可以）
#cythonize -3 -a -i yourmod.pyx

import pyximport

from distutils.core import setup
from Cython.Build import cythonize
pyximport.install(pyimport=True,language_level =3)

import os,shutil

'''
该文件的执行需要的在Terminal中输入   python setup.py build_ext --inplace ！！！
使用Cpython 编译python文件，关键函数编译成pyd文件（相当于dll）
'''
currdir = os.path.abspath('.') + '\\'
build_dir="build"
build_tmp_dir=build_dir + "/temp"
filter_dir_set = {'templates', 'plugin', 'static', 'orm\\data'}
except_files = {
    __file__,
}
# 针对多文件情况设置，单文件就只写一个就行
#so
toso_list=["./NGS/BasicUtil","./NGS/RUtil"]
#pyc(py)
topyc_list=["./pipelinecontrol","./compass","./NGS/Slave","./NGS/SwissArmyKnife","./NGS/Analysis"]

def filter_file(file_name):
    if file_name.__contains__(currdir):
        file_name = file_name.replace(currdir, '')
    if file_name in except_files:  # 过滤文件
        return True
    file_path = file_name.split("\\")
    if len(file_path) > 1:
        file_dir = ""
        for i in range(len(file_path)-1):
            file_dir = os.path.join(file_dir, file_path[i])
            if file_dir in filter_dir_set:
                return True
    return file_path[0] in filter_dir_set
def getpy(basepath=os.path.abspath('.'), parentpath='', name='',
          copyOther=False, delC=False):
    """
    获取py文件的路径
    :param basepath: 根路径
    :param parentpath: 父路径
    :param name: 文件/夹
    :param copy: 是否copy其他文件
    :return: py文件的迭代器
    """
    fullpath = os.path.join(basepath, parentpath, name)
    print(fullpath,basepath,parentpath,name)
    for fname in os.listdir(fullpath):
        ffile = os.path.join(fullpath, fname)
        if os.path.isdir(ffile) and fname != build_dir and not fname.startswith('.'):
            for f in getpy(basepath, os.path.join(parentpath, name), fname,
                          copyOther, delC):
                yield f
        elif os.path.isfile(ffile):
            ext = os.path.splitext(fname)[1]
            # 删除.c 临时文件
            if ext == ".c":
                if delC:
                    os.remove(ffile)
            elif not filter_file(ffile) and (ext not in ('.pyc', '.pyx')
                                             and ext in ('.py', '.pyx')
                                             and not fname.startswith('__')):
                yield os.path.join(parentpath, name, fname)
            elif copyOther and ext not in ('.pyc', '.pyx'):  # 复制其他文件到./build 目录下
                dstdir = os.path.join(basepath, build_dir, parentpath, name)
                if not os.path.isdir(dstdir):
                    os.makedirs(dstdir)
                shutil.copyfile(ffile, os.path.join(dstdir, fname))
        else:
            print(ffile)
            pass
#so_list=list(getpy(basepath=currdir, parentpath="/home/lrui/life/src/NGS/BasicUtil/"))
key_funs = ["./web/Service.py","./web/forms.py", "./web/dba.py", "./web/DBA.py","./web/views.py","./NGS/Service/Ancestralallele.py"]#["./pipelinecontrol/Util.py", "./NGS/Analysis/Detectsignalacrossgenome_master.py", "./NGS/BasicUtil/VCFutil.py","./NGS/BasicUtil/Util.py",
            #"./NGS/BasicUtil/DBManager.py", "./NGS/BasicUtil/geneUtil.py", "./NGS/BasicUtil/Caculators.py",
            #"./web/views.py", "./web/Service.py","./web/Entity.py"]
for parentpath in toso_list:
    key_funs+=list(getpy(parentpath=parentpath))
print(key_funs)
setup(
    name="life project", 
    ext_modules = cythonize(key_funs,annotate=True,exclude_failures=True,compiler_directives={"language_level":3}),
    script_args=["build_ext", "-b", build_dir, "-t", build_tmp_dir]
    )

'''
1、将编译后的so文件的命名更改成与原py文件一致
2、删除编译后得到的c文件和原py文件
'''
rootDir=os.getcwd()
print("——————", rootDir, "——————")
path="~/"

def list_all_files(rootdir):
#     import os
    _files = []
    list = os.listdir(rootdir) #列出文件夹下所有的目录与文件
    for i in range(0,len(list)):
           path = os.path.join(rootdir,list[i])
           if os.path.isdir(path):
              _files.extend(list_all_files(path))
           if os.path.isfile(path):
              _files.append(path)
    return _files
files = list_all_files(rootDir)

for fi in files:
    if fi.endswith(".so"):
        re_name = fi.split(".")[0] + ".so"
        print("rename",fi,re_name)
        os.rename(fi, re_name)
    elif fi.endswith(".c") or fi in key_funs:
        print("remove:",fi)
        os.remove(fi)
    elif not fi.endswith(".html") and "./"+os.path.relpath(fi, "/home/lrui/life/src") not in key_funs:
        distfile=os.path.join(os.path.abspath('.'), build_dir, os.path.relpath(fi, "/home/lrui"))
        if not os.path.exists(os.path.dirname(distfile)):
            os.makedirs(os.path.dirname(distfile))
        shutil.copy(fi,os.path.join(os.path.abspath('.'), build_dir, os.path.relpath(fi, "/home/lrui")))

print(files)
pyc_file_list=["./NGS/Service/fillancestraltemp_usingBAM.py","./NGS/Service/fillancestraltemp_usingdepth.py","./NGS/Service/makingjointabletemp.py"]
for parentpath in topyc_list:
    pyc_file_list+=list(getpy(parentpath=parentpath))
print("python -m compileall -b ",pyc_file_list)
for pyc_fi in pyc_file_list:
    os.system("python -m compileall -b "+pyc_fi)
    distfile=os.path.join(os.path.abspath('.'), build_dir, os.path.relpath(pyc_fi, "/home/lrui"))
    print("distfile",distfile)
    if not os.path.exists(os.path.dirname(distfile)):
        os.makedirs(os.path.dirname(distfile))
    os.system("mv "+pyc_fi+"c "+distfile)
"python -m compileall -b CaculateKaKsusingMusclePaml.py"#可以生成pyc文件，pyc文件可以直接想py文件一样使用
"gcc -pthread -shared /lustre/home/liurui/software/life/src/build/temp.linux-x86_64-3.6/NGS/BasicUtil/VCFutil.o -o /lustre/home/liurui/software/life/src/NGS/BasicUtil/VCFutil.so"