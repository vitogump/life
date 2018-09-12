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

import os

'''
该文件的执行需要的在Terminal中输入   python setup.py build_ext --inplace ！！！
使用Cpython 编译python文件，关键函数编译成pyd文件（相当于dll）
'''
# 针对多文件情况设置，单文件就只写一个就行
key_funs = ["./pipelinecontrol/Util.py", "test.py", "./NGS/BasicUtil/VCFutil.py","./NGS/BasicUtil/Util.py",
            "./NGS/BasicUtil/DBManager.py", "./NGS/BasicUtil/geneUtil.py", "./NGS/BasicUtil/Caculators.py",
            "./web/views.py", "./web/Service.py","./web/forms.py", "./web/dba.py", "./web/DBA.py","./web/Entity.py"]

setup(
    name="life project", 
    ext_modules = cythonize(key_funs,annotate=True,exclude_failures=True,compiler_directives={"language_level":3})
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
    if fi.__contains__(".so"):
        re_name = fi.split(".")[0] + ".so"
        print(re_name)
        os.rename(fi, re_name)
    elif fi.__contains__(".c") or fi in key_funs:
        os.remove(fi)

print(files)
"gcc -pthread -shared /lustre/home/liurui/software/life/src/build/temp.linux-x86_64-3.6/NGS/BasicUtil/VCFutil.o -o /lustre/home/liurui/software/life/src/NGS/BasicUtil/VCFutil.so"