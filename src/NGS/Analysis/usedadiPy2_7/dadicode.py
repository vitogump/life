import dadi,numpy,pylab,re
from numpy import array
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-n", "--popname", dest="popname",action="append",help="fs file name ")
parser.add_option("-f", "--fsfile", dest="fsfile",help="fs file name ")
parser.add_option("-m", "--model", dest="model",help="1,model1 2,model2 ....")
parser.add_option("-p","--parameters",dest="parameters",action="append",nargs=4,help="""parametername initvalue lower upper
                                                                                                                red   blue""")
parser.add_option("-t", "--tag",
                   dest="tag", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
# fsdata=dadi.Spectrum.from_file(options.fsfile)
dd=dadi.Misc.make_data_dict(options.fsfile)
print(options.popname)
fsdata=dadi.Spectrum.from_data_dict(dd,pop_ids=options.popname,polarized=True,projections=[26]*len(options.popname))
def split_mig_1_w_bottleneck(params,ns,pts):
    nuA,nuM,nuB,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA-TS,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_1_IM(params,ns,pts):
    nuA,s,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA-TS,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=s*nuA
    nuB=(1-s)*nuA
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def bottleneckafter_split_mig_1_IM(params,ns,pts):
    nuA,s,nuP,TA,TS,TBP,m12,m21=params
    if TA<TS:
        TA=TS+0.0000000001
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA-TS,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
#         print 'adjust',

    phi=dadi.Integration.two_pops(phi,xx,TS-TBP,nu1=nuM0,nu2=nuP0,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM0,nu2=nuP,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def bottleneckdom_splitwild_3d(params,ns,pts):
    nuA,nuBP,nuM,nuS,TA,TB,TMS,m12,m21,m13,m31,m23,m32=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    if TA<TB:
        TA=TB+0.0000000001
    phi=dadi.Integration.one_pop(phi,xx,TA-TB,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    if TB<TMS:
        TB=TMS+0.0000000001
    phi=dadi.Integration.two_pops(phi,xx,TB-TMS,nu1=nuA,nu2=nuBP,m12=m12,m21=m21)
    phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
    phi=dadi.Integration.three_pops(phi,xx,TMS,nu1=nuM,nu3=nuS,nu2=nuBP,m12=m12,m21=m21,m13=m13,m31=m31,m23=m23,m32=m32)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs
def splitdom_splitwild_3d(params,ns,pts):
    nuA,s1,s2,TA,TB,TMS,m12,m21,m13,m31,m23,m32=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    if TA<TB:
        TA=TB+0.0000000001
    phi=dadi.Integration.one_pop(phi,xx,TA-TB,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    if TB<TMS:
        TB=TMS+0.0000000001
    nuP=(1-s1)*nuA
    nuMS=s1*nuA
    phi=dadi.Integration.two_pops(phi,xx,TB-TMS,nu1=nuMS,nu2=nuP,m12=m12,m21=m21)
    phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
    nuM=s2*nuMS
    nuS=(1-s2)*nuMS
    phi=dadi.Integration.three_pops(phi,xx,TMS,nu1=nuM,nu3=nuS,nu2=nuP,m12=m12,m21=m21,m13=m13,m31=m31,m23=m23,m32=m32)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs
def expand_splitwithbottleneck_mig_1(params,ns,pts):
    nuA0,nuA1,nuM,nuB,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    nuA_func= lambda t: nuA0*(nuA1/nuA0)**(t/TA)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA_func)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_2(params,ns,pts):
#             nuA, nuM, nuB, TA,  TS,  m
# params=array([2,   1,   1,   0.5  ,0.1  ,1   ])
# upper_bound=[100 , 50,  50,  10   ,2    ,10]
# lower_bound=[1e-3, 1e-3,1e-3,1e-6,1e-6 ,0.01]
    nuA,nuM,nuM0,nuB,nuB0,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM_func= lambda t: nuM0*(nuM/nuM0)**(t/TS)
    nuB_func= lambda t: nuB0*(nuB/nuB0)**(t/TS)
    phi=dadi.Integration.two_pops(phi,xx,TS,nuM_func,nuB_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
ns=fsdata.sample_sizes
pts_1=[40,50,60]
if options.model=="split_mig_1_w_bottleneck":
    func=split_mig_1_w_bottleneck
elif options.model=="split_mig_2":
    func=split_mig_2
elif options.model=="split_mig_1_IM":
    func=split_mig_1_IM
elif options.model=="expand_splitwithbottleneck_mig_1":
    func=expand_splitwithbottleneck_mig_1
elif options.model=="bottleneckdom_splitwild_3d":
    func=bottleneckdom_splitwild_3d
elif options.model=="splitdom_splitwild_3d":
    func=splitdom_splitwild_3d
elif options.model=="bottleneckafter_split_mig_1_IM":
    func=bottleneckafter_split_mig_1_IM
paramslist=[]
upper_boundlist=[]
lower_boundlist=[]
paramsname=[]
for n,v,l,u in options.parameters:
    paramsname.append(n)
    paramslist.append(float(v))
    lower_boundlist.append(float(l))
    upper_boundlist.append(float(u))
params=array(paramslist)
upper_bound=upper_boundlist
lower_bound=lower_boundlist
#             nuA, nuM, nuB, TA,  TS,  m
# params=array([2,   1,   1,   0.5  ,0.1  ,1   ])
# upper_bound=[100 , 50,  50,  10   ,2    ,10]
# lower_bound=[1e-3, 1e-3,1e-3,1e-6,1e-6 ,0.01]
# params=array([10,10,1,0.8,0.1,1])
# upper_bound=[50,50,10,10,2,50]
# lower_bound=[1e-3,1e-3,1e-3,1e-3,1e-6,0.1]
print 'upper_bound',repr(upper_bound)
print 'lower_bound',repr(lower_bound)
print 'params',repr(params)
print 'paramsname',paramsname
func_ex=dadi.Numerics.make_extrap_func(func)
p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
popt=dadi.Inference.optimize_log(p0,fsdata,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
model=func_ex(popt,ns,pts_1)
theta=dadi.Inference.optimal_sfs_scaling(model,fsdata)
print theta
ll_opt=dadi.Inference.ll_multinom(model,fsdata)
Nref=theta/(4*9.97e-10*277919.59937328956)
print 'Nref',Nref
for i in range(len(popt)):
    if re.search(r"^T",paramsname[i])!=None :
        print "paramname",paramsname[i],popt[i],"generation",Nref*popt[i]*2
    elif re.search(r"^m",paramsname[i])!=None:
        print "paramname",paramsname[i],popt[i],"migration rate",popt[i]/(Nref*2)
    else:
        print "paramname",paramsname[i],popt[i],"effective pop size",Nref*popt[i]
print 'title:theta,ll_opt',paramsname
print 'Optimized parameters', repr([theta,ll_opt,popt])
 
namestr=""
for name in options.popname:
    namestr+=name

if len(options.popname)==2:
    pylab.figure()
    dadi.Plotting.plot_single_2d_sfs(fsdata,vmin=1)
    pylab.show()
    pylab.savefig('fsdata_split'+namestr+options.tag+options.model+'.png', dpi=100)
    
    pylab.figure()
    dadi.Plotting.plot_single_2d_sfs(model,vmin=1)
    pylab.show()
    pylab.savefig('model_split'+namestr+options.tag+options.model+'.png', dpi=100)
    
    pylab.figure()
    dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1)
    pylab.show()
    pylab.savefig('compare_split'+namestr+options.tag+options.model+'.png', dpi=100)
elif len(options.popname)==3:
    pylab.figure()
    dadi.Plotting.plot_3d_comp_Poisson(model,fsdata,vmin=1)
    pylab.show()
    pylab.savefig('compare_Poisson_split'+namestr+options.tag+options.model+'.png', dpi=100)
    pylab.figure()
    dadi.Plotting.plot_3d_comp_multinom(model,fsdata,vmin=1)
    pylab.show()
    pylab.savefig('compare_multinom_split'+namestr+options.tag+options.model+'.png', dpi=100)
    pylab.figure()
    dadi.Plotting.plot_3d_spectrum(fsdata,vmin=1)
    pylab.show()
    pylab.savefig('fsdata_split'+namestr+options.tag+options.model+'.png', dpi=100)
    
