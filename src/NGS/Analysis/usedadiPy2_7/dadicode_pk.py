import dadi,numpy,re
from numpy import array
import pylab
# fsdata=dadi.Spectrum.from_file("peking_mallard_dilutetodensity0.02.fs")
from optparse import OptionParser
parser = OptionParser()
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
# fsdata=dadi.Spectrum.from_file(options.fsfile)
dd=dadi.Misc.make_data_dict(options.fsfile)
fsdata=dadi.Spectrum.from_data_dict(dd,pop_ids=['mallard14',"pekingduck27"],polarized=True,projections=[26,26])
def split_mig_1(params,ns,pts):
    nuA,nuM,nuB,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_1_IM(params,ns,pts):
    nuA,s,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=nuA*s
    nuB=nuA*(1-s)
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def ancestral_cycle_split_mig_1_IM(params,ns,pts):
    nuA0,nuA1,nuA2,s,TA1,TA2,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    nuA_exp= lambda t: nuA0*(nuA1/nuA0)**(t/TA1)
    #expend
    phi=dadi.Integration.one_pop(phi,xx,TA1,nu=nuA_exp)
    nuA_dec= lambda t: nuA1*(nuA2/nuA1)**(t/TA2)
    phi=dadi.Integration.one_pop(phi,xx,TA2,nu=nuA_dec)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=nuA2*s
    nuB=nuA2*(1-s)
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def domestic_Instantaneous(params,ns,pts):
	nuM,nuP,TA,TB,mM_P,mP_M=params
	xx=dadi.Numerics.default_grid(pts)
	phi=dadi.PhiManip.phi_1D(xx)
	phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuM)
	
	phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi=dadi.Integration.two_pops(phi,xx,TB,nu1=nuM,nu2=nuP,m12=mM_P,m21=mP_M)
	fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs
def domestic_expgrowth(params,ns,pts):
	nuA,nuM,nuP,nuP0,TA,Tsplit,m12,m21=params
	xx=dadi.Numerics.default_grid(pts)
	phi=dadi.PhiManip.phi_1D(xx)
# 	nuM_func= lambda t: nuM0*(nuM/nuM0)**(t/Tsplit)
	phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
	phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
	
	nuP_func= lambda t: nuP0*(nuP/nuP0)**(t/Tsplit)
	phi=dadi.Integration.two_pops(phi,xx,Tsplit,nuM,nuP_func,m12=m12,m21=m21)
	fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs
def domestic_decline_growth(params,ns,pts):
    nuA,nuM0,nuP0,nuM,nuP1,nuPb,nuP,TA,TS,TBM,TBP,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/TS-TBP)
    
    phi=dadi.Integration.two_pops(phi,xx,TS-TBP,nuM0,nuP1_func,m12=m12,m21=m21)
    nuP0=nuP1_func(TS-TBP)
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/TBP)
#     phi=dadi.Integration.two_pops(phi,xx,TBM-TBP,nuM,nuP1_func,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,TBM,nuM,nuP1_func,m12=m12,m21=m21)
    nuP2_func=lambda t: nuPb*(nuP/nuPb)**(t/TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nuM,nuP2_func,m12=m12,m21=m21)
    
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs

def domestic_IM_decline_growth(params,ns,pts):
    nuA,s,nuM,nuP1,nuPb,nuP,TA,Tsplit,TBM,TBP,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
#     nuM_func= lambda t: nuM0*(nuM/nuM0)**(t/Tsplit)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuP0=(1-s)*nuA
    nuM0=s*nuA
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/Tsplit+TBP)
    
    phi=dadi.Integration.two_pops(phi,xx,Tsplit,nuM0,nuP1_func,m12=m12,m21=m21)
    nuP0=nuP1_func(Tsplit)
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/TBP)
#     phi=dadi.Integration.two_pops(phi,xx,TBM-TBP,nuM,nuP1_func,m12=m12,m21=m21)
# M bottle,P continue decreas
    phi=dadi.Integration.two_pops(phi,xx,TBM,nuM,nuP1_func,m12=m12,m21=m21)
    nuP2_func=lambda t: nuPb*(nuP/nuPb)**(t/TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nuM,nuP2_func,m12=m12,m21=m21)
    
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def ancestral_decline_domestic_IM_decline_growth(params,ns,pts):
    nuA0,nuA1,s,nuM1,nuM,nuP1,nuPb,nuP,TA,TS,TBM,TBP,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    nuA_func= lambda t: nuA0*(nuA1/nuA0)**(t/TA-TS)
    #stage1
    phi=dadi.Integration.one_pop(phi,xx,TA-TS,nu=nuA_func)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    #stage 2 
    nuP0=(1-s)*nuA1
    nuM0=s*nuA1
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/(TS-TBM))
    nuM1_func= lambda t:nuM0*(nuM0/nuM1)**(t/(TS-TBP))#especially here differnet time length
    phi=dadi.Integration.two_pops(phi,xx,TS-TBM,nuM1_func,nuP1_func,m12=m12,m21=m21)
    # stage 4 ,M bottle,P continue decreas
    nuP0=nuP1_func(TS-TBM)
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/TBM-TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBM,nuM,nuP1_func,m12=m12,m21=m21)
    #stage 5
    nuP2_func=lambda t: nuPb*(nuP/nuPb)**(t/TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nuM,nuP2_func,m12=m12,m21=m21)
    
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def ancestral_cycle_domestic_IM_decline_growth(params,ns,pts):
    nuA0,nuA1,nuA2,s,nuM1,nuM,nuP1,nuPb,nuP,TA1,TA2,Tsplit,TBM,TBP,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
#     nuA_func= lambda t: nuA0*(nuA1/nuA0)**(t/TA)
    #stage1
    nuA_exp= lambda t: nuA0*(nuA1/nuA0)**(t/TA1)
    #expend
    phi=dadi.Integration.one_pop(phi,xx,TA1,nu=nuA_exp)
    nuA_dec= lambda t: nuA1*(nuA2/nuA1)**(t/TA2)
    phi=dadi.Integration.one_pop(phi,xx,TA2,nu=nuA_dec)
#     phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA_func)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    #stage 2 
    nuP0=(1-s)*nuA1
    nuM0=s*nuA1
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/Tsplit+TBP)
    nuM1_func= lambda t:nuM0*(nuM0/nuM1)**(t/Tsplit+TBP)
    phi=dadi.Integration.two_pops(phi,xx,Tsplit,nuM1_func,nuP1_func,m12=m12,m21=m21)
    # stage 4 ,M bottle,P continue decreas
    nuP0=nuP1_func(Tsplit)
#     phi=dadi.Integration.two_pops(phi,xx,TBM-TBP,nuM,nuP1_func,m12=m12,m21=m21)
    nuP1_func= lambda t: nuP0*(nuP1/nuP0)**(t/TBM+TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBM,nuM,nuP1_func,m12=m12,m21=m21)
    #stage 5
    nuP2_func=lambda t: nuPb*(nuP/nuPb)**(t/TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nuM,nuP2_func,m12=m12,m21=m21)
    
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
ns=fsdata.sample_sizes
pts_1=[40,50,60]
if options.model=="domestic_Instantaneous":
    func=domestic_Instantaneous
elif options.model=="domestic_expgrowth":
    func=domestic_expgrowth
elif options.model=="split_mig_1":
    func=split_mig_1
elif options.model=="domestic_decline_growth":
    func=domestic_decline_growth
elif options.model=="domestic_IM_decline_growth":
    func=domestic_IM_decline_growth
elif options.model=="split_mig_1_IM":
    func=split_mig_1_IM
elif options.model=="ancestral_decline_domestic_IM_decline_growth":
    func=ancestral_decline_domestic_IM_decline_growth
elif options.model=="ancestral_cycle_split_mig_1_IM":
    func=ancestral_cycle_split_mig_1_IM
elif options.model=="ancestral_cycle_domestic_IM_decline_growth":
    func=ancestral_cycle_domestic_IM_decline_growth
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
Nref=theta/(4*9.97e-10*277944.66019938875)
ll_opt=dadi.Inference.ll_multinom(model,fsdata)
print 'Nref',Nref
print theta
for i in range(len(popt)):
    if re.search(r"^T",paramsname[i])!=None or re.search(r"^m",paramsname[i])!=None:
        print "paramname",paramsname[i],params[i],lower_bound[i],upper_bound[i],popt[i],"generation",Nref*popt[i]*2
    else:
        print "paramname",paramsname[i],params[i],lower_bound[i],upper_bound[i],popt[i],"effective pop size",Nref*popt[i]
print 'title:theta,ll_opt',paramsname
print 'Optimized parameters', repr([theta,ll_opt,popt])



print 'title:theta,ll_opt',paramsname
print 'Optimized parameters', repr([theta,ll_opt,popt])
pylab.figure()
dadi.Plotting.plot_single_2d_sfs(fsdata,vmin=1)
pylab.show()
pylab.savefig('fsdata_dom'+options.tag+options.model+'.png', dpi=100)

pylab.figure()
dadi.Plotting.plot_single_2d_sfs(model,vmin=1)
pylab.show()
pylab.savefig('model_dom'+options.tag+options.model+'.png', dpi=100)

pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1)
pylab.show()
pylab.savefig('compare_dom'+options.tag+options.model+'.png', dpi=100)
