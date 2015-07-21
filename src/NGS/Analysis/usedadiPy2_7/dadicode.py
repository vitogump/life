import dadi,numpy,sys,pylab
from numpy import array
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--fsfile", dest="fsfile",help="fs file name ")
parser.add_option("-m", "--model", dest="model",help="1,model1 2,model2 ....")
parser.add_option("-p","--parameters",dest="parameters",action="append",nargs=4,help="""parametername initvalue lower upper
                                                                                                                red   blue""")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
dd=dadi.Misc.make_data_dict(options.fsfile)
fsdata=dadi.Spectrum.from_data_dict(dd,pop_ids=['mallard14',"spotbilled13"],polarized=True,projections=[26,26])
# fsdata=dadi.Spectrum.from_file(options.fsfile)
for parametername,initvalue,lower,upper in options.parameters:
	if "nuM"==parametername.strip():
		nuM=float(initvalue)
		nuM_u=float(upper)
		nuM_l=float(lower)
	elif "nuB"==parametername.strip():
		nuB=float(initvalue)
		nuB_u=float(upper)
		nuB_l=float(lower)
	elif "nuM0"==parametername.strip():
		nuM0=float(initvalue)
		nuM0_u=float(upper)
		nuM0_l=float(lower)		
	elif "nuB0"==parametername.strip():
		nuB0=float(initvalue)
		nuB0_u=float(upper)
		nuB0_l=float(lower)
	elif "Tsplit"==parametername.strip():
		Tsplit=float(initvalue)
		Tsplit_u=float(upper)
		Tsplit_l=float(lower)
	elif "m"==parametername.strip():
		m=float(initvalue)
		m_u=float(upper)
		m_l=float(lower)
def split_mig(params,ns,pts):
	nuM,nuB,nuM0,nuB0,Tsplit,m=params
	xx=dadi.Numerics.default_grid(pts)
	phi=dadi.PhiManip.phi_1D(xx)
	
	phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
	nuM_func= lambda t: nuM0*(nuM/nuM0)**(t/Tsplit)
	nuB_func= lambda t: nuB0*(nuB/nuB0)**(t/Tsplit)
	phi=dadi.Integration.two_pops(phi,xx,Tsplit,nuM_func,nuB_func,m12=m,m21=m)
	fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs
ns=fsdata.sample_sizes
print 'ns:',repr(ns)
pts_1=[40,50,60]
if options.m=='2':
	func=split_mig
params=array([nuM,nuB,nuM0,nuB0,Tsplit,m])
upper_bound=[nuM_u,nuB_u,nuM0_u,nuB0_u,Tsplit_u,m_u]
lower_bound=[nuM_l,nuB_l,nuM0_l,nuB0_l,Tsplit_l,m_l]
print 'upper_bound',repr(upper_bound)
print 'lower_bound',repr(lower_bound)
print 'params',repr(params)
func_ex=dadi.Numerics.make_extrap_func(func)
p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
popt=dadi.Inference.optimize_log(p0,fsdata,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
model=func_ex(popt,ns,pts_1)
theta=dadi.Inference.optimal_sfs_scaling(model,fsdata)
ll_opt=dadi.Inference.ll_multinom(model,fsdata)

print 'Optimized parameters,theta,ll_opt,popt', repr([theta,ll_opt,popt])
pylab.figure()
dadi.Plotting.plot_single_2d_sfs(fsdata,vmin=1)
pylab.show()
pylab.savefig('fsdata.png', dpi=100)

pylab.figure()
dadi.Plotting.plot_single_2d_sfs(model,vmin=1)
pylab.show()
pylab.savefig('model.png', dpi=100)

pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1)
pylab.show()
pylab.savefig('compare.png', dpi=100)