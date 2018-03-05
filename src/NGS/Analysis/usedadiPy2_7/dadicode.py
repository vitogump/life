from numpy import array, savetxt, loadtxt, std
from optparse import OptionParser
import pickle,math,numpy

import dadi, pylab, re, signal
import matplotlib


matplotlib.use('Agg')


parser = OptionParser()
parser.add_option("-n", "--popname", dest="popname",action="append",nargs=2,help="fs_file_name projection")
parser.add_option("-f", "--fsfile", dest="fsfile",help="fs file name ")
parser.add_option("-l", "--genomelengthwhichsnpfrom", dest="genomelengthwhichsnpfrom",help="fs file name ")
parser.add_option("-b", "--bootstrap", dest="bootstrap",default=False,nargs=2,help="randomstr timetoout")
parser.add_option("-m", "--model", dest="model",help="1,model1 2,model2 ....")
parser.add_option("-p","--parameters",dest="parameters",action="append",nargs=4,help="""parametername initvalue lower upper
                                                                                                                red   blue""")
parser.add_option("-T", "--tag",
                   dest="tag", default="TAG",help="don't print status messages to stdout")
(options, args) = parser.parse_args()
# fsdata=dadi.Spectrum.from_file(options.fsfile)
print float(options.genomelengthwhichsnpfrom) 
dd=dadi.Misc.make_data_dict(options.fsfile)
popnamelist=[]
projectionlist=[]
for popname,projection in options.popname:
    popnamelist.append(popname)
    projectionlist.append(int(projection))
print(options.popname)
fsdata=dadi.Spectrum.from_data_dict(dd,pop_ids=popnamelist,polarized=True,projections=projectionlist)
def split_mig_w_bottleneck(params,ns,pts):
    nuA,nuM,nuP,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuP,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_1_w_bottleneck(params,ns,pts):
    nuA,s1,s2,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=s1*nuA
    nuP=(1-s1)*s2*nuA
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuP,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_bottleneck(params,ns,pts):
    s,Ts,Tb,nuW,nuD,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuD1=s+(nuD-s)*(Ts/(Tb+Ts))
    nuD1_func=lambda t:s+(nuD1-s)*t/Ts
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=(1-s),nu2=nuD1_func,m12=m12,m21=m21)
    nuD2_func=lambda t:nuD1+(nuD-nuD1)*(t/Tb)
    phi=dadi.Integration.two_pops(phi,xx,Tb,nu1=nuW,nu2=nuD2_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_WexpgrothBottleneck_Dlineardecline(params,ns,pts):
    s,Ts,Tb,nuW1,nuW,nuD,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuW1_func=lambda t:(1-s)*(nuW1/(1-s))**(t/Ts)
    nuD1=s+(nuD-s)*(Ts/(Tb+Ts))
    nuD1_func=lambda t:s+(nuD1-s)*(t/Ts)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=nuW1_func,nu2=nuD1_func,m12=m12,m21=m21)
    nuD2_func=lambda t:nuD1+(nuD-nuD1)*(t/Tb)
    phi=dadi.Integration.two_pops(phi,xx,Tb,nu1=nuW,nu2=nuD2_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_WexpgrothBottleneck_Dexpdecline(params,ns,pts):
    s,Ts,Tb,nuW1,nuW,nuD,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuW1_func=lambda t:(1-s)*(nuW1/(1-s))**(t/Ts)
    nuD1=s*(nuD/s)**(Ts/(Tb+Ts))
    nuD1_func=lambda t:s*(nuD1/s)**(t/Ts)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=nuW1_func,nu2=nuD1_func,m12=m12,m21=m21)
    nuD2_func=lambda t:nuD1*(nuD/nuD1)**(t/Tb)
    phi=dadi.Integration.two_pops(phi,xx,Tb,nu1=nuW,nu2=nuD2_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_W_Bottleneckegrowth_Dexpdecline(params,ns,pts):
    s,Ts,Tb,nuWb,nuW,nuD,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
#     nuW1_func=lambda t:(1-s)*(nuW1/(1-s))**(t/Ts)
    nuD1=s*(nuD/s)**(Ts/(Tb+Ts))
    nuD1_func=lambda t:s*(nuD1/s)**(t/Ts)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=(1-s),nu2=nuD1_func)
    nuD2_func=lambda t:nuD1*(nuD/nuD1)**(t/Tb)
    nuW2_func=lambda t:nuWb*(nuW/nuWb)**(t/Tb)
    phi=dadi.Integration.two_pops(phi,xx,Tb,nu1=nuW2_func,nu2=nuD2_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs  
def split_mig(params,ns,pts):
    nuW,nuD,Ts,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=nuW,nu2=nuD,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs

def split_nm(params,ns,pts):
    s,Ts=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=s,nu2=(1-s))
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_m(params,ns,pts):
    s,Ts,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi=dadi.Integration.two_pops(phi,xx,Ts,nu1=s,nu2=(1-s),m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs  
def split_mig_1_IM(params,ns,pts):
    nuA,s,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=s*nuA
    nuB=(1-s)*nuA
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_Afterbottleneck(params,ns,pts):
    nuA,nuAb,s,TA,TAb,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.Integration.one_pop(phi,xx,TAb,nu=nuAb)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM=s*nuA
    nuB=(1-s)*nuA
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def both_bottleneck_aftersplit(params,ns,pts):
    nuA,s,TA,TS,nuM,nuB,TB,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nu1=s*nuA
    nu2=(1-s)*nuA
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nu1,nu2=nu2,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,TB,nu1=nuM,nu2=nuB,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_mig_1_w_bottleneck_split_domDecrease_bottle_increase_wildbottle_IM(params,ns,pts):
    nuA,nuM0,nuP0,nuP1,nuPb,nuP,nuM,TA,TS,TBM,TBP,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
#     nuM0=s*nuA
#     nuP0=(1-s)*nuA
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TS+TBM))
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m12,m21=m21)
    nuP0=nuP_d_func(TS)
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TBM))
    phi=dadi.Integration.two_pops(phi,xx,TBM,nu1=nuM,nu2=nuP_d_func,m12=m12,m21=m21)
    nuP_i_func=lambda t: nuPb*(nuP/nuPb)**(t/(TBP))
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM,nu2=nuP_i_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs

def bottleneckafter_split_mig_1_IM(params,ns,pts):
    nuA,s,nuP,TA,TS,TBP,m12,m21=params
#     if TA<TS:
#         TA=TS+0.0000000001
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
#         print 'adjust',

    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP0,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM0,nu2=nuP)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_2(params,ns,pts):
    s,nu1,nu2,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nu1_func=lambda t: s *(nu1/s)**(t/TS)
    nu2_func=lambda t: (1-s) *(nu2/(1-s))**(t/TS)
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nu1_func,nu2=nu2_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_3(params,ns,pts):
    s1,s2,TS1,TS2,m12,m21,mMB,mBM,mBP,mPB,mMP,mPM=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi=dadi.Integration.two_pops(phi,xx,TS1,nu1=s1,nu2=1-s1,m12=m12,m21=m21)
    phi=dadi.Integration.three_pops(phi,xx,TS2,nu1=s1*s2,nu3=s1*(1-s2),nu2=1-s1,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def increader_split_domDecrease_bottle_increase_wildbottle_IM(params,ns,pts):
    nuA0,nuA1,s,nuM,nuP1,nuPb,nuP,TA,TS,TBM,TBP,m121,m211,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    nuA_func=lambda t: nuA0*(nuA1/nuA0)**(t/TA)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA_func)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA1
    nuP0=(1-s)*nuA1
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TS+TBM))
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m121,m21=m211)
    nuP0=nuP_d_func(TS)
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TBM))
    phi=dadi.Integration.two_pops(phi,xx,TBM,nu1=nuM,nu2=nuP_d_func,m12=m12,m21=m21)
    nuP_i_func=lambda t: nuPb*(nuP/nuPb)**(t/(TBP))
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM,nu2=nuP_i_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def split_domDecrease_bottle_increase_wildbottle_IM(params,ns,pts):
    nuA,s,nuM,nuP1,nuPb,nuP,TA,TS,TBM,TBP,m121,m211,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TS+TBM))
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m121,m21=m211)
    nuP0=nuP_d_func(TS)
    nuP_d_func=lambda t: nuP0*(nuP1/nuP0)**(t/(TBM))
    phi=dadi.Integration.two_pops(phi,xx,TBM,nu1=nuM,nu2=nuP_d_func,m12=m12,m21=m21)
    nuP_i_func=lambda t: nuPb*(nuP/nuPb)**(t/(TBP))
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM,nu2=nuP_i_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def splitdom_domlinerDecrease_IncreaseAfterBottle_wildbottle_mig_1_IM(params,ns,pts):
    nuA,s,nuP1,nuP,TA,Td,Ti,m12,m21=params


    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
    nuP_d_func = lambda t: nuP0 + (nuP1-nuP0)*t/(Td)
#     nuP_g_func= lambda t: nuP0 + (nuP1-nuP0)
    
    phi=dadi.Integration.two_pops(phi,xx,Td,nu1=nuM0,nu2=nuP_d_func,m12=m12,m21=m21)
#     nuP0=nuP_d_func(TS+TBP-TBM)
    nuP_i_func = lambda t: nuP1 * (nuP/nuP1)**(t/Ti)

#     T1=TS+TBP-TBM
#     phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM0,nu2=nuP)
#     nuP_g_func= lambda t: nuP2 + (nuP-nuP2)*t/(TBP)
    phi=dadi.Integration.two_pops(phi,xx,Ti,nu1=nuM0,nu2=nuP_i_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def splitdom_domlinerDecrease_bottledom_mig_1_IM(params,ns,pts):
    nuA,s,nuP1,nuP,TA,TS,TBP,m12,m21=params


    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
    nuP_d_func = lambda t: nuP0 + (nuP1-nuP0)*t/(TS)
    
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m12,m21=m21)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM0,nu2=nuP,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def splitdom_domlinerDecrease_bottledom_expIncrease_mig_1_IM(params,ns,pts):
    nuA,s,nuP1,nuP,TA,TS,TBP,m12,m21=params


    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
    nuP_d_func = lambda t: nuP0 + (nuP1-nuP0)*t/(TS)
    
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m12,m21=m21)
    nuP_i_func = lambda t: nuP1 * (nuP/nuP1)**(t/TBP)
    phi=dadi.Integration.two_pops(phi,xx,TBP,nu1=nuM0,nu2=nuP_i_func,m12=m12,m21=m21)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_growth_bottlesplit(params,ns,pts):
    nuA,nu1,nu2,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    nu_func=lambda t: numpy.exp(numpy.log(nuA)*t/TA)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nu_func)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nu1,nu2=nu2,m12=m12,m21=m21)

    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def IM_bottlesplitegrowth(params,ns,pts):
    nuA,nu1b,nu2b,nu1,nu2,TA,TS,m12,m21=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
#     nu_func=lambda t: numpy.exp(numpy.log(nuA)*t/TA)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nu1_func=lambda t:nu1b*(nu1/nu1b)**(t/TS)
    nu2_func=lambda t:nu2b*(nu2/nu2b)**(t/TS)
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nu1_func,nu2=nu2_func,m12=m12,m21=m21)

    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    return fs
def splitdom_splitwild_3d(params,ns,pts):
    nuA,s1,s2,TA,TS1,TS2,m12,m21,mMB,mBM,mBP,mPB,mMP,mPM=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)

    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuW=nuA*s1
    nuP=nuA*(1-s1)
    phi=dadi.Integration.two_pops(phi,xx,TS1,nu1=nuW,nu2=nuP,m12=m12,m21=m21)
    phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
    nuM=nuW*s2
    nuB=nuW*(1-s2)
    phi=dadi.Integration.three_pops(phi,xx,TS2,nu1=nuM,nu3=nuB,nu2=nuP,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs
def splitdom_splitwild_3d_domlineDecrease(params,ns,pts):
    nuA,nu2,s1,s2,TA,TS1,TS2,m12,m21,mMB,mBM,mBP,mPB,mMP,mPM=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)

    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuW=nuA*s1
    nuP0=nuA*(1-s1)
    nuP_d_func = lambda t: nuP0 + (nu2-nuP0)*t/(TS1+TS2)
    phi=dadi.Integration.two_pops(phi,xx,TS1,nu1=nuW,nu2=nuP_d_func,m12=m12,m21=m21)
    phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
    nuM=nuW*s2
    nuB=nuW*(1-s2)
    nuP0=nuP_d_func(TS1)
    nuP_d_func = lambda t: nuP0 + (nu2-nuP0)*t/(TS2)
    phi=dadi.Integration.three_pops(phi,xx,TS2,nu1=nuM,nu3=nuB,nu2=nuP_d_func,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
    fs=dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs

def splitdom_splitwild_bottledom_3d(params,ns,pts):
    nuA,nuP,s1,s2,TA,TS1,TS2,TBP,m12,m21,mMB,mBM,mBP,mPB,mMP,mPM=params
    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
#     if TS1<TS2:
#         TS1=TS2+0.0000001        
#     if TA<TS1:
#         TA=TS1+0.0000001
    
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuP0=(1-s1)*nuA
    nuMS=s1*nuA
    if TBP>TS2:
        T1=TS1
        T2=TBP-TS2
        T3=TS2
        phi=dadi.Integration.two_pops(phi,xx,T1,nu1=nuMS,nu2=nuP0,m12=m12,m21=m21)
        phi=dadi.Integration.two_pops(phi,xx,T2,nu1=nuMS,nu2=nuP,m12=m12,m21=m21)
        phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
        nuM=nuMS*s2
        nuB=nuMS*(1-s2)
        phi=dadi.Integration.three_pops(phi,xx,T3,nu1=nuM,nu3=nuB,nu2=nuP,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
    else:
        T1=TS1+TBP-TS2
        T2=TS2-TBP
        T3=TBP
        phi=dadi.Integration.two_pops(phi,xx,T1,nu1=nuMS,nu2=nuP0,m12=m12,m21=m21)
        phi=dadi.PhiManip.phi_2D_to_3D_split_1(xx,phi)
        nuM=nuMS*s2
        nuB=nuMS*(1-s2)
        phi=dadi.Integration.three_pops(phi,xx,T2,nu1=nuM,nu3=nuB,nu2=nuP0,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
        phi=dadi.Integration.three_pops(phi,xx,T3,nu1=nuM,nu3=nuB,nu2=nuP,m12=mMP,m21=mPM,m13=mMB,m31=mBM,m23=mPB,m32=mBP)
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
# pts_1=[40,50,60]
pts_1=[50,60,70]
if options.model=="split_mig_1_w_bottleneck":
    func=split_mig_1_w_bottleneck
elif options.model=="split_mig_w_bottleneck":
    func=split_mig_w_bottleneck
elif options.model=="IM_W_Bottleneckegrowth_Dexpdecline":
    func=IM_W_Bottleneckegrowth_Dexpdecline
elif options.model=="IM_bottlesplitegrowth":
    func=IM_bottlesplitegrowth
elif options.model=="split_mig_2":
    func=split_mig_2
elif options.model=="split_mig_1_IM":
    func=split_mig_1_IM
elif options.model=="expand_splitwithbottleneck_mig_1":
    func=expand_splitwithbottleneck_mig_1
elif options.model=="splitdom_splitwild_3d":
    func=splitdom_splitwild_3d
elif options.model=="splitdom_splitwild_bottledom_3d":
    func=splitdom_splitwild_bottledom_3d
elif options.model=="bottleneckafter_split_mig_1_IM":
    func=bottleneckafter_split_mig_1_IM
elif options.model=="splitdom_domlinerDecrease_IncreaseAfterBottle_wildbottle_mig_1_IM":
    func=splitdom_domlinerDecrease_IncreaseAfterBottle_wildbottle_mig_1_IM
elif options.model=="split_mig":
    func=split_mig
elif options.model=="IM_bottleneck":
    func=IM_bottleneck
elif options.model=="IM_WexpgrothBottleneck_Dexpdecline":
    func=IM_WexpgrothBottleneck_Dexpdecline
elif options.model=="IM_WexpgrothBottleneck_Dlineardecline":
    func=IM_WexpgrothBottleneck_Dlineardecline
elif options.model=="IM_growth_bottlesplit":
    func=IM_growth_bottlesplit
elif options.model=="splitdom_domlinerDecrease_bottledom_mig_1_IM":
    func=splitdom_domlinerDecrease_bottledom_mig_1_IM
elif options.model=="splitdom_domlinerDecrease_bottledom_expIncrease_mig_1_IM":
    func=splitdom_domlinerDecrease_bottledom_expIncrease_mig_1_IM
elif options.model=="split_domDecrease_bottle_increase_wildbottle_IM":
    func=split_domDecrease_bottle_increase_wildbottle_IM
elif options.model=="increader_split_domDecrease_bottle_increase_wildbottle_IM":
    func=increader_split_domDecrease_bottle_increase_wildbottle_IM
elif options.model=="splitdom_splitwild_3d_domlineDecrease":
    func=splitdom_splitwild_3d_domlineDecrease
elif options.model=="split_mig_Afterbottleneck":
    func=split_mig_Afterbottleneck
elif options.model=="IM_2":
    func=IM_2
elif options.model=="split_nm":
    func=split_nm
elif options.model=="split_m":
    func=split_m
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
ll_param_MAPlist={}
namestr=""
for name in popnamelist:
    namestr+=name
class TimeOutException(Exception):
    print
    pass
def myHandler_exit(signum, frame):
    print "Now, it's the time,exit"
    exit(-1)
def myHandler_continue(signum, frame):
    print "Now, it's the time,continue"
    raise TimeOutException("time out")
if options.bootstrap!=False:

    #######
#     signal.signal(signal.SIGALRM, myHandler)
#     signal.alarm(int(options.bootstrap[2]))
    
    for name in paramsname:
        ll_param_MAPlist[name]=[]
    ll_param_MAPlist["likelihood"]=[]
    ll_param_MAPlist["theta"]=[]
#     bootstraps=[]
    residualarraylist=[]
    residualhistlist=[]  
#     for cycle in range(int(options.bootstrap[0])):
    ll_param_MAP={}
    for name in paramsname:
        ll_param_MAP[name]=[]
    ll_param_MAP["likelihood"]=[]
    ll_param_MAP["theta"]=[]
    of=open(namestr+options.tag+options.model+options.bootstrap[0]+".parameter","w")
#         print(cycle)
    bootstrap_data=fsdata.sample()
    nsbootstrap=bootstrap_data.sample_sizes
    func_ex=dadi.Numerics.make_extrap_func(func)
    p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
    termitetime=int(options.bootstrap[1])*70
    
    try:
        signal.signal(signal.SIGALRM, myHandler_continue)
        signal.alarm(int(termitetime*2))
        func_ex=dadi.Numerics.make_extrap_func(func)
        popt=dadi.Inference.optimize(p0,bootstrap_data,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
        model=func_ex(popt,nsbootstrap,pts_1)
        signal.alarm(0)
    except :
        print("continue optimize_log_lbfgsb")
        try:
            signal.signal(signal.SIGALRM, myHandler_continue)
            signal.alarm(int(termitetime*1.5))
            func_ex=dadi.Numerics.make_extrap_func(func)
            popt=dadi.Inference.optimize_log_lbfgsb(p0,bootstrap_data,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
            model=func_ex(popt,nsbootstrap,pts_1)
            signal.alarm(0)
        except:
            popt=array([0,0,0])
    
    if (math.isinf(popt[0]) or int(popt[0])==0 or math.isnan(popt[0])) and (math.isinf(popt[-1]) or int(popt[-1])==0 or math.isnan(popt[-1])):
        pts_2=[80,90,100]
        signal.signal(signal.SIGALRM, myHandler_exit)
        signal.alarm(termitetime)
        popt=dadi.Inference.optimize_log_fmin(p0,bootstrap_data,func_ex,pts_2,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
        signal.alarm(0)
        model=func_ex(popt,nsbootstrap,pts_2)
        if (math.isinf(popt[0]) or int(popt[0])==0 or math.isnan(popt[0])) and (math.isinf(popt[-1]) or int(popt[-1])==0 or math.isnan(popt[-1])):
            exit(-1)
#             pts_1=[]
    
    theta=dadi.Inference.optimal_sfs_scaling(model,bootstrap_data)
    ll =dadi.Inference.ll(model*theta,bootstrap_data)
    
    Nref=theta/(4*9.97e-10*float(options.genomelengthwhichsnpfrom))
    for i in range(len(popt)):
        if re.search(r"^T",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[Nref*popt[i]*2,popt[i]]
            print "paramname",paramsname[i],popt[i],"generation",Nref*popt[i]*2
        elif re.search(r"^m",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[popt[i]/(Nref*2),popt[i]]
            
            print "paramname",paramsname[i],popt[i],"migration rate",popt[i]/(Nref*2)
        elif re.search(r"^nu",paramsname[i])!=None or re.search(r"^s",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[popt[i]*Nref,popt[i]]
            print "paramname",paramsname[i],popt[i],"effective pop size",Nref*popt[i]
        else:
            ll_param_MAP[paramsname[i]]=[popt[i],popt[i]]
    ll_param_MAP["likelihood"]=[ll,ll]
    ll_param_MAP["theta"]=[theta,theta]
    btstrap=[ll,theta]+[e for e in popt]
    bof=open(options.fsfile+namestr+options.tag+options.model+options.bootstrap[0]+"btstrap.temp","w")
    for e in btstrap:
        print >>bof,e,
    bof.close()
#     pickle.dump({"ss":btstrap},open(options.fsfile+namestr+options.tag+options.model+options.bootstrap[0]+"btstrap.pickle","wb"),protocol=0)
#     bootstraps.append(btstrap)
    for a in sorted(ll_param_MAP.keys()):
        print >>of,a,ll_param_MAP[a][0],ll_param_MAP[a][1]
    else:
        pass
    of.close()
    if len(popnamelist)==2:
        try:
            fig=pylab.figure(1)
            fig.clear()
            dadi.Plotting.plot_single_2d_sfs(bootstrap_data,vmin=1)
    #         pylab.show()
            fig.savefig('fsdata_split'+namestr+options.tag+options.model+options.bootstrap[0]+'.png', dpi=600)
            fig=pylab.figure(1)
            fig.clear()
            dadi.Plotting.plot_2d_comp_multinom(model,bootstrap_data,vmin=1,residualfilenamepre=options.fsfile+namestr+options.tag+options.model)
            fig.savefig('compare_split'+namestr+options.tag+options.model+options.bootstrap[0]+'.png', dpi=600)
        except:
            print("plot exception")
            exit(-1)
        
        #collection result 
#         residualarray=pickle.load(open(options.fsfile+namestr+options.tag+options.model+"array.pickle","rb"))
# #         u.encoding='latin1'
# #         residualarray=u.load()
#         residualhis=pickle.load(open(options.fsfile+namestr+options.tag+options.model+"hist.pickle","rb"))
# #         u.encoding='latin1'
# #         residualhis=u.load()
#         residualarraylist.append(residualarray)
#         residualhistlist.append(residualhis)
#         inf=open(namestr+options.tag+options.model+options.bootstrap[1]+".parameter","r")
#         for resultline in inf:
#             linelist=re.split(r"\s+",resultline.strip())
#             if len(linelist)>=2:
#                 name=linelist[0]
#                 convert_value=linelist[1]
#                 value=linelist[2]
#                 ll_param_MAPlist[name].append((convert_value,value))
#         inf.close()
#         print ll_param_MAPlist
#     pickle.dump(residualarraylist,open(options.fsfile+namestr+options.tag+options.model+"arraylist.pickle",'wb'))
#     pickle.dump(residualhistlist,open(options.fsfile+namestr+options.tag+options.model+"histlist.pickle",'wb'))
#     fof=open(options.fsfile[:20]+"_"+namestr+options.model+options.tag+".final_parameter","w")

#         sigma_boot = std(bootstraps, axis=0)[1:]
#         fig=pylab.figure(1)
#         fig.clear()
#         fig.hist(bootstraps[:,1], bins=20, normed=True)
#         fig.savefig('hist'+namestr+options.tag+options.model+'.png', dpi=600)
    
#     for i in range(len(ll_param_MAP["likelihood"])):
#     for a in sorted(ll_param_MAPlist.keys()):
#         print >>fof,a+"_convertvalue\t"+a+"_value\t",
#     else:
#         print >>fof
#     for i in range(len(ll_param_MAPlist["likelihood"])):
#         for a in sorted(ll_param_MAPlist.keys()):
#             print a,i
#             print >>fof,ll_param_MAPlist[a][i],
#         else:
#             print >>fof
#     fof.close()
#     if bootstraps is not None:
#         print(bootstraps)
#         bootstraps = array(bootstraps)
#         savetxt(namestr+options.tag+options.model+'2Dboots.npy', bootstraps)
#         bootstraps = loadtxt(namestr+options.tag+options.model+'2Dboots.npy')    
    #     exit()
else:
    func_ex=dadi.Numerics.make_extrap_func(func)
    p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
    popt=dadi.Inference.optimize_log(params,fsdata,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
    model=func_ex(popt,ns,pts_1)
    theta=dadi.Inference.optimal_sfs_scaling(model,fsdata)
    print theta
    ll_opt=dadi.Inference.ll_multinom(model,fsdata)#equal to ll =dadi.Inference.ll(model*theta,fsdata)
    ll=dadi.Inference.ll(model,fsdata)
    Nref=theta/(4*9.97e-10*float(options.genomelengthwhichsnpfrom))
    print 'Nref',Nref
    for i in range(len(popt)):
        if re.search(r"^T",paramsname[i])!=None :
            print "paramname",paramsname[i],popt[i],"generation",Nref*popt[i]*2
        elif re.search(r"^m",paramsname[i])!=None:
            print "paramname",paramsname[i],popt[i],"migration rate",popt[i]/(Nref*2)
        else:
            print "paramname",paramsname[i],popt[i],"effective pop size",Nref*popt[i]
    print 'title:theta,ll,ll_opt',paramsname
    print 'Optimized parameters', repr([theta,ll,ll_opt,popt])
     
    
    
    if len(popnamelist)==2:
        print("print figure")
        fig=pylab.figure(1)
        fig.clear()
        dadi.Plotting.plot_single_2d_sfs(fsdata,vmin=1)
#         pylab.show()
        fig.savefig('fsdata_split'+namestr+options.tag+options.model+'.png', dpi=600)
        
        pylab.figure()
        dadi.Plotting.plot_single_2d_sfs(model,vmin=1)
#         pylab.show()
        pylab.savefig('model_split'+namestr+options.tag+options.model+'.png', dpi=600)
        
        pylab.figure()
        dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1,residualfilenamepre=options.fsfile+namestr+options.tag+options.model)
#         pylab.show()
        pylab.savefig('compare_split'+namestr+options.tag+options.model+'.png', dpi=600)
        pylab.figure()
        dadi.Plotting.plot_2d_comp_Poisson(model,fsdata,vmin=1)
#         pylab.show()
        pylab.savefig('compare_Poisson_split'+namestr+options.tag+options.model+'.png', dpi=100)
    elif len(popnamelist)==3:
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
    
