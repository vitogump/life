from optparse import OptionParser

from numpy import array

import dadi, numpy, pylab, re,signal


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
def splitdom_domlinerDecrease_mig_1_IM(params,ns,pts):
    nuA,s,nuP1,TA,TS,m12,m21=params


    xx=dadi.Numerics.default_grid(pts)
    phi=dadi.PhiManip.phi_1D(xx)
    phi=dadi.Integration.one_pop(phi,xx,TA,nu=nuA)
    phi=dadi.PhiManip.phi_1D_to_2D(xx,phi)
    nuM0=s*nuA
    nuP0=(1-s)*nuA
    nuP_d_func = lambda t: nuP0 + (nuP1-nuP0)*t/(TS)
    
    phi=dadi.Integration.two_pops(phi,xx,TS,nu1=nuM0,nu2=nuP_d_func,m12=m12,m21=m21)

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
# def split_mig_1_w_bottleneck_split_domedcrease_increaseAfterBottle_wild_bottle_split_3d(params,ns,pts):
#     nuA,nuP,nuW0,nuWb,nuM,nuS,nuP0,nuP1,nuPb,nuP,nuS,TA,TS1,TBW,TS2,TBP,m12,m21,mMB,mBM,mBP,mPB,mMP,mPM=params
#     xx=dadi.Numerics.default_grid(pts)
#     phi=dadi.PhiManip.phi_1D(xx)
#     unfinished
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
pts_1=[40,50,60]
if options.model=="split_mig_1_w_bottleneck":
    func=split_mig_1_w_bottleneck
elif options.model=="split_mig_w_bottleneck":
    func=split_mig_w_bottleneck
elif options.model=="split_mig_1_w_bottleneck_split_domDecrease_bottle_increase_wildbottle_IM":
    func=split_mig_1_w_bottleneck_split_domDecrease_bottle_increase_wildbottle_IM
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
elif options.model=="IM_2":
    func=IM_2
elif options.model=="IM_3":
    func=IM_3
elif options.model=="splitdom_domlinerDecrease_mig_1_IM":
    func=splitdom_domlinerDecrease_mig_1_IM
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
elif options.model=="both_bottleneck_aftersplit":
    func=both_bottleneck_aftersplit
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
ll_param_MAP={}
namestr=""
for name in popnamelist:
    namestr+=name
if options.bootstrap!=False:
    def myHandler(signum, frame):
        print("Now, it's the time")
        exit(-1)
    #######
    signal.signal(signal.SIGALRM, myHandler)
    signal.alarm(int(options.bootstrap[1]))
    of=open(namestr+options.tag+options.bootstrap[0]+".parameter","w")
    for name in paramsname:
        ll_param_MAP[name]=[]
    ll_param_MAP["likelihood"]=[]
    ll_param_MAP["theta"]=[]
    func_ex=dadi.Numerics.make_extrap_func(func)
    p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
    popt=dadi.Inference.optimize_log(p0,fsdata,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
    model=func_ex(popt,ns,pts_1)
    theta=dadi.Inference.optimal_sfs_scaling(model,fsdata)
    ll_opt=dadi.Inference.ll_multinom(model,fsdata)
    
    Nref=theta/(4*9.97e-10*float(options.genomelengthwhichsnpfrom))
    for i in range(len(popt)):
        if re.search(r"^T",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[Nref*popt[i]*2,popt[i]]
        elif re.search(r"^m",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[popt[i]/(Nref*2),popt[i]]
        elif re.search(r"^nu",paramsname[i])!=None:
            ll_param_MAP[paramsname[i]]=[popt[i]*Nref,popt[i]]
        else:
            ll_param_MAP[paramsname[i]]=[popt[i],popt[i]]
    ll_param_MAP["likelihood"]=[ll_opt,ll_opt]
    ll_param_MAP["theta"]=[theta,theta]
    for a in sorted(ll_param_MAP.keys()):
        print >>of,a,ll_param_MAP[a][0],ll_param_MAP[a][1]
    else:
        pass
#         print >>of,""
#     for i in range(len(ll_param_MAP["likelihood"])):
#         for a in sorted(ll_param_MAP.keys()):
#             print >>of,ll_param_MAP[a],
#         else:
#             print >>of,""
    if len(popnamelist)==2:
        dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1,residualfilenamepre=options.fsfile+namestr+options.tag+options.model)
    of.close()
#     exit()
else:
    func_ex=dadi.Numerics.make_extrap_func(func)
    p0=dadi.Misc.perturb_params(params,lower_bound=lower_bound,upper_bound=upper_bound)
    popt=dadi.Inference.optimize_log(p0,fsdata,func_ex,pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(params))
    model=func_ex(popt,ns,pts_1)
    theta=dadi.Inference.optimal_sfs_scaling(model,fsdata)
    print theta
    ll_opt=dadi.Inference.ll_multinom(model,fsdata)
    Nref=theta/(4*9.97e-10*float(options.genomelengthwhichsnpfrom))
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
     
    
    
    if len(popnamelist)==2:
        print("print figure")
        pylab.figure()
        dadi.Plotting.plot_single_2d_sfs(fsdata,vmin=1)
        pylab.show()
        pylab.savefig('fsdata_split'+namestr+options.tag+options.model+'.png', dpi=600)
        
        pylab.figure()
        dadi.Plotting.plot_single_2d_sfs(model,vmin=1)
        pylab.show()
        pylab.savefig('model_split'+namestr+options.tag+options.model+'.png', dpi=600)
        
        pylab.figure()
        dadi.Plotting.plot_2d_comp_multinom(model,fsdata,vmin=1,residualfilenamepre=options.fsfile+namestr+options.tag+options.model)
        pylab.show()
        pylab.savefig('compare_split'+namestr+options.tag+options.model+'.png', dpi=600)
        pylab.figure()
        dadi.Plotting.plot_2d_comp_Poisson(model,fsdata,vmin=1)
        pylab.show()
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
    
