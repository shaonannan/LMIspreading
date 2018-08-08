import numpy as np
import itertools
import matplotlib.pyplot as plt
import Inpututilities_beta01 as inp
from internalpopulation import RecurrentPopulation
from externalpopulation import ExternalPopulation
from simulation import Simulation
from connection import Connection as Connection


"""
02/02/2018 version
edited by SYX

initialization  and configuration
vflag for dawson function 
fin for \rho_{Eq}
"""


# intuitive network structure
Net_settings = {'hyp_num': 1,
                'xhyp':1,
                'yhyp':1,
                'xn':4,
                'yn':4, # 25 subpopulation each 20 cells 20 * 25 = 500 cells per hyp!
                'dxx_hyp':500,
                'nmax': 0,
                'Final_time':2000,
                'dt':0.1,
                'dv':1e-3}
# here we use orientation and phase
AmpHyp = Net_settings['xn'] * Net_settings['yn'] /2.0
Fun_map_settings = {'ori_num':4,
                   'phase_num':4}
# cell numbers within identical orientation hypercolumn
#Cell_type_num = {'e':int(512/AmpHyp),
#                'i':int(498/AmpHyp)}
NPHA = Fun_map_settings ['phase_num']
Cell_type_num = {'e':int(1024/AmpHyp/NPHA),
                'i':int(1024/AmpHyp/NPHA)}

print(Net_settings['hyp_num'])
ori_2d = np.zeros((Net_settings['xhyp'] * Net_settings['xn'],Net_settings['yhyp'] * Net_settings['yn']))
ori    = np.zeros(Net_settings['hyp_num'] * Net_settings['xn'] * Net_settings['yn'])
interval_ori = 180.0/Fun_map_settings['ori_num']
"""
this part is for generating orientation map (pin-wheel structure)
if np.mod(Net_settings['nx'],2) == 0,upper side cut off = np.floor(Net_settings['nx']/2), cen = cut-off - 0.5 
[0:cut-off] upper side; while [cut-off:] bottom side
if np.mod(Net_settings['nx'],2) == 1, upper side cut off = np.floor((Net_settings['nx']-1)/2), cen = cut-off
[0:cut-off] upper side; maintain unchanged [cut-off]; while bottom side [cut-off+1:]
"""
nbin = Net_settings['xn']
if np.mod(nbin,2) == 0:
    cut_off = np.floor(nbin/2.0)
    cen     = cut_off - 0.5
    (uc,up,bottom) = (int(cut_off),int(cut_off),int(cut_off))
if np.mod(nbin,2) == 1:
    cut_off = np.floor((nbin-1)/2.0)
    cen     = cut_off
    (uc,up,bottom) = (int(cut_off),int(cut_off),int(cut_off+1))
for i in np.arange(1):
    for p in np.arange(Net_settings['xn']):
        for q in np.arange(Net_settings['yn']):
            ttindex = q + p * Net_settings['yn'] + i * Net_settings['yn'] * Net_settings['xn']
            oriabs  = (np.arctan2(p-cen,q-cen)) * 180.0 / np.pi # ttindex
            if oriabs<0.0:
                oriabs = oriabs + 360.0
            oriabs = oriabs * 0.5
            oriabs = np.floor(oriabs/interval_ori)  # for integer not real DEGREE
            if oriabs == Fun_map_settings['ori_num']:
                oriabs = oriabs - 1
            ori[ttindex] = oriabs
            ori_2d[p,q]  = ori[ttindex] # only half(or even smaller) of the total 
"""
"""
if(Net_settings['xhyp'] > 1):
    # 2nd hypercolum! maybe next line~
    # unchanged
    # because of the next line, so we should use x from 'xn' to 2*'xn' y stay in 0:yn
    ori_2d[Net_settings['xn']+uc,:Net_settings['yn']] = ori_2d[uc,:Net_settings['yn']]
    start = 1*Net_settings['xn'] * Net_settings['yn']
    ori_2d[Net_settings['xn']*1:Net_settings['xn']*1+up,:Net_settings['yn']] = ori_2d[Net_settings['xn']-1:bottom-1:-1,:Net_settings['yn']]
    ori_2d[Net_settings['xn']*1+bottom:Net_settings['xn']*2,:Net_settings['yn']] = ori_2d[up-1::-1,:Net_settings['yn']]
    #ori[start:start+Net_settings['xn']*Net_settings['yn']] = np.squeeze(np.reshape(ori_2d[Net_settings['xn']*1:Net_settings['xn']*2,:],
                                                                    #   (1,Net_settings['xn']**2)))
    
    # then here are two templates for 2*k / 2*k+1
    # for 2*k
    ori_template_2k  = np.squeeze(ori_2d[:Net_settings['xn'],:Net_settings['yn']])
    ori_template_2k1 = np.squeeze(ori_2d[Net_settings['xn']:2*Net_settings['xn'],:Net_settings['yn']])
    
    # only take effects on x(long) axis, means being constrainted in the first row
    
    for ixhyp in range(2,Net_settings['xhyp']):
        if np.mod(ixhyp,2) == 0:
            ori_2d[Net_settings['xn']*ixhyp:Net_settings['xn']*(ixhyp+1),:Net_settings['yn']] = ori_template_2k#np.reshape(ori_template_2k,(3,1))
        if np.mod(ixhyp,2) == 1:
            ori_2d[Net_settings['xn']*ixhyp:Net_settings['xn']*(ixhyp+1),:Net_settings['yn']] = ori_template_2k1#np.reshape(ori_template_2k1,(3,1))

if Net_settings['yhyp']>1:
    # 2nd hypercolum! maybe next column~
    # unchanged
    # because of the next line, so we should use x from 'xn' to 2*'xn' y stay in 0:yn
    ori_2d[:,Net_settings['yn']+uc] = ori_2d[:,uc]
    ori_2d[:,Net_settings['yn']*1:Net_settings['yn']*1+up] = ori_2d[:,Net_settings['yn']-1:bottom-1:-1]
    ori_2d[:,Net_settings['yn']*1+bottom:Net_settings['yn']*2] = ori_2d[:,up-1::-1]
    
    
    # then here are two templates for 2*k / 2*k+1
    # for 2*k
    ori_template_2k  = np.squeeze(ori_2d[:,:Net_settings['yn']])
    ori_template_2k1 = np.squeeze(ori_2d[:,Net_settings['yn']:2*Net_settings['yn']])
    
    # only take effects on y(short) axis
    for iyhyp in range(2,Net_settings['yhyp']):
        if np.mod(iyhyp,2) == 0:
            ori_2d[:,Net_settings['yn']*iyhyp:Net_settings['yn']*(iyhyp+1)] = ori_template_2k
        if np.mod(iyhyp,2) == 1:
            ori_2d[:,Net_settings['yn']*iyhyp:Net_settings['yn']*(iyhyp+1)] = ori_template_2k1
ori = np.squeeze(np.reshape(ori_2d,(1,Net_settings['xn']*Net_settings['yn']*Net_settings['hyp_num'])))

''' simple version orientation map'''
'''
ori_interval = np.pi/Net_settings['xn']/Net_settings['yn']
Net_settings['nmax'] = Net_settings['hyp_num'] * Net_settings['xn'] * Net_settings['yn']
ori = np.zeros(Net_settings['nmax'])
for ihyp in range(Net_settings['hyp_num']):
    for ipatch in range(Net_settings['xn'] * Net_settings['yn']):
        idx_tt = ihyp * Net_settings['xn'] * Net_settings['yn'] + ipatch
        ori[idx_tt] = int(ipatch)
'''

def create_functional_columns(Structure_Net,Functional_Map,ori):
    nmax = Structure_Net['hyp_num'] * Structure_Net['xn'] * Structure_Net['yn']
    Structure_Net['nmax'] = nmax

    hypindx = np.arange(Structure_Net['xhyp'])
    hypindy = np.arange(Structure_Net['yhyp'])
    CGindx  = np.arange(Structure_Net['xn'])
    CGindy  = np.arange(Structure_Net['yn'])
    Orientation_map = {}
    Hypercol_map    = {}
    Hypercol        = np.zeros((nmax,3))
    Phase_map       = {}
    index_2_loc     = np.zeros((nmax,3))
    """
    # this is for single column!!! ~ if there are more than one (physical) column, things will be
    # different
    # orientation map
    for ihyp in hypind:
        for ix in CGindx:
            for iy in CGindy:
                ttindex = iy + ix * Structure_Net['yn'] + ihyp * Structure_Net['xn'] * Structure_Net['yn']
                Hypercol_map[(ttindex,'e')] = ihyp
                Hypercol_map[(ttindex,'i')] = ihyp
                Hypercol[ttindex]           = ihyp
                Phase_map[(ttindex,'e')] = np.random.randint(Fun_map_settings['phase_num'], size=1)
                Phase_map[(ttindex,'i')] = np.random.randint(Fun_map_settings['phase_num'], size=1)
                Orientation_map[(ttindex,'e')] = ori[ttindex]
                Orientation_map[(ttindex,'i')] = ori[ttindex]
                index_2_loc[ttindex,0] = ihyp
                index_2_loc[ttindex,1] = ihyp * Structure_Net['xn'] + ix
                index_2_loc[ttindex,2] = iy
    """
    tt_hypx = Structure_Net['yhyp'] * Structure_Net['yn'] * Structure_Net['xn']
    tt_linx = Structure_Net['yhyp'] * Structure_Net['yn']
    tt_lxhy = Structure_Net['yn']
    for ihypx in hypindx:
        for ix in CGindx:
            for ihypy in hypindy:
                for iy in CGindy:
                    ttindex = iy + ihypy * tt_lxhy + ix * tt_linx + ihypx * tt_hypx
                    ttindex_xax = ix + ihypx * Structure_Net['xn']
                    ttindex_yax = iy + ihypy * Structure_Net['yn']
                    # checking ?
                    ttindex_check = ttindex_xax * Structure_Net['yn'] * Structure_Net['yhyp'] + ttindex_yax
                    if (ttindex_check!=ttindex):
                        print('This chessboard-grid is wrong .')
                    #else:
                        #print('Good chessboard-grid.')
                    ttindex_hyp = ihypx * Structure_Net['yhyp'] + ihypy

                    Hypercol_map[(ttindex,'e')] = ttindex_hyp
                    Hypercol_map[(ttindex,'i')] = ttindex_hyp
                    Hypercol[ttindex,2]         = ttindex_hyp
                    Hypercol[ttindex,0]         = ihypx
                    Hypercol[ttindex,1]         = ihypy

                    Phase_map[(ttindex,'e')] = np.random.randint(Fun_map_settings['phase_num'],size = 1)
                    Phase_map[(ttindex,'i')] = np.random.randint(Fun_map_settings['phase_num'],size = 1)
                    Orientation_map[(ttindex,'e')] = ori[ttindex]
                    Orientation_map[(ttindex,'i')] = ori[ttindex]
                    index_2_loc[ttindex,0] = ttindex_hyp
                    index_2_loc[ttindex,1] = ttindex_xax
                    index_2_loc[ttindex,2] = ttindex_yax
    return (Structure_Net,Hypercol,Phase_map,Orientation_map,index_2_loc)

#******************** period distance **********************************************
def preprocess_for_distance_period(index_2_loc,dxx,dyy,Net_settings,Fun_map_settings):
    nmax = Net_settings['nmax']
    # for X distance
    period_x = Net_settings['xn']
    global_x = np.reshape(np.squeeze(index_2_loc[:,1]),(nmax,1))
    global_xmatrix = np.repeat(global_x,nmax,axis=1)
    global_xdis    = np.abs(global_xmatrix - global_xmatrix.T)
    global_xdis    = np.minimum(global_xdis, period_x - global_xdis)**2*dxx
    # for Y distance
    period_y = Net_settings['yn']
    global_y = np.reshape(np.squeeze(index_2_loc[:,2]),(nmax,1))
    global_ymatrix = np.repeat(global_y,nmax,axis=1)
    global_ydis    = np.abs(global_ymatrix - global_ymatrix.T)
    global_ydis    = np.minimum(global_ydis, period_y - global_ydis)**2*dyy
    global_dist    = np.sqrt(global_xdis + global_ydis)
    
    return (global_xdis,global_ydis,global_dist)

# DISTANCE MATRIX ! BAC
#******************** cut-off distance **********************************************
def preprocess_for_distance(index_2_loc,dxx,dyy,Net_settings,Fun_map_settings):
    nmax = Net_settings['nmax']
    # for X distance
    global_x = np.reshape(np.squeeze(index_2_loc[:,1]),(nmax,1))
    global_xmatrix = np.repeat(global_x,nmax,axis=1)
    global_xdis    = np.abs(global_xmatrix - global_xmatrix.T)**2*dxx
    # for Y distance
    global_y = np.reshape(np.squeeze(index_2_loc[:,2]),(nmax,1))
    global_ymatrix = np.repeat(global_y,nmax,axis=1)
    global_ydis    = np.abs(global_ymatrix - global_ymatrix.T)**2*dyy
    global_dist    = np.sqrt(global_xdis + global_ydis)
    
    return (global_xdis,global_ydis,global_dist)

def normal_function(x, mean=0, sigma=1.0):
    """
    Returns the value of probability density of normal distribution N(mean,sigma) at point `x`.
    """
    _normalization_factor = np.sqrt(2 * np.pi)

#    return np.exp(-np.power((x - mean)/sigma, 2)/2) / (sigma * _normalization_factor)
    return np.exp(-np.power((x - mean)/sigma, 2)/1.0) / (sigma * _normalization_factor) # do not divide 2.0

def circular_dist(a, b, period):
    """
    Returns the distance between a and b (scalars) in a domain with `period` period.
    """
    return np.minimum(np.abs(a - b), period - np.abs(a - b))

# Create visual stimuli
def utilhvsvec(x):
    f = np.maximum(x,0)
    return f

def visual_temporal_kernel(t1,t2,dt,tfinal):
    tt = np.arange(0,tfinal,dt)
    f  = tt/(t1**2)*np.exp(-tt/t1) - tt/(t2**2)*np.exp(-tt/t2)
    return f



Net_settings,Hypm,Pham,Orim,index_2_loc = create_functional_columns(Net_settings,Fun_map_settings,ori)            
Net_settings['nmax'] = Net_settings['xn'] * Net_settings['yn'] * Net_settings['hyp_num']
# MODIFYING 'nmax' in Net_settings, which is used in latter code
#Net_settings,Hypm,Pham,Orim,index_2_loc = create_functional_columns(Net_settings,Fun_map_settings,ori)

"""
using in real python code 2018/04/02 version 0
"""
# Simulation settings:
t0 = 0.0
dt = Net_settings['dt']
tf = Net_settings['Final_time']
dv = Net_settings['dv']
verbose = True

update_method = 'approx'
approx_order = 1
tol  = 1e-14

(dt, tfinal) = (dt, tf)



def input_signal_stream(active_resp,t_properties,Net_settings):
    
    (dt, tfinal) = t_properties

    ntt = int(tfinal/dt) # 20 to escape from error
    # Create base-visual stimuli
    External_stimuli_dict = np.ones(ntt) * 0.07#0.45#

    # Create populations:
    background_population_dict = {}
    internal_population_dict   = {}
    CG_range = np.arange(Net_settings['nmax'])
    NPHA = 4
    # base activities
    for layer, celltype in itertools.product(CG_range, ['e', 'i']):  
        for iphase in range(NPHA):
            tmp_Input = 0.9*np.squeeze(active_resp[layer,iphase,:])
            ExtInp = External_stimuli_dict + tmp_Input[:ntt]
            background_population_dict[layer,iphase,celltype] = ExternalPopulation(ExtInp,dt, record=False)
            # internal_population_dict[layer, celltype]  = RecurrentPopulation(dt = dt,v_min=-1.0, v_max=1.0, dv=dv, update_method=update_method, approx_order=approx_order, tol=tol) old version internal population
            internal_population_dict[layer, iphase,celltype]  = RecurrentPopulation(dt = dt,v_min=-1.0, v_max=1.0, dv=dv, update_method=update_method,
                                    approx_order=approx_order, tol=tol, hyp_idx = 0,ei_pop = celltype,NumCell = Cell_type_num[celltype])

    population_list = list(background_population_dict.values()) + list(internal_population_dict.values())
    return (background_population_dict,internal_population_dict,population_list)

omega        = 1/2.0/10 * np.pi # 4Hz
sigma,gamma  = 8, 0.8
Lambda       = 12
ft           = 4.0/1000.0 * np.pi
GaborBank_yx = inp.MomentGabor(Net_settings['xn'],Net_settings['xn'],sigma,gamma,Lambda)

feedflgn_small,feedflgn_yx,dg  = inp.drifting_spatio_grating(sigma, Lambda, 0, gamma,omega,ft,dt,1,tfinal,GaborBank_yx)
# drifting_grating(sigma,Lambda,0,gamma,omega,dt,1,tfinal,GaborBank_yx)
background_population_dict,internal_population_dict,population_list = input_signal_stream(feedflgn_yx,[dt,tfinal],Net_settings)






# PARAMETERS FOR CONNECTIVITY
(denexc,deninh,axnexc,axninh) = (50.0,50.0,200.0,100.0)
sr_exc2exc = np.sqrt(denexc * denexc + axnexc * axnexc)
sr_inh2exc = np.sqrt(denexc * denexc + axninh * axninh)
sr_exc2inh = np.sqrt(deninh * deninh + axnexc * axnexc)
sr_inh2inh = np.sqrt(deninh * deninh + axninh * axninh)

sr_sigma = {'short_range_exc_exc':sr_exc2exc,
           'short_range_exc_inh':sr_exc2inh,
           'short_range_inh_exc':sr_inh2exc,
           'short_range_inh_inh':sr_inh2inh} # source2target
lr_sigma = {'long_range_exc_exc': 1.00*1000,
            'long_range_exc_inh': 1.00*1000}
NE_source = Cell_type_num['e']
NI_source = Cell_type_num['i']

AmpHyp = Net_settings['xn'] * Net_settings['yn'] /2.0   # 2 is original CG patches in identical orientation hypercolumn
AmpOri = Fun_map_settings['ori_num'] / 2.0
# /Fun_map_settings['phase_num'] 
g_norm = {'short_range_exc_exc':1.25/NE_source/Net_settings['xn']/Net_settings['yn']/Fun_map_settings['phase_num'] ,#0.24/NE_source,
           'short_range_exc_inh':1.30/NE_source/Net_settings['xn']/Net_settings['yn']/Fun_map_settings['phase_num'] ,#0.33/NE_source,
           'short_range_inh_exc':-1.32/NI_source/Net_settings['xn']/Net_settings['yn']/Fun_map_settings['phase_num'] ,#-0.44/NI_source,
           'short_range_inh_inh':-0.65/NI_source/Net_settings['xn']/Net_settings['yn']/Fun_map_settings['phase_num'] ,#-0.09/NI_source,
           
           'long_range_exc_exc_fast':0.84/NE_source/AmpHyp * AmpOri/Fun_map_settings['phase_num'] ,#0.364/NE_source,#0.426/NE_source,#0.24/NE_source,
           'long_range_exc_inh_fast':0.82/NE_source/AmpHyp * AmpOri/Fun_map_settings['phase_num'] ,##0.396/NE_source,#0.33/NE_source,
           
           'long_range_exc_exc':0.664/NE_source/AmpHyp * AmpOri * 6.64/Fun_map_settings['phase_num'] ,#2.46,#1.64,#2.56,#2.46,#1.26,#1.28,
           'long_range_exc_inh':0.664/NE_source/AmpHyp * AmpOri * 6.62/Fun_map_settings['phase_num'] }#2.86}#1.96}#2.86}#2.86}#1.64}#1.96}

# FUNCTIONS FOR CONNECTIVITY
def cortical_to_cortical_connection_normalization(global_dist,sr_sigma,lr_sigma,ori_map,hyp_map_tt,Net_settings,Fun_map_settings):
    # for short-range connections normalization
    norm_cortical_connection = {}
    nmax    = Net_settings['nmax']
    # for long range connextions in different hypercolumns
    nori = Fun_map_settings['ori_num']
    # initialize
    norm_cortical_connection['lr','exc2exc'] = np.zeros((nmax,nmax))
    norm_cortical_connection['lr','exc2inh'] = np.zeros((nmax,nmax))
    """
    """
    # the 1st algorithm, using old code  hyp_map == np.squeeze(hyp_map_tt[:,2])
    # the second algorithm, using identical long-range connections hypx_map = np.squeeze(hyp_map_tt[:,0])
    # hypy_map = np.squeeze(hyp_map_tt[:,1])
    hyp_map = np.squeeze(hyp_map_tt[:,2])
    hypx_map = np.squeeze(hyp_map_tt[:,0])
    hypy_map = np.squeeze(hyp_map_tt[:,1])
    # total normalized value
    lrhyp_xdis,lrhyp_ydis,lrhyp_dist,cen_amp = create_longrange_connection(Net_settings,lr_sigma)
    
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>> the 2nd Algorithm >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for ihyp in range(Net_settings['hyp_num']): # within identical orientation hyper-column
        for iori in range(nori):
            index_hyper = np.where(hyp_map == ihyp)
            index_orien = np.where(ori_map == iori)
            index_target = np.intersect1d(index_hyper,index_orien)
            index_source = np.setdiff1d(index_orien,index_hyper)
            ee_curr = np.zeros((nmax,nmax))
            lentarget = len(index_target)
            lensource = len(index_source)
            print('line 679 >>>>>>>>>>>>>>>>>>>>>>>>>>>')
            print('length target:',lentarget,'; length source:',lensource)
            [x,y] = np.meshgrid(index_target,index_source)
            
            A_temp  = (np.squeeze(global_dist[x,y])).T
            A_temp  = np.reshape(A_temp,(lentarget,lensource))
            print('line 555 >>>>>>>>>>>>>>>>>>>>>>>>>>>')
#            print('558 Atemp:',(A_temp))
#            print('559 target, source: ',lentarget,lensource)
            print('size dist-matrix: [',np.shape(A_temp)[0],',',np.shape(A_temp)[1],']')
            
            A_temp  = normal_function(A_temp,0,lr_sigma['long_range_exc_exc'])
            A_sum   = np.reshape(np.sum(A_temp,axis = 1),(lentarget,1))
            #print('sum',A_sum)
            A_sum   = np.ones_like(A_sum)/A_sum
            A_sum   = np.repeat(A_sum,lensource,axis = 1)
            
            ee_curr[x.T,y.T] = A_temp * A_sum 
            norm_cortical_connection['lr','exc2exc'] = norm_cortical_connection['lr','exc2exc'] + ee_curr
            
            ei_curr = np.zeros((nmax,nmax))
            lentarget = len(index_target)
            lensource = len(index_source)
            [x,y] = np.meshgrid(index_target,index_source)
            
            A_temp  = (np.squeeze(global_dist[x,y])).T
            A_temp  = np.reshape(A_temp,(lentarget,lensource))
            
            A_temp  = normal_function(A_temp,0,lr_sigma['long_range_exc_inh'])
            A_sum   = np.reshape(np.sum(A_temp,axis = 1),(lentarget,1))
            A_sum   = 1.0/np.repeat(A_sum,lensource,axis = 1)
            
            ei_curr[x.T,y.T] = A_temp * A_sum 
            norm_cortical_connection['lr','exc2inh'] = norm_cortical_connection['lr','exc2inh'] + ei_curr
    """        
    """               
            
    # for short range connextions within identical hypercolumn
    nmax    = Net_settings['nmax']
    nmaxhyp = Net_settings['xn'] * Net_settings['yn']
    
    norm_cortical_connection['sr','exc2exc'] = np.zeros((nmax,nmax))
    norm_cortical_connection['sr','exc2inh'] = np.zeros((nmax,nmax))
    norm_cortical_connection['sr','inh2exc'] = np.zeros((nmax,nmax))
    norm_cortical_connection['sr','inh2inh'] = np.zeros((nmax,nmax))
    for ihyp in range(Net_settings['hyp_num']):
        index_start = nmaxhyp * ihyp 
        index_end   = nmaxhyp * (ihyp + 1)
        
        # exc to exc
        ee_curr = np.zeros((nmax,nmax))
        A_temp  = np.squeeze(global_dist[index_start:index_end,:])
        A_temp  = normal_function(A_temp,0,sr_sigma['short_range_exc_exc'])
        A_sum   = np.reshape(np.sum(A_temp,axis = 1),(nmaxhyp,1))
        A_sum   = 1.0/np.repeat(A_sum,nmax,axis = 1)
        A_norm  = (A_temp * A_sum) # x - target cell y source cell
        ee_curr[index_start:index_end,:] = A_norm
        norm_cortical_connection['sr','exc2exc'] = norm_cortical_connection['sr','exc2exc'] + ee_curr
        
        # exc to inh
        ei_curr = np.zeros((nmax,nmax))
        A_temp  = np.squeeze(global_dist[index_start:index_end,:])
        A_temp  = normal_function(A_temp,0,sr_sigma['short_range_exc_inh'])
        A_sum   = np.reshape(np.sum(A_temp,axis = 1),(nmaxhyp,1))
        A_sum   = 1.0/np.repeat(A_sum,nmax,axis = 1)
        A_norm  = (A_temp * A_sum) # x - target cell y source cell
        ei_curr[index_start:index_end,:] = A_norm
        norm_cortical_connection['sr','exc2inh'] = norm_cortical_connection['sr','exc2inh'] + ei_curr
        
        # inh to exc
        ie_curr = np.zeros((nmax,nmax))
        A_temp  = np.squeeze(global_dist[index_start:index_end,:])
        A_temp  = normal_function(A_temp,0,sr_sigma['short_range_inh_exc'])
        A_sum   = np.reshape(np.sum(A_temp,axis = 1),(nmaxhyp,1))
        A_sum   = 1.0/np.repeat(A_sum,nmax,axis = 1)
        A_norm  = (A_temp * A_sum) # x - target cell y source cell
        ie_curr[index_start:index_end,:] = A_norm
        norm_cortical_connection['sr','inh2exc'] = norm_cortical_connection['sr','inh2exc'] + ie_curr
        
        # inh to inh
        ii_curr = np.zeros((nmax,nmax))
        A_temp  = np.squeeze(global_dist[index_start:index_end,:])
        A_temp  = normal_function(A_temp,1,sr_sigma['short_range_inh_inh'])
        A_sum   = np.reshape(np.sum(A_temp,axis = 1),(nmaxhyp,1))
        A_sum   = 1.0/np.repeat(A_sum,nmax,axis = 1)
        A_norm  = (A_temp * A_sum) # x - target cell y source cell
        ii_curr[index_start:index_end,:] = A_norm
        norm_cortical_connection['sr','inh2inh'] = norm_cortical_connection['sr','inh2inh'] + ii_curr
        
    return norm_cortical_connection
    
def cortical_to_cortical_connection(background_population_dict, internal_population_dict, g_norm, delay, orientations,
                                    phases,cortical_norm):
    connection_list = []
    """
    Could be used to create DEE/DIE/DEI/DII/(DEEwithLEEf/DIEwithLIEf)
    """
    nmax = Net_settings['nmax']
    NPHA = 4
    DEE,DEI,DIE,DII = np.zeros((nmax,nmax,NPHA,NPHA)),np.zeros((nmax,nmax,NPHA,NPHA)),np.zeros((nmax,nmax,NPHA,NPHA)),np.zeros((nmax,nmax,NPHA,NPHA))

    # feedforward connections
    for ihyp in range(Net_settings['hyp_num']):
        num_in_hyp = Net_settings['xn'] * Net_settings['yn']
        start_id = ihyp * num_in_hyp
        end_id   = (ihyp + 1) * num_in_hyp
        
        for ixs in range(Net_settings['xn']):
            for iys in range(Net_settings['yn']):
                sub_index_seed = iys + ixs * ( Net_settings['yn'])
                g_index_seed   = sub_index_seed + start_id
                for iphase in range(NPHA):
                    # source_exc, target_exc
                    source_population = background_population_dict[g_index_seed,iphase,'e']
                    target_population = internal_population_dict[g_index_seed,iphase,'e']
                    # curr_connection = Connection(source_population,target_population,1.0, weights = 0.118,probs = 1.0,conn_type = 'ShortRange')
                    curr_connection = Connection(source_population,target_population,nsyn = 1.0,nsyn_post = Cell_type_num['e'], 
                                                 weights = 0.132,probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)#0.030#
                    connection_list.append(curr_connection)
                    
                    # source_inh, target_inh
                    source_population = background_population_dict[g_index_seed,iphase,'i']
                    target_population = internal_population_dict[g_index_seed,iphase,'i']
                    # curr_connection = Connection(source_population,target_population,1.0, weights = 0.117,probs = 1.0,conn_type = 'ShortRange')
                    curr_connection = Connection(source_population,target_population,nsyn = 1.0,nsyn_post = Cell_type_num['i'], 
                                                 weights = 0.1314,probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)#0.0290#
                    connection_list.append(curr_connection)
    """
    """
    # short-range connectivity nomalized 
    ge2e    = g_norm['short_range_exc_exc']
    ge2i    = g_norm['short_range_exc_inh']
    gi2e    = g_norm['short_range_inh_exc']
    gi2i    = g_norm['short_range_inh_inh']  
    
    gle2ef   = g_norm['long_range_exc_exc_fast']
    gle2if   = g_norm['long_range_exc_inh_fast']   
    
    gle2e   = g_norm['long_range_exc_exc']
    gle2i   = g_norm['long_range_exc_inh']              

    
    short_range_exc2exc = cortical_norm['sr','exc2exc']
    short_range_exc2inh = cortical_norm['sr','exc2inh']
    short_range_inh2exc = cortical_norm['sr','inh2exc']
    short_range_inh2inh = cortical_norm['sr','inh2inh']
    
    long_range_exc2exc  = cortical_norm['lr','exc2exc']
    long_range_exc2inh  = cortical_norm['lr','exc2inh']
    """
    """
    # self-connectivity (only when g_index_seed == target_index)
    for ihyp in range(Net_settings['hyp_num']):
        num_in_hyp = Net_settings['xn'] * Net_settings['yn']
        start_id = ihyp * num_in_hyp
        end_id   = (ihyp + 1) * num_in_hyp
        
        for ixs in range(Net_settings['xn']):
            for iys in range(Net_settings['yn']):
                sub_index_seed = iys + ixs * ( Net_settings['yn'])
                g_index_seed   = sub_index_seed + start_id
                i_up_triangle  = g_index_seed

                # source_exc, target_exc
                if(short_range_exc2exc[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                    for  isphase in range(NPHA):
                        for itphase in range(NPHA):
                            source_population = internal_population_dict[g_index_seed,isphase,'e']
                            target_population = internal_population_dict[i_up_triangle,itphase,'e']
                            # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2e * 
                            #                              short_range_exc2exc[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')
                            curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['e'], 
                                                     weights = ge2e * short_range_exc2exc[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                            DEE[i_up_triangle,g_index_seed,isphase,itphase] += ge2e * short_range_exc2exc[i_up_triangle,g_index_seed]

                            connection_list.append(curr_connection)

                # source_exc, target_inh
                if(short_range_exc2inh[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                    for  isphase in range(NPHA):
                        for itphase in range(NPHA):
                            source_population = internal_population_dict[g_index_seed,isphase,'e']
                            target_population = internal_population_dict[i_up_triangle,itphase,'i']
                            # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2i * 
                            #                              short_range_exc2inh[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                            curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['i'], 
                                                     weights = ge2i * short_range_exc2inh[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                            DIE[i_up_triangle,g_index_seed,isphase,itphase] += ge2i  * short_range_exc2inh[i_up_triangle,g_index_seed]

                            connection_list.append(curr_connection)

                # source_inh, target_exc
                if(short_range_inh2exc[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                    for  isphase in range(NPHA):
                        for itphase in range(NPHA):                    
                            source_population = internal_population_dict[g_index_seed,isphase,'i']
                            target_population = internal_population_dict[i_up_triangle,itphase,'e']
                            # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2e *
                            #                              short_range_inh2exc[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                            curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['e'], 
                                                     weights = gi2e * short_range_inh2exc[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                            DEI[i_up_triangle,g_index_seed,isphase,itphase] += gi2e  * short_range_inh2exc[i_up_triangle,g_index_seed]

                            connection_list.append(curr_connection)

                # source_inh, target_inh
                if(short_range_inh2inh[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                    for  isphase in range(NPHA):
                        for itphase in range(NPHA):   
                            source_population = internal_population_dict[g_index_seed,isphase,'i']
                            target_population = internal_population_dict[i_up_triangle,itphase,'i']
                            # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2i * 
                            #                              short_range_inh2inh[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                            curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['i'], 
                                                     weights = gi2i * short_range_inh2inh[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                            DII[i_up_triangle,g_index_seed,isphase,itphase] += gi2i  * short_range_inh2inh[i_up_triangle,g_index_seed]

                            connection_list.append(curr_connection)
  
    """
    # long-range connectivity 
    # choose excitatory as !!! source neuron !!!
    for ihyp in range(Net_settings['hyp_num']):
        num_in_hyp = Net_settings['xn'] * Net_settings['yn']
        start_id = ihyp * num_in_hyp
        end_id   = (ihyp + 1) * num_in_hyp
        
        for ixs in range(Net_settings['xn']):
            for iys in range(Net_settings['yn']):
                sub_index_seed = iys + ixs * ( Net_settings['yn'])
                g_index_seed   = sub_index_seed + start_id
                
                long_range_exc2exc  = cortical_norm['lr','exc2exc']
                long_range_exc2inh  = cortical_norm['lr','exc2inh']
                
                # long-range excitatory 2 excitatory
                # already know seed cell index
                # long-range connectivity [target,source]
                # choose excitatory subpopulation as source neuron
                total_exc2exc = np.squeeze(long_range_exc2exc[:,g_index_seed])
                total_exc2inh = np.squeeze(long_range_exc2inh[:,g_index_seed])
                # effective targetneuron
                effective_target_exc = np.where(total_exc2exc)
                effective_target_inh = np.where(total_exc2inh)
                
                for itarget in range(len(effective_target_exc[0])):
                    icell = effective_target_exc[0][itarget]
                    checkicell = effective_target_inh[0][itarget]
                    if(icell!=checkicell):
                        print('Cell-index for Long-range connectivity is wrong!\n')
                    if(long_range_exc2exc[icell,g_index_seed]>(4.7*1e-6)):
                        source_population = internal_population_dict[g_index_seed,'e']
                        # excitatory target
                        target_population = internal_population_dict[icell,'e']

                        # curr_connection = Connection(source_population,target_population,Cell_type_num['e'],
                        #                              weights = gle2e * long_range_exc2exc[icell,g_index_seed], probs = 1.0, conn_type = 'LongRange' )

                        curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['e'], 
                                             weights = gle2e * long_range_exc2exc[icell,g_index_seed],probs = 1.0,conn_type = 'LongRange',v_min = -1.0,v_max = 1.0,dv = dv)

                        connection_list.append(curr_connection)
#                        print('; No.',itarget,' target excitatory cell: ',icell)
                        
                        # long-range fast synaptic input

                        # curr_connection = Connection(source_population,target_population,Cell_type_num['e'],
                        #                              weights = gle2ef * long_range_exc2exc[icell,g_index_seed], probs = 1.0, conn_type = 'ShortRange' )

                        curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['e'], 
                                             weights = gle2ef * long_range_exc2exc[icell,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                        DEE[icell,g_index_seed] += gle2ef * long_range_exc2exc[icell,g_index_seed]

                        connection_list.append(curr_connection)
                        
                    # inhibitory target
                    if(long_range_exc2inh[checkicell,g_index_seed]>(4.7*1e-6)):
                        source_population = internal_population_dict[g_index_seed,'e']
                        target_population = internal_population_dict[checkicell,'i']
                        # curr_connection = Connection(source_population,target_population,Cell_type_num['e'],
                        #                              weights = gle2i * long_range_exc2inh[checkicell,g_index_seed], probs = 1.0, conn_type = 'LongRange')

                        curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['i'], 
                                             weights = gle2i * long_range_exc2inh[checkicell,g_index_seed],probs = 1.0,conn_type = 'LongRange',v_min = -1.0,v_max = 1.0,dv = dv)

                        connection_list.append(curr_connection)
#                        print(' inhibitory cell: ',checkicell)
                        
                        # long-range fast synaptic input
                        # curr_connection = Connection(source_population,target_population,Cell_type_num['e'],
                        #                              weights = gle2if * long_range_exc2inh[checkicell,g_index_seed], probs = 1.0, conn_type = 'ShortRange')

                        curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['i'], 
                                             weights = gle2if * long_range_exc2inh[checkicell,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                        DIE[checkicell,g_index_seed] += gle2if * long_range_exc2inh[checkicell,g_index_seed]

                        connection_list.append(curr_connection)
#                print(' \n')
    """

    nmax = Net_settings['hyp_num'] * Net_settings['xn'] * Net_settings['yn']
    # short-range connectivity, different patches, not 'delta' self-connectivity
    for ihyp in range(Net_settings['hyp_num']):
        num_in_hyp = Net_settings['xn'] * Net_settings['yn']
        start_id = ihyp * num_in_hyp
        end_id   = (ihyp + 1) * num_in_hyp
        
        for ixs in range(Net_settings['xn']):
            for iys in range(Net_settings['yn']):
                sub_index_seed = iys + ixs * ( Net_settings['yn'])
                g_index_seed   = sub_index_seed + start_id
                for i_up_triangle in range(g_index_seed+1,nmax):#end_id):

                    # source_exc, target_exc
                    if(short_range_exc2exc[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):
                                source_population = internal_population_dict[g_index_seed,isphase,'e']
                                target_population = internal_population_dict[i_up_triangle,itphase,'e']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2e * 
                                #                              short_range_exc2exc[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['e'], 
                                                     weights = ge2e * short_range_exc2exc[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DEE[i_up_triangle,g_index_seed,isphase,itphase] += ge2e * short_range_exc2exc[i_up_triangle,g_index_seed]

                                connection_list.append(curr_connection)
                    
                    # inverse
                    if(short_range_exc2exc[g_index_seed,i_up_triangle]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[i_up_triangle,isphase,'e']
                                target_population = internal_population_dict[g_index_seed,itphase,'e']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2e * 
                                #                              short_range_exc2exc[g_index_seed,i_up_triangle], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['e'], 
                                                     weights = ge2e * short_range_exc2exc[g_index_seed,i_up_triangle],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DEE[g_index_seed,i_up_triangle,isphase,itphase] += ge2e * short_range_exc2exc[g_index_seed,i_up_triangle]

                                connection_list.append(curr_connection)
                    
                    # source_exc, target_inh
                    if(short_range_exc2inh[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[g_index_seed,isphase,'e']
                                target_population = internal_population_dict[i_up_triangle,itphase,'i']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2i * 
                                #                              short_range_exc2inh[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['i'], 
                                                     weights = ge2i * short_range_exc2inh[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)


                                DIE[i_up_triangle,g_index_seed,isphase,itphase] += ge2i * short_range_exc2inh[i_up_triangle,g_index_seed]

                                connection_list.append(curr_connection)
                    
                    if(short_range_exc2inh[g_index_seed,i_up_triangle]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[i_up_triangle,isphase,'e']
                                target_population = internal_population_dict[g_index_seed,itphase,'i']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['e'], weights = ge2i * 
                                #                              short_range_exc2inh[g_index_seed,i_up_triangle], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['e'],nsyn_post = Cell_type_num['i'], 
                                                     weights = ge2i * short_range_exc2inh[g_index_seed,i_up_triangle],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DIE[g_index_seed,i_up_triangle,isphase,itphase] += ge2i * short_range_exc2inh[g_index_seed,i_up_triangle]

                                connection_list.append(curr_connection)
                    
                    # source_inh, target_exc
                    if(short_range_inh2exc[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[g_index_seed,isphase,'i']
                                target_population = internal_population_dict[i_up_triangle,itphase,'e']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2e * 
                                #                              short_range_inh2exc[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['e'], 
                                                     weights = gi2e * short_range_inh2exc[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DEI[i_up_triangle,g_index_seed,isphase,itphase] += gi2e * short_range_inh2exc[i_up_triangle,g_index_seed]

                                connection_list.append(curr_connection)

                    if(short_range_inh2exc[g_index_seed,i_up_triangle]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[i_up_triangle,isphase,'i']
                                target_population = internal_population_dict[g_index_seed,itphase,'e']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2e * 
                                #                              short_range_inh2exc[g_index_seed,i_up_triangle], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['e'], 
                                                     weights = gi2e * short_range_inh2exc[g_index_seed,i_up_triangle],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DEI[g_index_seed,i_up_triangle,isphase,itphase] += gi2e * short_range_inh2exc[g_index_seed,i_up_triangle]

                                connection_list.append(curr_connection)
                    
                    # source_inh, target_inh
                    if(short_range_inh2inh[i_up_triangle,g_index_seed]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):                        
                                source_population = internal_population_dict[g_index_seed,isphase,'i']
                                target_population = internal_population_dict[i_up_triangle,itphase,'i']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2i * 
                                #                              short_range_inh2inh[i_up_triangle,g_index_seed], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['i'], 
                                                     weights = gi2i * short_range_inh2inh[i_up_triangle,g_index_seed],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DII[i_up_triangle,g_index_seed,isphase,itphase] += gi2i * short_range_inh2inh[i_up_triangle,g_index_seed]

                                connection_list.append(curr_connection)

                    if(short_range_inh2inh[g_index_seed,i_up_triangle]>(4.7*1e-6)):
                        for isphase in range(NPHA):
                            for itphase in range(NPHA):
                                source_population = internal_population_dict[i_up_triangle,isphase,'i']
                                target_population = internal_population_dict[g_index_seed,itphase,'i']
                                # curr_connection = Connection(source_population,target_population,Cell_type_num['i'], weights = gi2i * 
                                #                              short_range_inh2inh[g_index_seed,i_up_triangle], probs = 1.0,conn_type = 'ShortRange')

                                curr_connection = Connection(source_population,target_population,nsyn = Cell_type_num['i'],nsyn_post = Cell_type_num['i'], 
                                                     weights = gi2i * short_range_inh2inh[g_index_seed,i_up_triangle],probs = 1.0,conn_type = 'ShortRange',v_min = -1.0,v_max = 1.0,dv = dv)

                                DII[g_index_seed,i_up_triangle,isphase,itphase] += gi2i * short_range_inh2inh[g_index_seed,i_up_triangle]


                                connection_list.append(curr_connection)
    """
    """        
    return connection_list,DEE,DIE,DEI,DII


def create_longrange_connection(Net_settings,sigma_hyplr):
    lrtarget_xaxis, lrtarget_yaxis = 0.0,0.0
    # center 0,0
    xhyp_tt = Net_settings['xhyp'] * 1
    lenx    = xhyp_tt * 2 + 1
    yhyp_tt = Net_settings['yhyp'] * 1 # 0#
    leny    = yhyp_tt * 2 + 1 #1#
    # here we only care about x-axis, not y-axis perpendicular to exact x-axis 

    dxx_hyp2 = Net_settings['dxx_hyp']**2

    xhyp_range = np.reshape(np.arange(-xhyp_tt,xhyp_tt+1),(lenx,1))
    xhyp_range = np.repeat(xhyp_range,leny,axis = 1)

    yhyp_range = np.reshape(np.arange(-yhyp_tt,yhyp_tt+1),(1,leny))
    yhyp_range = np.repeat(yhyp_range,lenx,axis = 0)

    lrhyp_xdis = np.abs(xhyp_range - lrtarget_xaxis * np.ones_like(xhyp_range))
    lrhyp_xdis = lrhyp_xdis**2 * dxx_hyp2

    lrhyp_ydis = np.abs(yhyp_range - lrtarget_yaxis * np.ones_like(yhyp_range))
    lrhyp_ydis = lrhyp_ydis**2 * dxx_hyp2

    lrhyp_dist = np.sqrt(lrhyp_xdis + lrhyp_ydis) # squart

    # normalization long-range connections for ei and ee
    _normalization_factor = np.sqrt(2*np.pi)

    cen_amp = {}
    cen_amp['lr','exc2exc'] = np.zeros((xhyp_tt+1,yhyp_tt+1))
#   print('xhyp:',xhyp_tt,'yhyp:',yhyp_tt)
    cen_amp['lr','exc2inh'] = np.zeros((xhyp_tt+1,yhyp_tt+1))
    # long-range excitatory to excitatory
    mean,sigma = 0.0,sigma_hyplr['long_range_exc_exc']
    exp_amp    = np.exp(-np.power((lrhyp_dist - mean)/sigma, 2)/2.0) / (sigma * _normalization_factor)
    exp_amp[xhyp_tt,yhyp_tt] = 0.0
    sum_amp = np.squeeze(np.reshape(exp_amp,(lenx*leny,1)))
    sum_amp = np.sum(sum_amp)
    exp_amp /= sum_amp
    cen_amp['lr','exc2exc'] = exp_amp[xhyp_tt:,yhyp_tt:]

    # long-range excitatory to inhibitory
    mean,sigma = 0.0,sigma_hyplr['long_range_exc_inh']
    exp_amp    = np.exp(-np.power((lrhyp_dist - mean)/sigma, 2)/2.0) / (sigma * _normalization_factor)
    exp_amp[xhyp_tt,yhyp_tt] = 0.0
    sum_amp = np.squeeze(np.reshape(exp_amp,(lenx*leny,1)))
    sum_amp = np.sum(sum_amp)
    exp_amp /= sum_amp
    cen_amp['lr','exc2inh'] = exp_amp[xhyp_tt:,yhyp_tt:]

    return lrhyp_xdis,lrhyp_ydis,lrhyp_dist,cen_amp



Net_settings,Hypm,Pham,Orim,index_2_loc = create_functional_columns(Net_settings,Fun_map_settings,ori)
(dxx,dyy) = (500.0/Net_settings['xn'],500.0/Net_settings['yn'])
dxx = dxx**2
dyy = dyy**2
global_x,global_y,global_dist = preprocess_for_distance(index_2_loc,dxx,dyy,Net_settings,Fun_map_settings)
cortical_norm = cortical_to_cortical_connection_normalization(global_dist,sr_sigma,lr_sigma,ori,Hypm,Net_settings,Fun_map_settings)

connection_list,DEE,DIE,DEI,DII = cortical_to_cortical_connection(background_population_dict,internal_population_dict,g_norm,0.0,Orim,Pham,
                                             cortical_norm)


"""
"""
simulation = Simulation(population_list, connection_list,Net_settings,Cell_type_num,DEE,DIE,DEI,DII,verbose=True)
(mEbin_ra,mIbin_ra,xEbin_ra,xIbin_ra,VEavg_ra,VIavg_ra,VEstd_ra,VIstd_ra,VEavgbin_ra,VIavgbin_ra,VEstdbin_ra,VIstdbin_ra,nEbin_ra,nIbin_ra) = simulation.update(t0=t0, dt=dt, tf=tf)

import time
import scipy.io as scio 
ISOTIMEFORMAT='%Y%m%d%H%M%S'
filename=str(time.strftime(ISOTIMEFORMAT)) + '_large.mat'
scio.savemat(filename,{'mEbin_ra':mEbin_ra,'mIbin_ra':mIbin_ra,'xEbin_ra':xEbin_ra,'xIbin_ra':xIbin_ra,'VEavg_ra':VEavg_ra,'VIavg_ra':VIavg_ra,'VEstd_ra':VEstd_ra,'VIstd_ra':VIstd_ra,'nEbin_ra':nEbin_ra,'nIbin_ra':nIbin_ra}) 


#tbin = 150
#for i in range(Net_settings['nmax']):
#    ttt = np.arange(150) * 1.0
#    plt.figure(12)
#    if (np.ceil(i/3)<3):
#        plt.plot(ttt,mEbin_ra[:tbin,i],'r')
#    if (np.ceil(i/3)>2)&(np.ceil(i/3)<5):
#        plt.plot(ttt,mEbin_ra[:tbin,i],'m')
#    if np.ceil(i/3)>4:
#        plt.plot(ttt,mEbin_ra[:tbin,i],'y')
#    plt.xlim([0,tbin])
#    plt.ylim([0,20.0])
#    plt.pause(0.1)
#    plt.figure(23)
#    if (np.ceil(i/3)<3):
#        plt.plot(ttt,mIbin_ra[:tbin,i],'r')
#    if (np.ceil(i/3)>2)&(np.ceil(i/3)<5):
#        plt.plot(ttt,mIbin_ra[:tbin,i],'m')
#    if (np.ceil(i/3)>4):
#        plt.plot(ttt,mIbin_ra[:tbin,i],'y')
#    plt.xlim([0,tbin])
#    plt.ylim([0,20.0])
#    plt.pause(0.1)

