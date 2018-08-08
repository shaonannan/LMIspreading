import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


def MomentGabor(nxhyp,nyhyp,sigma,gamma,Lambda):
    ''' should be Even '''
    NPHA = 4
    if np.mod(nxhyp,2)==1:
        print('re-choose,should be Even!!!')
    else:
        GaborBank = {}
        center_x,center_y = nxhyp/2-0.5, nyhyp/2-0.5
        for iraw in range(nxhyp):
            for icol in range(nyhyp):
                idxtt = icol + iraw * nyhyp
                tgv = iraw - center_x
                tgh = icol - center_y
                tgtheta = 0.5*np.arctan2(tgh,tgv)
                
                for iphase in range(NPHA):
                    randphase = 2*np.pi*iphase/NPHA/1.0
                    GaborBank[idxtt,iphase] = gabor_fn(sigma, tgtheta,Lambda,randphase,gamma)
                    '''
                    plt.imshow(GaborBank[idxtt,iphase])
                    plt.pause(0.5)
                    plt.show()
                    '''
    return GaborBank


def gabor_fn(sigma, theta, Lambda, psi, gamma):
    sigma_x = sigma
    sigma_y = float(sigma) / gamma

    # ------这部分内容是为了确定卷积核的大小------
    # Bounding box
    nstds = 1.2 # Number of standard deviation sigma
    thetao = np.pi/4
    xmax = max(abs(nstds * sigma_x * np.cos(thetao)), abs(nstds * sigma_y * np.sin(thetao)))
    xmax = np.ceil(max(1, xmax))
    ymax = max(abs(nstds * sigma_x * np.sin(thetao)), abs(nstds * sigma_y * np.cos(thetao)))
    ymax = np.ceil(max(1, ymax))
    xmin = -xmax
    ymin = -ymax
    (y, x) = np.meshgrid(np.arange(ymin, ymax + 1), np.arange(xmin, xmax + 1))
    # ------这部分内容是为了确定卷积核的大小------

    # Rotation 
    x_theta = x * np.cos(theta) + y * np.sin(theta)
    y_theta = -x * np.sin(theta) + y * np.cos(theta)

    # ------这部分正是上面的公式（1）------
    gb = np.exp(-.5 * (x_theta ** 2 / sigma_x ** 2 + y_theta ** 2 / sigma_y ** 2)) * np.cos(2 * np.pi / Lambda * x_theta + psi)
    return gb

# '''
# First Step!
# Generating Gabor-filter Bank
# '''
# sigma,gamma = 8, 0.8
# Lambda      = 12
# GaborBank_yx = MomentGabor(4,4,sigma,gamma,Lambda)
def drifting_spatio_grating(sigma, Lambda, psi, gamma,omega,ft,dt,dtl,tfinal,GaborBank):
    NPHA = 4
    '''
    ft: time frequency
    lambda: space frequency
    sigma: Bounding box
    '''
    ntt      = int(tfinal/dtl)
    ntt_real = int(tfinal/dt)
    bin_nt   = int(dtl/dt)

    sigma_x = sigma
    sigma_y = float(sigma) / gamma

    nstds = 1.2 # Number of standard deviation sigma
    thetao = np.pi/4.0
    xmax = max(abs(nstds * sigma_x * np.cos(thetao)), abs(nstds * sigma_y * np.sin(thetao)))
    xmax = np.ceil(max(1, xmax))
    ymax = max(abs(nstds * sigma_x * np.sin(thetao)), abs(nstds * sigma_y * np.cos(thetao)))
    ymax = np.ceil(max(1, ymax))
    xmin = -xmax
    ymin = -ymax
    (y, x) = np.meshgrid(np.arange(ymin, ymax + 1), np.arange(xmin, xmax + 1))

    nmax = 16
    

    feedflgn      = np.zeros((nmax,NPHA,ntt))
    feedflgn_real = np.zeros((nmax,NPHA,ntt_real))
    clock_dg      = {}

    omega  = 20/180.0*np.pi/1000

    for itorien in range(ntt):
        '''
        omega-theta : orientation
        ft - phase : phase
        '''
        t_curr = dtl * itorien
        theta  = t_curr * omega
        phase  = t_curr * ft
        psi    = phase
        
        start_nt = max(0,itorien*bin_nt)
        end_nt   = min((itorien+1)*bin_nt,ntt_real)

        x_theta = x * np.cos(theta) + y * np.sin(theta)
        y_theta = -x * np.sin(theta) + y * np.cos(theta)

        gb_dg   = np.exp(-.5 * (x_theta ** 2 / sigma_x ** 2 + y_theta ** 2 / sigma_y ** 2)) * np.cos(2 * np.pi / Lambda * x_theta + psi)
        clock_dg[itorien] = gb_dg

        for imoment in range(nmax):
            Resp_const = 0
            for iphase in range(NPHA):
                # print(np.shape(GaborBank[imoment,iphase]))
                GaborIndividual = (GaborBank[imoment,iphase])    
    
                # ''' fast fake turning dg '''        
                # RespConv        = max(np.reshape(signal.convolve(GaborIndividual,gb_dg,"full"),[-1]))
                # feedflgn[imoment,itorien]              = RespConv
                # feedflgn_real[imoment,start_nt:end_nt] = RespConv
                # # print('size:',nmax,np.shape(GaborIndividual),np.shape(gb_dg),np.shape(RespConv))
    
    
                ''' normal turning dg '''
                RespConv        = (GaborIndividual*gb_dg)
                Resp_const      = max(0,np.sum(np.reshape(RespConv,[-1])))
                feedflgn[imoment,iphase,itorien]              = Resp_const
                feedflgn_real[imoment,iphase,start_nt:end_nt] = Resp_const
    ttfeedf = np.reshape(feedflgn,[-1])
    maxttfeedf = max(ttfeedf)

    for imoment in range(nmax):
        for iphase in range(NPHA):
            feedflgn[imoment,iphase,:] = feedflgn[imoment,iphase,:]/maxttfeedf#/np.max(np.squeeze(feedflgn[imoment,iphase,:]))
            feedflgn_real[imoment,iphase,:] = feedflgn_real[imoment,iphase,:]/maxttfeedf#/np.max(np.squeeze(feedflgn_real[imoment,iphase,:]))

    return feedflgn,feedflgn_real,clock_dg




def drifting_grating(sigma, Lambda, psi, gamma,omega,dt,dtl,tfinal,GaborBank):
    '''
    ft: time frequency
    lambda: space frequency
    sigma: Bounding box
    '''
    ntt      = int(tfinal/dtl)
    ntt_real = int(tfinal/dt)
    bin_nt   = int(dtl/dt)

    sigma_x = sigma
    sigma_y = float(sigma) / gamma

    nstds = 1.2 # Number of standard deviation sigma
    thetao = np.pi/4.0
    xmax = max(abs(nstds * sigma_x * np.cos(thetao)), abs(nstds * sigma_y * np.sin(thetao)))
    xmax = np.ceil(max(1, xmax))
    ymax = max(abs(nstds * sigma_x * np.sin(thetao)), abs(nstds * sigma_y * np.cos(thetao)))
    ymax = np.ceil(max(1, ymax))
    xmin = -xmax
    ymin = -ymax
    (y, x) = np.meshgrid(np.arange(ymin, ymax + 1), np.arange(xmin, xmax + 1))

    nmax = len(GaborBank)
    

    feedflgn      = np.zeros((nmax,ntt))
    feedflgn_real = np.zeros((nmax,ntt_real))

    for itorien in range(ntt):
        t_curr = dtl * itorien
        theta = t_curr * omega
        
        start_nt = max(0,itorien*bin_nt)
        end_nt   = min((itorien+1)*bin_nt,ntt_real)

        x_theta = x * np.cos(theta) + y * np.sin(theta)
        y_theta = -x * np.sin(theta) + y * np.cos(theta)

        gb_dg   = np.exp(-.5 * (x_theta ** 2 / sigma_x ** 2 + y_theta ** 2 / sigma_y ** 2)) * np.cos(2 * np.pi / Lambda * x_theta + psi)

        for imoment in range(nmax):
            GaborIndividual = (GaborBank[imoment])
            
            RespConv        = max(np.reshape(signal.convolve(GaborIndividual,gb_dg,"full"),[-1]))
            # print('size:',nmax,np.shape(GaborIndividual),np.shape(gb_dg),np.shape(RespConv))
            feedflgn[imoment,itorien]              = RespConv
            
            feedflgn_real[imoment,start_nt:end_nt] = RespConv

    for imoment in range(nmax):
        feedflgn[imoment,:] = feedflgn[imoment,:]/np.max(np.squeeze(feedflgn[imoment,:]))
        feedflgn_real[imoment,:] = feedflgn_real[imoment,:]/np.max(np.squeeze(feedflgn_real[imoment,:]))

    return feedflgn,feedflgn_real

# '''
# Second Step!
# Generating Temporal Response
# '''
# omega       = 1/4.0/1000.0 * np.pi # 4Hz
# ft = 1/4.0/1000.0 * 2 * np.pi
# feedflgn,feedflgn_yx,dg_stm = drifting_spatio_grating(sigma, Lambda, 0, gamma,omega,ft,0.1,1,6000,GaborBank_yx) # drifting_grating(sigma,Lambda,0,gamma,omega,0.1,1,2000,GaborBank_yx)


