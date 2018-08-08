from connectiondistributioncollection import ConnectionDistributionCollection
import time
import numpy as np
import utilities as util
import matplotlib.pyplot as plt



class Simulation(object):
    """
    Parameters:
    list :
        All sub-population (cluster)
        All connection (cluster)
        [type of both is 'List', which is changable variable, and could be changed]
        
    generate after initiate(by hand)
        connection_distribution
        connection_distribution_list
        [the differences between connection, connection_distribution and connection_distribution_list are
        connection ---> the component of 'connection_list', record all information and related information and object,like source and others
        connection_distribution --> this variable is a preparation variable for further processing, each 'connection' could generate a 
        class 'connecton_distribution' and then, using weight,syn,prob, could calculate flux_matrix and threshold
        each 'connection_distribution' item is defined by 'weight''syn ''prob', items with identical symbol will be classified to the same
        distribution
        connection_distribution_list --> this is a 'basket', store all unique connections(definition of unique: unique symbol
        'weight','syn','prob' no matter the target/source population)
    """
    def __init__(self,population_list,connection_list,Net_settings,Cell_type_num,DEE,DIE,DEI,DII,verbose=True):
        
        self.verbose = verbose
        self.population_list = population_list
        self.connection_list = [c for c in connection_list if c.nsyn!=0.0]
        self.Net_settings    = Net_settings
        tfinal = Net_settings['Final_time']
        dt     = Net_settings['dt']
        self.ntt = int(tfinal/dt)
        self.m_record = None
        ''' all for MFE '''
        self.VE,self.VI = None,None
        self.Vedges,self.Vbins = None,None
        self.NE,self.NI = Cell_type_num['e'],Cell_type_num['i']
        self.MFE_num  = 0
        self.MFE_flag = 0
        
        self.DEE = DEE
        self.DEI = DEI
        self.DIE = DIE
        self.DII = DII

    
    def initialize(self,t0=0.0):
        """
        initialize by hand, first put all sub-population and connection-pair
        !!! put them on the same platform!!! simulationBridge
        """
        ''' initialize P_MFE '''
        ''' at first, we only use NHYP as NPATCH '''
        self.iteration_max = self.ntt+100
        iteration_max = self.iteration_max
        self.tbin_tmp = 0 # initial
        self.tbinsize = 1.0
        dtperbin = int(self.tbinsize/self.dt)
        self.dtperbin = dtperbin
        iteration_bin = int(iteration_max/dtperbin)
        NPATCH = self.Net_settings['nmax']
        NE,NI  = self.NE,self.NI
        NPHA   = 4
        self.VE,self.VI = np.zeros((NE,NPHA*NPATCH)),np.zeros((NI,NPHA*NPATCH))
        # DTBIN_RECORD_FLAG
        self.tbin_ra = np.zeros((iteration_max,1))
        self.mE_ra   = np.zeros((iteration_max,NPHA*NPATCH))
        self.mI_ra   = np.zeros((iteration_max,NPHA*NPATCH))
        self.mEbin_ra = np.zeros((iteration_bin,NPHA*NPATCH))
        self.mIbin_ra = np.zeros_like(self.mEbin_ra)
        self.xEbin_ra = np.zeros_like(self.mEbin_ra)
        self.xIbin_ra = np.zeros_like(self.xEbin_ra)
        
        self.nEbin_ra = np.zeros_like(self.xEbin_ra)
        self.nIbin_ra = np.zeros_like(self.xEbin_ra)


        self.VEavgbin_ra = np.zeros((iteration_bin,NPHA*NPATCH))
        self.VIavgbin_ra = np.zeros_like(self.VEavgbin_ra)
        self.VEstdbin_ra = np.zeros_like(self.VIavgbin_ra)
        self.VIstdbin_ra = np.zeros_like(self.VEstdbin_ra)
        
        self.VEavg_ra = np.zeros((iteration_max,NPHA*NPATCH))
        self.VIavg_ra = np.zeros_like(self.VEavg_ra)
        self.VEstd_ra = np.zeros_like(self.VIavg_ra)
        self.VIstd_ra = np.zeros_like(self.VEstd_ra)
        self.rE,self.rI  = None,None
        self.NPATCH = NPATCH
        self.NPHA   = NPHA
    



        dv = self.Net_settings['dv']
        self.Vedges = util.get_v_edges(-1.0,1.0,dv)
        ''' bins = edges - 1'''
        # in internal , rhov length len(self.Vbins), len(Vedges)-1
        self.Vbins = 0.5*(self.Vedges[0:-1] + self.Vedges[1:]) 
        Vedges = self.Vedges 
        Vbins  = self.Vbins

        self.rE = np.zeros((len(self.Vbins),self.NPATCH*self.NPHA))
        self.rI = np.zeros_like(self.rE)

        
        # An connection_distribution_list (store unique connection(defined by weight,syn,prob))
        self.connection_distribution_collection = ConnectionDistributionCollection() # this is 
        self.t = t0

        # Matrix to record 
        numCGPatch = self.Net_settings['nmax'] * 2 # excitatory and inhibitory
        # 2 * numCGPatch = External Population and Recurrent Population
        # set Matrix to record only Internal Population
        self.m_record = np.zeros((NPHA*(numCGPatch+1), self.ntt + 10))
        self.v_record = np.zeros_like(self.m_record)
        
        # put all subpopulation and all connections into the same platform
        for subpop in self.population_list:
            subpop.simulation = self    # .simulation = self(self is what we called 'simulation')
        for connpair in self.connection_list:
            connpair.simulation = self
            
        # initialize population_list, calculate         
        for p in self.population_list:
            p.initialize()      # 2   
        
        for c in self.connection_list:
            #print 'initialize population'
            c.initialize()      # 1
  
    def update(self,t0,dt,tf):
        self.dt = dt
        self.tf = tf   
        # initialize:
        start_time = time.time()
        self.initialize(t0)
        self.initialize_time_period = time.time()-start_time
        
        # start_running
        start_time = time.time()
        counter = 0
        numCGPatch = self.Net_settings['nmax']*2
        print('nET:',self.Net_settings['nmax']*2)
        '''
        at first 
        '''
        Vbins,Vedges,NPATCH,NPHA = self.Vbins,self.Vedges,self.NPATCH,self.NPHA
        while self.t < self.tf:
            # refresh current time as well as current time-step
            self.t+=self.dt
            self.tbin_tmp = int(np.floor(self.t/self.tbinsize))
#            if (self.t > 100.0) &(flag_change==1):
#                ind_rec = 0
#                for p in self.population_list:
#                    ind_rec += 1
#                    if ind_rec <= numCGPatch:
#                        p.cu
            #if self.verbose: print ('time: %s' % self.t)
            ind_rec,idxE,idxI = 0,0,0   # start to accumulate index of hypercolumn

            for p in self.population_list:
                ''' updating OP 2 modes: updating under Moment/updating under MFE '''
                # updating under Moment- full
                p.USUALorMFE = 1
                ind_rec += 1
                '''
                Recording at first, before p.update(),
                rE and rI purely from(after) MFE should be recorded in rE/I(bin)_ra, rather
                than RvE from Moment
                '''
                # before Moment iteration
                if(ind_rec>numCGPatch*self.NPHA): # means internal-population, not external-population
                    if p.ei_pop == 'e':
                        ''' Voltage distribution should be recorded each dt as well as each dtbin'''
                        # dt-recording
                        self.VEavg_ra[counter,idxE] = p.v1
                        self.VEstd_ra[counter,idxE] = np.sqrt(p.v2-p.v1**2)
                        # dtbin-recording
                        self.VEavgbin_ra[self.tbin_tmp,idxE] += p.v1*dt
                        self.VEstdbin_ra[self.tbin_tmp,idxE] += np.sqrt(p.v2-p.v1**2)*dt

                        idxE += 1

                            

                    else:
                        ''' Voltage distribution should be recorded each dt as well as each dtbin'''
                        # dt-recording VE/Iavg
                        self.VIavg_ra[counter,idxI] = p.v1
                        self.VIstd_ra[counter,idxI] = np.sqrt(p.v2-p.v1**2)
                        # dtbin-recording
                        self.VIavgbin_ra[self.tbin_tmp,idxI] += p.v1*dt
                        self.VIstdbin_ra[self.tbin_tmp,idxI] += np.sqrt(p.v2-p.v1**2)*dt
                        
                        
                        idxI += 1

                        
                p.update()
                '''
                when using USUALorMFE==1
                updating rhov as well as firing rate
                
                next, should record firing rate mE/I in mE/I(bin)_ra
                [but not rE/I(bin)_ra]
                
                and also, RvE/I were extracted out from p-list, which were used
                to calculate MFE probability                   
                
                '''
                if(ind_rec>numCGPatch*self.NPHA):
                    # print('num: ',numCGPatch,np.shape(self.v_record))
                    #print('ind_rec:',ind_rec,'firing rate:',p.curr_firing_rate)
                    self.v_record[ind_rec-numCGPatch*self.NPHA,counter] = p.v1
                    self.m_record[ind_rec-numCGPatch*self.NPHA,counter] = p.curr_firing_rate#
                if(counter>0):
                    if(ind_rec>numCGPatch*self.NPHA):
                        if p.ei_pop == 'e': 
                            continue
                            #print('excite : %.5f'%p.local_pevent)
            ind_rec,idxE,idxI = 0,0,0                
            for p in self.population_list:
                ind_rec += 1
                if(counter>0):
                    if(ind_rec>numCGPatch*self.NPHA):
                        if p.ei_pop == 'e': 
                            '''
                            and also extract curr_rhov to calculate PMFE
                            '''
                            self.rE[:,idxE] = p.curr_rhov
#                            print('pre: ',p.firing_rate)
#                            p.firing_rate = 0.0
#                            print('pos: ',p.curr_firing_rate)
                            '''
                            here, recording new firing rate 
                            mE/I_ra
                            and also extract curr_rhov to calculate PMFE
                            '''
                            self.mE_ra[counter,idxE] = p.curr_firing_rate
                           
                            idxE += 1
                        else:
                            self.mI_ra[counter,idxI] = p.curr_firing_rate
                            self.rI[:,idxI] = p.curr_rhov                           
                            
                            idxI += 1
                
            

            for c in self.connection_list:
                c.update()
            counter +=1
            
            ''' recording ! '''
                    # DTBIN_RECORD_FLAG
            self.tbin_ra[counter] = np.floor(self.t/self.tbinsize)
            tbin = int(np.floor(self.t/self.tbinsize))
            ind_rec,idxE,idxI   = 0,0,0
            NE = self.NE
            NI = self.NI
            for p in self.population_list:
                ind_rec +=1
                if(counter>0):
                    if(ind_rec>numCGPatch*self.NPHA):
                        if p.ei_pop == 'e':
                            self.mEbin_ra[tbin,idxE] +=  p.curr_firing_rate * NE * dt 
                            self.xEbin_ra[tbin,idxE] += util.psample(p.curr_firing_rate * NE * dt) 
                            
                            self.nEbin_ra[tbin,idxE] += p.total_Inmda * NE * dt
                            idxE += 1
                        else:
                            self.mIbin_ra[tbin,idxI] += p.curr_firing_rate * NE * dt 
                            self.xIbin_ra[tbin,idxI] += util.psample(p.curr_firing_rate * NE * dt) 
                            
                            self.nIbin_ra[tbin,idxI] += p.total_Inmda * NI * dt
                            idxI += 1
                            
            ''' visualizing '''
            if np.mod(counter,50) < 1:
                if np.mod(counter,50) == 0:
                    print("t_sum: ",counter * self.dt)
                for i in range(0,NPATCH*self.NPHA,self.NPHA):
                    print('Excitatory pop %d :%.4f'%(i,self.VEavgbin_ra[tbin,i]))
                    print('Inhibitory pop %d :%.4f'%(i,self.VIavgbin_ra[tbin,i]))
                    ttt = np.arange(tbin) * 1.0
                    plt.figure(220)
                    if (i<4*self.NPHA):
                        plt.plot(ttt,self.VEavgbin_ra[:tbin,i],'r')
                    if (i>=4*self.NPHA)&(i<4*self.NPHA*2):
                        plt.plot(ttt,self.VEavgbin_ra[:tbin,i],'m')
                    if (i>=4*self.NPHA*2)&(i<4*self.NPHA*3):
                        plt.plot(ttt,self.VEavgbin_ra[:tbin,i],'g')
                    if (i>=4*self.NPHA*3):
                        plt.plot(ttt,self.VEavgbin_ra[:tbin,i],'b')
                    plt.xlim([0,int(self.tf)])
                    plt.ylim([0,1.2])
                    plt.pause(0.1)
                    plt.figure(222)
                    if (i<4*self.NPHA):
                        plt.plot(ttt,self.VIavgbin_ra[:tbin,i],'r')
                    if (i>=4*self.NPHA)&(i<4*self.NPHA*2):
                        plt.plot(ttt,self.VIavgbin_ra[:tbin,i],'m')
                    if (i>=4*self.NPHA*2)&(i<4*self.NPHA*3):
                        plt.plot(ttt,self.VIavgbin_ra[:tbin,i],'g')
                    if (i>=4*self.NPHA*3):
                        plt.plot(ttt,self.VIavgbin_ra[:tbin,i],'b')
                    plt.xlim([0,int(self.tf)])
                    plt.ylim([0,1.2])
                    plt.pause(0.1)
                    
                           
        return self.mEbin_ra,self.mIbin_ra,self.xEbin_ra,self.xIbin_ra,self.VEavg_ra,self.VIavg_ra,self.VEstd_ra,self.VIstd_ra,self.VEavgbin_ra,self.VIavgbin_ra,self.VEstdbin_ra,self.VIstdbin_ra,self.nEbin_ra,self.nIbin_ra