
# Written by Chris Calderon 2015 (Chris.Calderon@UrsaAnalytics.com)
#
#
# Copyright 2015 Ursa Analytics, Inc.

   # Licensed under the Apache License, Version 2.0 (the "License");
   # you may not use this file except in compliance with the License.
   # You may obtain a copy of the License at

   #     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__license__ = "Apache License, Version 2.0"
__author__  = "Chris Calderon, Ursa Analytics, Inc. [www.UrsaAnalytics.com]"
__status__  = "Development"

import numpy as np
import scipy.special as spspecial


def langevinScaling(pars,betakT310=139.28553022983505):
    """
    assumes input 1D (k,sig,epslocStd,v) from dXt= (v + kX_t)dt + sigdBt
    w/ betaKT310 = 1/(69.479/10./1000.*310/300)

    """
    parsNew=pars[:] #make implicit copy of vector
#    scale = pars[1]**2/2.*betakT310
#    parsNew[0]=parsNew[0]*scale
#    parsNew[-1]=parsNew[-1]*scale

    return parsNew


class ModifiedKalmanFilter1DwithCrossCorr(object):
    def __init__(self,tsData,tData,StaticErrorEstSeq=[],dt=[]): 
        """
        setup class for "Motion Blur Filter" [MotBnB version] for 1D  case where state evolves in continuous time and "frame gaps" may exist.
        SDE driven by standard Brownian motion and discrete measurements are "blurred"  
        (blurring occurs both from "dynamic" motion blur and "static" point spread function fitting errors).   
        filtering code ignores computational issues inherent to multivariate models.  
        also assume stationary distribution for initial cov vs. information filter (latter is problematic in MLE anyway)

        Note: although code can handle estimation cases with \kappa is near zero, recommended to use MA1 code provided if \kappa is within
        parameter uncertainty of zero since likelihood expansion used for limiting case (analytic form exists, this numerical implementation 
        just switches to taylor likelihood proxy to avoid potential numerical overflow (this introduces some slight bias in estimates obtained for very small kappa).    
        if true  \kappa < 0 in DGP, model isn't "confined" and current code should not be used for this "unstable" case).   improved algorithm could be made to "take limits" vs. taylor proxy

        class setup to allow "traditional" kalman filter by redefining some methods and using some aux variables (modifications for 1D case illustrated below) 

        input:
        tsData:  list [each element is an np.array of time series data]
        tData:   list [each element is list of time stamps (assumed in seconds) for each of the observations].  
        StaticErrorEstSeq [optional]:  list or np.array of same length as tsData with STD DEV est of static errors.  allows fusing in estimations (default zero)
                           if this parameter is input, the addition of the estimated \hat{sigma}^{loc} to the input stream gives a refined estimate of the empirical "static noise" 
        dt    [optional]:   exposure time of camera.  if not provided, smallest time spacing between observations assumed to be exposure time.  can be scalar float.
        output: 
        likelihood: log likelihood of innovations
        xfilt:      filtered state for given parameters  
        pit:        PIT associated with likelihood function evaluated at given parameters
        Shist:      history of innovation covariance evaluated at given parameters.  
        """
        
        #error check input for more stringent list of numpy array type in "batch" version 
        if all(isinstance(x,np.ndarray) for x in tsData) and all(isinstance(x,np.ndarray) for x in tData):
            # self._T=max(tsData.shape) #compute length of time series
            self._T=[max(i.shape) for i in tsData] #compute length of each time series in batch and store as list
            # self._y=np.reshape(tsData,(self._T,1),order='F')#reshape data to column vector (scalar time series assumed)
            self._y=[np.reshape(i,(self._T[ix],1),order='F') for ix,i in enumerate(tsData)]
            # dtData = np.diff(tData) #find difference of time series.  within this framework, assume spacing between x_{1|0} and x_1 is dt
            dtData =[np.diff(i) for i in tData]
        else:
            print 'Input to BATCH version incorrect.  Need list of np.ndarrays for BATCH version'
            print 'Aborting initialization (error will result when calling object methods if proceed with input as is)'
            print 'To use this code with a single track encoded in np.ndarray "tsData", just pass in [tsData].'
            return 

        if StaticErrorEstSeq == []:
            self._Rbase = [[0]*i for i in self._T] #initialize with zeros.
            
        elif all(isinstance(x,np.ndarray) for x in StaticErrorEstSeq):
            if all(len(i) == self._T[ix] for ix, i in enumerate(StaticErrorEstSeq)): #check lengths are consistent
                self._Rbase =StaticErrorEstSeq
        else:
            print 'Input to BATCH version incorrect.  Need list of np.ndarrays for BATCH version'
            print 'Aborting initialization (error will result when calling object methods if proceed with input as is)'
            print 'To use this code with a single track encoded in np.ndarray "StaticErrorEstSeq", just pass in [StaticErrorEstSeq].'
            return 


        if dt ==[]:
            # self._dt = np.min(dtData) #find smallest gap in time and assume this is exposure time
            self._dt = [np.min(i) for i in dtData] #find smallest gap in time in each time series and assume this is exposure time
        elif isinstance(dt,float): #allow scalar pass in for exposure time since this will often be constant in across trajectories
            self._dt = [dt]*len(self._T) #make copies to facilitate iterating through object.  recall self_.T is a list of length "len(tsData)" [# batches] containing integer of time series lengths
        elif all(isinstance(x,float) for x in dt) and len(dt)==len(self._T):
            self._dt = dt
        else:
            print 'Option input dt not of correct type or length (input float or list of floats [length latter = number of trajectories]'

        
        
        #self._gapN = [int(np.round(i))-1 for i in dtData/self._dt]  #create an array of frame its.  rounding enforced, if frames are actually non-integer multiples of exposure time, 
        #code is not currently setup for this specialized situation.  if gapN = 0, implies no frame gaps encountered; gapN = 1 one frame gap, etc. 
        #self._gapN = self._gapN + [0] #pad last entry with no gap (not used to evaluate likelihood anyway)
        
        self._gapN = []
        for ix,dtiData in enumerate(dtData):
            tmp = [int(np.round(i))-1 for i in dtiData/self._dt[ix]]
            tmp += [0]
            self._gapN.append(tmp)

    def KFfilterOU1d(self,pars,evalInnovStats='True',P10=[]):
        #run naive simple filter where scalar is observed.  
        #note scalar case doesn't need to worry about symm pos. def. mats and other computational issues commonly encountered in multivariate case  
        #***** make computations easy to follow from a "textbook" classic KF implementation;  
        #***** a module executing the "standard" Kalman filter where correlation between process and measurement noise is zero \forall t is provided to illustrate connection
        #***** between motion blur filter and standard Kalman filter

        loglikelihood = 0.
        xfiltLIST =[] #track the filtered estimates
        pitLIST =[]  #return the pit random variables
        ShistLIST = []
        
        for ix, yi in enumerate(self._y):
            x10 =  pars[3]/pars[0] #for batch, short traj the norm.  smoothing basically meaningless, just use parametric stationary dist mean as prior
            delta = self._dt[ix]
            kappa = pars[0]
            sig   = pars[1]
            Rbasei=self._Rbase[ix]
            gapNi = self._gapN[ix]
            Ti = self._T[ix]
    
          
            
            F,Q,H,R,C,HF,A,HA = self.PttUpdateOU1dPars(pars,delta,0)
            
            if P10 == []:  # assume stationary dist instead of info filter for initial uncertainty if no user supplied uncertainty provided
                if np.abs(kappa)>1E-5:
                    P10 = sig**2/kappa/2. #standard stationary confined diffusion covariance
                else:
                    P10 = np.abs(delta*sig**2 - delta**2*kappa*sig**2) #if small kappa is passed in, use a small t expansion (close to diffusion) to permit 
                    #directed and pure diffusion models without divide by zero issues (taking kappa -> 0 gives exact directed or pure diffusion likelihood, 
                    #but don't trust computational routine to correctly  "infer" limits)
    
            P_Innov = P10
            x00 =  x10 #keep explicit copy of filtered estimate from previous iteration for special blurred filter 
         
     
            xfilt =[] #track the filtered estimates
            pit =[]  #return the pit random variables
            Shist = []
            
            
            #for first observation, assume no gap exists between "x00" and "y1" (i.e., prior state forecast has no gap).  
            #under this assumption initialize filter pars 
            Fi = F;Qi = Q;Hi = H;Ri = R;Ci = C;HFi = HF;Ai = A;HAi = HA
            
    
            for idx, y in enumerate(yi):
     
                Ri = Ri + ( Rbasei[idx] + pars[2] )**2 #use the input static error sequence (squared) to give time dependent R.
                
    
                if evalInnovStats:
                    #compute innovation statistics
                    Sinv = 1./(Hi*P_Innov*Hi+Ri) #use variable P_Innov to stress fundamental difference in innovation computation 
                    #between motion blur and classic kalman filter
                    Sinv = max((Sinv,np.finfo(float).resolution)) #time varying feature or roundoff can introduce negative innovation variance.  only permit values above machine precision
                                                                  #if MLE has Sinv < 0, GoF tests will readily identify this issue. 
                    z = (y-HFi*x00-HAi)*np.sqrt(Sinv) #HF and HA are other nuances of blurred filter formulation used 
                    piti = spspecial.erfc(z/np.sqrt(2.))*.5 #compute CDF of normal
                    # loglikelihood +=  1/2.*(np.log(Sinv)) + -z*z/2. -1/2.*np.log((2*np.pi))
                    loglikelihood += (1/2.*(np.log(Sinv)) + -z*z/2. -1/2.*np.log((2*np.pi)))/Ti
                    
                    pit.append(piti)
                    Shist.append(1./Sinv)
                   
                
                #compute gain and then fuse in information from current measurement 
                K = self.computeGain(P_Innov,Ci,Hi,Ri,Fi) #different gain computation for "classic" KF and MBF
                # K = self.computeGain(P10,C,H,R,F)
                x11 = x10 + K*(y-HFi*x00-HAi) #HF and HA are nuances of blurred filter formulation used 
                x00 = x11 #keep explicit copy of filtered estimate from current iteration for special blurred filter likelihood eval
                xfilt.append(x11) #store filtered estimate
    
                P00 = P10 -  K*(Hi*P_Innov*Hi+Ri)*K 
                #before updating forecast, check to see if gap exists using precomputed info
                if gapNi[idx]>0:
                    Fi,Qi,Hi,Ri,Ci,HFi,Ai,HAi = self.PttUpdateOU1dPars(pars,delta,gapNi[idx]) #compute quantities assuming uniform continuous illumination in next frame
                    #for QD's where blinking can occur within the camera's exposure time, recommended to ignore first time series entry after a blink so above assumption is more reasonable
                else:
                #use precomputed quantities with specified "dt" spacing
                    Fi = F;Qi = Q;Hi = H;Ri = R;Ci = C;HFi = HF;Ai = A;HAi = HA
    
                
                #update/forecast state for simple mean zero OU model
                x10=Fi*x11 + Ai
                P10 = Fi*P00*Fi  +  Qi
    
                P_Innov = self.P_Innov(P10,P00) #nuance of motion blur filter
     
                
            xfilt = np.array(xfilt)
            #store the filtered statistics
            xfiltLIST.append(xfilt)
            pitLIST.append(pit)
            ShistLIST.append(Shist)
            # loglikelihood = loglikelihood/Ti #time average normalization done earlier in BATCH version

        return loglikelihood,xfiltLIST,pitLIST,ShistLIST #return all stats (create wrapper to make more efficient feval calls)
        

    def evalCostFunc(self,pars):
        feval = self.KFfilterOU1d(pars,evalInnovStats='True')
        negloglike = -feval[0][0]
        return negloglike #return negative loglikehood for minimization routines (also set flags to make computation more efficient) 

    def computeGain(self,P10,C,H,R,F):
            K   = (C+F*P10*H)/(H*P10*H+R) #blur form for updating covariance of filtered state.
            return K

    def P_Innov(self,P10,P00): #simple switch function permitting both the classic and motion blur filter with one code base
            return P00 #blur form for updating covariance of filtered state.
            

    def PttUpdateOU1dPars(self,pars,delta,n):
        #par is assumed to contain (kappa,sigma,stdloc,v)
        #delta is the exposure time
        #n is the number of frame gaps (n=0 implies no gaps, n=1 one gaps, etc.)
        kappa = pars[0]
        sigma = pars[1]
        R     =  0 #in time varying code, permit negative parameters reducing input variance and assign 
                   #localization contribution to net measurement noise  in main code (here just assign "blur" contribution)
        
        deltaNet = (n+1)*delta

        if len(pars)>3:
            alpha = pars[3]
        else:
            alpha = 0

        F     = np.exp(-kappa*deltaNet) #standard res
        Q     = (sigma**2/2./kappa)*(1.-np.exp(-2.*kappa*deltaNet))
        fp    = alpha/kappa 
        A     = (1-F)*fp #form assumes kappa>0 implies stable linear system

        Qblur = sigma**2*(2*delta*kappa*np.exp(2*delta*kappa*(n + 1)) - np.exp(2*delta*kappa) + 2*np.exp(delta*kappa) - 2*np.exp(2*delta*kappa*(n + 1)) + 2*np.exp(delta*kappa*(2*n + 1)) - 1)*np.exp(-2*delta*kappa*(n + 1))/(2*delta**2*kappa**3)
        H     = (np.exp(delta*kappa) - 1)*np.exp(-delta*kappa*(n + 1))/(delta*kappa)
        #compute the exact cross correlation term of time integrated OU vs. discretely sampled state (S in notation of Anderson and Moore...I prefer using S for innovation covariance)
        C     = sigma**2*(1 - np.exp(-delta*kappa*(2*n + 1)) + np.exp(-2*delta*kappa*(n + 1)) - np.exp(-delta*kappa))/(2*delta*kappa**2)
        # compute Integral((1-exp(-kappa*(s)))*alpha/delta/kappa,(s,0,delta)) [form also assumes kappa>0 implies stable linear system]
        HA =  fp - fp*np.exp(-delta*kappa*n)/(delta*kappa) + fp*np.exp(-delta*kappa)*np.exp(-delta*kappa*n)/(delta*kappa)
            
 
        
        R     += Qblur #add blur contribution to effective measurement noise 
        HF    = H #special case for blur model
        
        return F,Q,H,R,C,HF,A,HA 


class ClassicKalmanFilter(ModifiedKalmanFilter1DwithCrossCorr):
    """
    generates parameters for using the "blur" version of the 1D KF filter with the "classic Kalman filter" where there is no
    statistical dependence / correlation between process and measurement noise.

    for I/O and notes, see parent class.  the methods redefined here show how to introduce variables and redefine quantities
    to implement the "classic" KF.
    """
    def __init__(self,tsData,tData,StaticErrorEstSeq=[],dt=[]):
        super(ClassicKalmanFilter, self).__init__(tsData,tData,StaticErrorEstSeq,dt)
        
    def computeGain(self,P10,C,H,R,F):
        K   = (P10*H)/(H*P10*H+R) #gain form required for using classic KF within "motion blur filter" formulation. 
        #note: C=0 required for "standard" classic KF (enforced in code)
        return K

    def P_Innov(self,P10,P00): #simple switch function permitting both the classic and motion blur filter with one code base
            return P10 #blur form for updating covariance of filtered state.

    def PttUpdateOU1dPars(self,pars,delta,n):
        #par is assumed to contain (kappa,sigma,stdloc,v)
        #delta is the exposure time
        #n is the number of frame gaps (n=0 implies no gaps, n=1 one gaps, etc.)

        #for classic KF, adjusting for frame gaps is simple since there is no accounting for motion blur in last frame
        delta = (n+1)*delta #only need to mod delta

        kappa = pars[0]
        sig   = pars[1]
        # R     = pars[2]**2 
        R     = 0 #in time varying code, permit negative parameters reducing input variance and assign 
                   #localization contribution to net measurement noise  in main code (here just assign "blur" contribution)
        if len(pars)>3:
            alpha = pars[3]
        else:
            alpha = 0
        #Keep expression below simple.  just note numerical issues may arise if kappa near machine zero is attempted (practically not too relevant since MA1 case and KF should give numerically identical/similar results)
        F     = np.exp(-kappa*delta) 
        Q     = (sig**2/2./kappa)*(1.-np.exp(-2.*kappa*delta))
        H     = 1.
        HF    = H*F
        C     = 0.
        fp    = alpha/kappa 
        A     = (1-F)*fp #assumes kappa>0 implies stable linear system
        HA    = H*A
        return F,Q,H,R,C,HF,A,HA

class ModifiedKalmanFilter1DTimeLapse(ModifiedKalmanFilter1DwithCrossCorr):
    """
        setup special case of "time lapse" movies.  easy to simulate by subsampling blur.  however, frame gaps are the norm here, so modify 
        the code to not recompute upon observation of each gap.  interface changes since exposure time and time spacing between observations are 
        needed in this new framework.
    """
    def __init__(self,tsData,tData,exposureTime,StaticErrorEstSeq=[]):

        super(ModifiedKalmanFilter1DTimeLapse, self).__init__(tsData,tData,StaticErrorEstSeq,exposureTime) #make optional input requirement in new init
        
        # tmp = self._gapN
        # self._lapseSize = max(set(tmp), key=tmp.count) #store the size of the most frequently observed "gap size" (to facilitate computing KF pars).  for uniformly sampled time lapse, only one parameter eval needed if exposure time and time stamp provided
        # self._gapN[-1] = self._lapseSize #replace dummy last arg by most frequently observed value (facilitates interpreting warning messages below)

        #batch version, don't put error checking in.

        self._lapseSize=[]
        for ix,i in enumerate(self._gapN):
            tmp = max(set(i), key= i.count) #store the size of the most frequently observed "gap size" (to facilitate computing KF pars).  for uniformly sampled time lapse, only one parameter eval needed if exposure time and time stamp provided
            self._lapseSize.append(tmp)
            self._gapN[ix][-1] = tmp #replace dummy last arg by most frequently observed value (facilitates interpreting warning messages below)


        for tmp in self._gapN:
            if len(set(tmp))>1:
                print 'Warning: In time lapse info, more than one gap size observed:',set(tmp)
            
        print 'Exposure Time in Time Lapse Object:',exposureTime
        print 'Spacing  of time lapse assumed in computations for (if zero, continuous illumination assumed): ',self._lapseSize

    def KFfilterOU1d(self,parsIn,evalInnovStats='True',P10=[]):
        """
        setup a special interface to make the most frequent gap the default for the parameters of the KF.  
        leverage new info in __init__ to modify parent class code.  changes relatively minor.
        """
        loglikelihood = 0.
        xfiltLIST =[] #track the filtered estimates
        pitLIST =[]  #return the pit random variables
        ShistLIST = []

        pars = langevinScaling(parsIn)
        # pars = parsIn[:]
        
        for ix, yi in enumerate(self._y):
            x10 =  pars[3]/pars[0] #for batch, short traj the norm.  smoothing basically meaningless, just use parametric stationary dist mean as prior
            delta = self._dt[ix]
            kappa = pars[0]
            sig   = pars[1]
            Rbasei=self._Rbase[ix]
            gapNi = self._gapN[ix]
            Ti = self._T[ix]
            lapseSizei = self._lapseSize[ix]
            
        
    
          
            
            F,Q,H,R,C,HF,A,HA = self.PttUpdateOU1dPars(pars,delta,lapseSizei) #compute parameters using the most frequently observed gap size.  should be uniform
            
            
            if P10 == []:  # assume stationary dist instead of info filter for initial uncertainty if no user supplied uncertainty provided
                if np.abs(kappa)>1E-5:
                    P10 = sig**2/kappa/2. #standard stationary confined diffusion covariance
                else:
                    P10 = np.abs(delta*sig**2 - delta**2*kappa*sig**2) #if small kappa is passed in, use a small t expansion (close to diffusion) to permit 
                    #directed and pure diffusion models without divide by zero issues (taking kappa -> 0 gives exact directed or pure diffusion likelihood, 
                    #but don't trust computational routine to correctly  "infer" limits)
    
            P_Innov = P10
            x00 =  x10 #keep explicit copy of filtered estimate from previous iteration for special blurred filter 
         
     
            xfilt =[] #track the filtered estimates
            pit =[]  #return the pit random variables
            Shist = []
            
            
            #initialize filter pars with most common gap size value
            Fi = F;Qi = Q;Hi = H;Ri = R;Ci = C;HFi = HF;Ai = A;HAi = HA
            
    
            for idx, y in enumerate(yi):
     
                Ri = Ri + ( Rbasei[idx] + pars[2] )**2 #use the input static error sequence (squared) to give time dependent R.
                
    
                if evalInnovStats:
                    #compute innovation statistics
                    Sinv = 1./(Hi*P_Innov*Hi+Ri) #use variable P_Innov to stress fundamental difference in innovation computation 
                    #between motion blur and classic kalman filter
                    Sinv = max((Sinv,np.finfo(float).resolution)) #time varying feature or roundoff can introduce negative innovation variance.  only permit values above machine precision
                                                                  #if MLE has Sinv < 0, GoF tests will readily identify this issue. 
                    z = (y-HFi*x00-HAi)*np.sqrt(Sinv) #HF and HA are other nuances of blurred filter formulation used 
                    piti = spspecial.erfc(z/np.sqrt(2.))*.5 #compute CDF of normal
                    # loglikelihood += 1/2.*(np.log(Sinv)) + -z*z/2. -1/2.*np.log((2*np.pi))
                    loglikelihood += (1/2.*(np.log(Sinv)) + -z*z/2. -1/2.*np.log((2*np.pi)))/Ti
                    pit.append(piti)
                    Shist.append(1./Sinv)
                   
                
                #compute gain and then fuse in information from current measurement 
                K = self.computeGain(P_Innov,Ci,Hi,Ri,Fi) #different gain computation for "classic" KF and MBF
                # K = self.computeGain(P10,C,H,R,F)
                x11 = x10 + K*(y-HFi*x00-HAi) #HF and HA are nuances of blurred filter formulation used 
                x00 = x11 #keep explicit copy of filtered estimate from current iteration for special blurred filter likelihood eval
                xfilt.append(x11) #store filtered estimate
    
                P00 = P10 -  K*(Hi*P_Innov*Hi+Ri)*K 
                #before updating forecast, check to see if gap exists using precomputed info
                if gapNi[idx]!=lapseSizei:
                    Fi,Qi,Hi,Ri,Ci,HFi,Ai,HAi = self.PttUpdateOU1dPars(pars,delta,gapNi[idx]) #compute quantities assuming uniform continuous illumination in next frame
                    #for QD's where blinking can occur within the camera's exposure time, recommended to ignore first time series entry after a blink so above assumption is more reasonable.  for time lapse, this isn't an issue
                else:
                #use precomputed quantities with specified "dt" spacing
                    Fi = F;Qi = Q;Hi = H;Ri = R;Ci = C;HFi = HF;Ai = A;HAi = HA
    
                
                #update/forecast state for simple mean zero OU model
                x10=Fi*x11 + Ai
                P10 = Fi*P00*Fi  +  Qi
    
                P_Innov = self.P_Innov(P10,P00) #nuance of motion blur filter
     
                
            xfilt = np.array(xfilt)
            #store the filtered statistics
            xfiltLIST.append(xfilt)
            pitLIST.append(pit)
            ShistLIST.append(Shist)
            
            # xfilt = np.array(xfilt)
            # loglikelihood = loglikelihood/self._T #return empirical time average of loglikelihood

        return loglikelihood,xfiltLIST,pitLIST,ShistLIST  #return all stats (create wrapper to make more efficient feval calls)


def simTimeIntegratedOUdeltaICvel(pars,delta,Tsub,Nsim,T,xIC=0):
    """
    inputs:

    par   : assumed to contain (kappa,sigma,stdloc,vel)
    delta : time interval between observations
    Tsub  : integer specifying granularity of division between delta (used for simple numerical time integration)
    Nsim  : number of paths to simulate
    T     : length of each path
    xIC   : Dirac point mass of ensemble.  code can be modified to handle distribution of ICs., though not relevant to this study.

    outputs:

    xraw  : TxN nparray of (exact) discretely sampled OU process 
    xblur : "     "   simple time quadrature of xraw with accuracy dictated by Tsub (avoid using analytical expression in order to show convergence to numerical results)

    [both outputs check out when compared to analytical expressions]    
    """
    

    kappa = pars[0]
    sig   = pars[1]
    R     = pars[2]**2
    if len(pars)>3:  #allow for varying length inputs (assume constant vel is zero typically)
        vel   = pars[3]
    else:
        vel   = 0.
 
    dt    = float(delta)/float(Tsub)
    F     = np.exp(-kappa*dt)

    if np.abs(kappa)>1E-5:    #this code chunk allows an SPT legacy model "pure directed diffusion" (if kappa^2 is near double precision zero, no computational difference)
        fp = vel/kappa 
        A  = (1-F)*fp #assumes kappa>0 implies stable linear system
        sqrtQ = np.sqrt((sig**2/2./kappa)*(1.-np.exp(-2.*kappa*dt)))
    else:
        A  = vel*(dt - kappa*dt**2/2.) # + O(delta^3) + O(kappa^2)
        sqrtQ = np.sqrt(dt*sig**2 - dt**2*kappa*sig**2) # + O(delta^3)++ O(kappa^2)

    

    
    x0    = np.ones(Nsim)*xIC
    
    xraw=[]
    xblur=[]
    xrawsampled=[]

    for i in range(T):
        xloc=[]
        for j in range(Tsub):
            noise = sqrtQ*np.random.randn(Nsim)
            x0=F*x0+noise+A
            xloc.append(x0)
            xraw.append(x0)

        xrawsampled.append(x0)
        xblur.append(np.mean(xloc,axis=0))

    xraw = np.array(xraw)
    xrawsampled = np.array(xrawsampled)
    xblur= np.array(xblur)
    

    return xraw,xrawsampled,xblur

def simTimeIntegratedOUdeltaICvel_RBM(pars,delta,Tsub,Nsim,T,xIC=0):
    """
    inputs:

    par   : assumed to contain (kappa,sigma,stdloc,vel)
    delta : time interval between observations
    Tsub  : integer specifying granularity of division between delta (used for simple numerical time integration)
    Nsim  : number of paths to simulate
    T     : length of each path
    xIC   : dirac point mass of ensemble.  code can be modified to handle distribution of ICs., though not relevant to this study.

    outputs:

    xraw  : TxN nparray of (exact) discretely sampled OU process 
    xblur : "     "   simple time quadrature of xraw with accuracy dictated by Tsub (avoid using analytical expression in order to show convergence to numerical results)

    [both outputs check out when compared to analytical expressions]    
    """
    

    kappa = pars[0]
    sig   = pars[1]
    R     = pars[2]**2
    
    L=np.sqrt(6.*sig**2/kappa)
    
    if len(pars)>3:  #allow for varying length inputs (assume constant vel is zero typically)
        vel   = pars[3]
    else:
        vel   = 0.
 
    dt    = float(delta)/float(Tsub)
    sqrtQ = sig*np.sqrt(dt) #RBM specific
#     A = vel #RBM specific
    fp = vel/kappa
    fp += -L/2. + xIC #shift to have box center at same location as corresponding OU process mean and have prescribed IC

    
    x0    = np.ones(Nsim)*xIC
    W     = np.zeros(Nsim)+L/2.
    
    xraw=[]
    xblur=[]
    xrawsampled=[]

    for i in range(T):
        xloc=[]
        for j in range(Tsub):
            noise = sqrtQ*np.random.randn(Nsim)
            W = W + noise
            
            
            #use a simple upcrossing scheme to impose barrier condition (not valid  if A!=0, so just enforce A=0)
            for ix in range(Nsim):
                if np.mod(np.floor(W[ix]/L),2)==0: #even case
                    x0i = W[ix]-np.floor(W[ix]/L)*L
                else: #odd case
                    x0i = L-(W[ix]-np.floor(W[ix]/L)*L)    
                x0[ix] = x0i
#             x0 = np.random.randn(Nsim) #DEBUG looking for strange mutable list "append" side effect.  if one accumulate nparrays
            #in "vectorized storage, things work as expected, but if x0 as a standard list is appended, mutable side effects cause problems.
            #why does np.ndarray guard against this silliness?  seems assigning a element pointing to a scalar alone is BAD (if scalar changes, all copies update
            #wheras if full array appended, "copy" implicit...very strange ) 
            
                
#             xloc.append(x0) 
#             xraw.append(x0)
            xloc.append(x0.copy()) 
            xraw.append(x0.copy())
        
#         xrawsampled.append(x0)
        xrawsampled.append(x0.copy()) #weird, without copy get strange mutable side effects...  document here. 
        
        xblur.append(np.mean(xloc,axis=0))

    xraw = np.array(xraw)  + fp
    xrawsampled = np.array(xrawsampled) + fp  
    xblur= np.array(xblur) + fp  
    

    return xraw,xrawsampled,xblur  








