#!/usr/bin/env python

#collection of routines for estimating diffusion + noise model from uniformly sampled observations.
#modification to non-uniform sampling straight-forward (just modify diagonal...this makes eigen-values harder to compute,
    # but efficient sparse solvers can be leveraged)...code just uses "ideal" time spacing as a proxy.  KF and MBF code can properly handle non-uniform observations.

#this versions allows for arbitrary correlation in the diagonal term by using a nonstandard parameterization of 
#the MA1 (introduced by Berglund in 2010 PRE).  only relevant if blur dominates all "measurement noise" 

import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as sla
import scipy.sparse as sparse
import timeit
import h5py

class CostFuncMA1Diff(object):
    def __init__(self,tsData,dt,blurCoef=1./6.): 
        """
        class for setting up costfunc of MA1 of differenced measurements in diffusion plus noise model
        dx_t= v dt+ sqrt(sqrtD*2)dBt
        y_ti = x_ti + \epsilon_i*sqrtR
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be...see comments for code changes required)
        """
        

        self._T=max(tsData.shape)-1 #compute length of differenced time series
        tsData=np.reshape(tsData,(self._T+1,1),order='F')#reshape data to column vector (scalar time series assumed)
        self._dy=np.diff(tsData,axis=0) #compute differenced time series (do not manipulate this, modify local copies via -pars[2]*dt)
        
        ii=np.arange(self._T) #construct and store coo sparse indices for MA1
        self._rowi=np.concatenate((ii,ii[:-1],ii[:-1]+1))
        self._coli=np.concatenate((ii,ii[:-1]+1,ii[:-1]))
        self._dt=dt #store internally to simplify interface to 3rd party opt functions
        self._blurCoef=blurCoef

    def evalCostFuncVel(self,pars):
        """
        interface method to par opt routine
        MA1 of differenced measurements in diffusion plus noise model
        dx_t= v dt+ sqrt(sqrtD*2)dBt
        y_ti = x_ti + \epsilon_i*sqrtR
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be...see comments for code changes required)
        """
        
        #use Berglund parameterization (allows for negative or positive correlation in MA1 covariance)
        #np.array([1.]*self._T)*2.*pars[0]**2*self._dt
        Reff   = pars[1]**2 - 2.*pars[0]**2*self._dt*self._blurCoef #use equation 7 from Berglund for R = blurCoef (permits negative evals) 
        Deff   = pars[0]**2*self._dt

        #compute MA1 covariance matrix given effective parameters above
        T = len(self._dy)
        S = self.constructMA1(Deff,Reff,Tsize=T)


        valsE  = np.array([1.]*T)*2.*Deff #formulation assume covariance of form: 2*Deff*Id + Reff*SecondDiffMat
        tmp    = np.arange(T)+1
        valsE2 = (4.*Reff)*  ((np.sin(tmp*np.pi/2./(T+1)))**2) #pure dirichlet boundary conditions
        vals   = valsE + valsE2 #
        #below sets up eigenvectors 
        # T=self._T;V=np.array([[np.sqrt(2/(T+1.))*np.sin(i*j*np.pi/(T+1.)) for j in np.arange(T)+1] for i in np.arange(T)+1])
        # #norm of V*diag(vals)*V.T - full(S) should be zero (e.g., eigen decomposition)
        # S2 = np.dot(V.T,np.diag(vals))
        # S2 = np.dot(S2,V)
        ############################################
        loglikelihood=sum(np.log(vals))/2.
        
        #compute quadratic form contribution to log likelihood
        dy=self._dy-pars[2]*self._dt #
        #execute solve required to compute quadratic form
        tmp = self.sparseSolveWrapper(S,dy)
        quadForm = np.dot(dy.T,tmp)

        loglikelihood+=quadForm/2. #TODO:  replace this line by adding in quadratic form.
        #note negative of (unormalized) loglikelihood  computed above 


        return loglikelihood

    def evalCostFunc(self,pars):
        """
        interface method to par opt routine
        MA1 of differenced measurements in diffusion plus noise model
        dx_t= 0 dt+ sqrt(sqrtD*2)dBt
        y_ti = x_ti + \epsilon_i*sqrtR
        pars:=(sqrtD,sqrtR)
        dt:= scalar giving time sampling (can modify to vector dt if need be...see comments for code changes required)
        """
        
   

   
        #use Berglund parameterization (allows for negative or positive correlation in MA1 covariance)
        #np.array([1.]*self._T)*2.*pars[0]**2*self._dt
        Reff   = pars[1]**2 - 2.*pars[0]**2*self._dt*1./6. #use equation 7 from Berglund for R = 1/6. (permits negative evals) 
        Deff   = pars[0]**2*self._dt

        #compute MA1 covariance matrix given effective parameters above
        #S = self.constructMA1(Deff,Reff,Tsize=len(self._dy))
        T = len(self._dy)
        S = self.constructMA1(Deff,Reff,Tsize=T)

        valsE  = np.array([1.]*T)*2.*Deff #formulation assume covariance of form: 2*Deff*Id + Reff*SecondDiffMat
        tmp    = np.arange(T)+1
        valsE2 = (4.*Reff)*  ((np.sin(tmp*np.pi/2./(T+1)))**2) #pure dirichlet boundary conditions
        vals   = valsE + valsE2 #
        #below sets up eigenvectors 
        # T=self._T;V=np.array([[np.sqrt(2/(T+1.))*np.sin(i*j*np.pi/(T+1.)) for j in np.arange(T)+1] for i in np.arange(T)+1])
        # #norm of V*diag(vals)*V.T - full(S) should be zero (e.g., eigen decomposition)
        # S2 = np.dot(V.T,np.diag(vals))
        # S2 = np.dot(S2,V)
        ############################################
        loglikelihood=sum(np.log(vals))/2.
        
        #compute quadratic form contribution to log likelihood
        # dy=self._dy-pars[2]*self._dt #
        dy=self._dy-0.*self._dt #
        #compute solve required to compute quadratic form
        tmp = self.sparseSolveWrapper(S,dy)
        quadForm = np.dot(dy.T,tmp)

        loglikelihood+=quadForm/2. #TODO:  replace this line by adding in quadratic form.
        #note negative of (unormalized) loglikelihood  computed above 


        return loglikelihood
        
    def sparseEigs(self,S):
        """
        compute eigenspectrum in parts for sparse SPD S of size nxn.  sparse symmetric eigen problem should be doable in one quick shot, but not currently possible in scipy.sparse.linalg    
        use krylov based eigensolver here to get full spectrum in two phases (built-in scipy funcs won't return full eigenspectrum)

        this routine is only needed for nonuniform time spacing case
        """ 
        k1 = int(np.ceil(self._T/2.))
        vals1 = sla.eigsh(S, k=k1,return_eigenvectors=False,which='LM') 
        k2 = int(np.floor(self._T/2.))
        vals2 = sla.eigsh(S, k=k2,return_eigenvectors=False,which='SM')
        vals=np.concatenate((vals1,vals2))

        return vals
    
    def constructMA1(self,Deff,Reff,Tsize=None):
        """
        precompute the coo sparse matrix indices of a tri-banded MA1  matrix (stored in rowi, coli) and return sparse coo mat
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be)
        Reff:=  effective measurement noise variance
        Deff:= effective measurement noise (blur adjustments made elsewhere)
        Tsize:= length of data set (allows for different size trajectories in a batch...check for this condition in calling routine and pass in size if need be)
        """
        #form sparse MA1 matrix
        if Tsize == None:
            Ti = self._T #use precomputed size of first traj
            rowi = self._rowi
            coli = self._coli

        else:
            Ti = Tsize 
            ii=np.arange(Ti) #construct and store coo sparse indices for MA1
            rowi=np.concatenate((ii,ii[:-1],ii[:-1]+1))
            coli=np.concatenate((ii,ii[:-1]+1,ii[:-1]))

        R=Reff
        mainDiag=(2*R+2*Deff)*(np.array([1.]*Ti)) #expression "np.array([1.]*N" like matlab ones(N,1) (with row/column left open)
        band=-R*(np.array([1.]*(Ti-1)))
    
        svals=np.concatenate((mainDiag,band,band))
        svals=np.array([float(i) for i in svals]) #crude approach to computing a array with shape (T,) vs (T,1).  difference required for sparse
        S=sparse.coo_matrix((svals,(rowi,coli)),shape=(Ti,Ti))
        return S
    
    def sparseSolveWrapper(self,S,RHS):
    ###############Solver 1 [Mystery Method...that is user is just accepting whatever scipy's soup de jour is...better to be explicit on solver]##################### 
        # start_time = timeit.default_timer()
        # tmp = sla.spsolve(S.tocsc(),RHS) #this doesn't give details on how solve is done.  should exploit SPD nature of S for solve...unsure if scipy does this by default.  
        # elapsed = timeit.default_timer() - start_time
        # print 'tmystery', elapsed
        # print tmp
    ###############
        
    
    ############### Solver 2 [Explicit Call to SupaLU] ##################### 
        # start_time = timeit.default_timer()
        supalu = sla.splu(S.tocsc())
        tmp = supalu.solve(RHS.reshape(-1))
        # elapsed = timeit.default_timer() - start_time
        # print 'tsuperLU',elapsed 
        # print tmp
    ###############
        return tmp

class CostFuncMA1Diff_Batch(CostFuncMA1Diff):

    def __init__(self,tsData,dt,blurCoef=1./6.):
        if isinstance(tsData,list):
        #initialize 
            print 'Note: This version of MA code assumes uniform time spacing (for nonuniform sampling result is known to be misspecified)'
            super(CostFuncMA1Diff_Batch,self).__init__(tsData[0],dt,blurCoef=blurCoef) #make optional input requirement in new init
            self._tsData=tsData #store list for later use
            print 'Using MA Blur Coefficient: ',blurCoef
        else:
            print 'Error in object init.  This version requires list of time series for input.'

    def evalCostFuncVel_Batch(self,pars):
        
        Nbatch = len(self._tsData)
        loglikelihood = 0
        for tsi in self._tsData:
            self._dy=np.diff(tsi,axis=0) 
            loglikelihoodi = self.evalCostFuncVel(pars)
            loglikelihood += loglikelihoodi/Nbatch
        return loglikelihood

def wouth5(floc,DATAMAT): #write a simple hdf5 file into "floc";  assumes numpy array passed as DATAMAT.  handling other datatypes with h5py is fairly easy
    """
    simple wrapper utility func for  
    writing an hdf5 file dumping everything into generic name useful for playing with other programs "/dataset0."  

    floc: name of output file (path) [can use "~" in path] 
    DATAMAT: numpy matrix or vector



    """
    if floc[0]=='~':  #replace standard unix home shortcut with explicit path
                floc='/Users/calderoc' + floc[1:]
    f=h5py.File(floc,'w')#default name for h5import dataset when no arg given (so i use this in my scripts as well)
    dset = f.create_dataset('dataset0', data=DATAMAT)
    f.close()
    return 0

def runOPT(sig1=.4,sig2=.4,sqrtR=35/100.,dt=10/1000.,N=10,T=50):
    #>> findMA1.runOPT(sig1=.4,sig2=.4,sqrtR=35/100.,dt=10/1000.,N=10,T=50)
    ts=np.arange(N)

    dW,W=simSDE.simDB(T=T,N=N,dt=dt,sig=sig1*np.sqrt(2))
    dW2,W2=simSDE.simDB(T=T,N=N,dt=dt,sig=sig2*np.sqrt(2))
    

    W2+=np.reshape(W[:,-1],(N,1)) #create a smooth transition by adding terminal value of 
    print W.shape
    W=np.hstack((W,W2))
    print W.shape

    Y=W+np.random.randn(W.shape[0],W.shape[1])*sqrtR
    fracsplit=.5 #adjust this parameter to reflect mix of sig1 and sig2 in sampled data
    sigEff = np.sqrt(fracsplit*2*sig1**2+fracsplit*2*sig2**2) 
    sigEff = sigEff/2. #make sigEff^2= D_Eff 
    Xtrue= np.array([sigEff,sqrtR,0])
    
    #iterate over paths and carry out optimization
    resH=[]
    for i,yi in enumerate(Y):
        print 'Iteration:',i
        costInstance = CostFuncMA1Diff(yi,dt)
        res = spo.minimize(costInstance.evalCostFunc, Xtrue/2., method='nelder-mead',options={'xtol': 1e-5, 'disp': False})
        print res.x
        resH.append(res.x)
    # res = spo.minimize(test2.evalCostFunc, Xtrue, method='nelder-mead',options={'xtol': 1e-4, 'disp': True})
    resH=np.asarray(resH)

    print '******* Result Summary ***********************'
    print ''
    print 'True (or Effective) Par:', Xtrue
    print ''
    print 'parameter means,medians, max, min of NxPar history vec:'
    print np.mean(np.abs(resH),axis=0) #takes abs value since optimization was unconstrained (cost function squares sig and sqrtR, so no diff;  physically both pars must be >0)
    print np.median(np.abs(resH),axis=0)
    print np.max(np.abs(resH),axis=0)
    print np.min(np.abs(resH),axis=0)
    print 'parameter STD of NxPar history vec:'
    print np.std(np.abs(resH),axis=0)

    print 'Ddt/2R:', sig1**2*dt/(2.*sqrtR**2)
    fname = '~/optRes_dt_'+str(int(dt*1000)) + '_sig1_' + str(sig1) +  '_sig2_' + str(sig1) +'.h5'
    print "saving  results to file: ", fname
    wouth5(fname,resH)

def simDB(N=1000,T=20,dt=2./100.,sig=1): #write out an ASCII file to "floc" given a 2D numpy array "DATAMAT"
    """
      simulate batch of increments and paths of pure BM.  
    """
    sdt=np.sqrt(dt)
    dW=np.random.randn(N,T)*sdt*sig;
    W=np.cumsum(dW,axis=1)
    return dW,W


 
