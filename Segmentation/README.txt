========================================================================
Copyright & License
========================================================================
Copyright (C) 2016, Ursa Analytics, Inc.

Permission is granted for anyone to copy, use, or modify these
programs and accompanying documents for purposes of research or
education, provided this copyright notice is retained, and note is
made of any changes that have been made.

These programs and documents are distributed without any warranty,
express or implied.  As the programs were written for research
purposes only, they have not been tested to the degree that would be
advisable in any important application.  All use of these programs is
entirely at the user's own risk.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

The pfiles in this directory contain modifications made in 2016 of open-source code downloaded from:
https://www.stat.washington.edu/~ebfox/software/HDPHMM_HDPSLDS_toolbox.zip.

As required by the derivative work license, the original license is retained and the changes
are noted below:

* added functionality to map continuous time SDE parameter priors to discrete state space models enableing exact likelihoods for Kalman Filter (no approximations)

* updated corresponding inverse wishart function calls

* updated label assignment scripts to save and plot trajectory segmentation results based on new MCMC output

* added new driver to account for single particle tracking models driven by standard diffusion (classic free, confined, and directed diffusion)

* added new option to runstuff.m in source code in "runstuffUrsaBatch.p".  this allows for "inverse Wishart on (A,Sigma) and normal on mu ('MNIW-N'), matrix normal" option for  AR(1) KF model (with random level offset) using continuous time SDE parameters driving priiors discussed in Refs. below.

* added python batch run routine call MATLAB routines (tested on R2015b)  with wrapper script "multiJobSubmitJLS.py".  this code depends on MATLAB and lightspeed addon toolbox.  data in paper was produced with lightspeed compiled for OS X 10.12 with MATLAB  R2015b.

Technical details provided in: 

Calderon, C. P. (2014). Data-Driven Techniques for Detecting Dynamical State Changes in Noisily Measured 3D Single-Molecule Trajectories. Molecules, 19(11), 18381–18398. http://doi.org/10.3390/molecules191118381

Calderon, C. P., & Bloom, K. (2015). Inferring Latent States and Refining Force Estimates via Hierarchical Dirichlet Process Modeling in Single Particle Tracking Experiments. PloS ONE, 10(9), e0137633. http://doi.org/10.1371/journal.pone.0137633

Derivative Works and Associated Copyright / License:
========================================================================
Copyright & License
========================================================================

Copyright (C) 2009, Emily B. Fox and Erik B. Sudderth.

http://web.mit.edu/ebfox/www/

Permission is granted for anyone to copy, use, or modify these
programs and accompanying documents for purposes of research or
education, provided this copyright notice is retained, and note is
made of any changes that have been made.

These programs and documents are distributed without any warranty,
express or implied.  As the programs were written for research
purposes only, they have not been tested to the degree that would be
advisable in any important application.  All use of these programs is
entirely at the user's own risk.

========================================================================
Acknowledgments
========================================================================

Portions of the package were adapted from Yee Whye Teh's
"Nonparametric Bayesian Mixture Models" package, release 1.
Available from:  http://www.gatsby.ucl.ac.uk/~ywteh
