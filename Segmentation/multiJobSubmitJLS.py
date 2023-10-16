#!/usr/bin/python

##################################################################################################
#Simple python batch script for running modified HDP-SLDS segmentation 
#Copyright 2017 Ursa Analytics, Inc. All Rights Reserved. 
# 
#
# Disclaimer of Warranty. Unless required by applicable law or agreed to in writing, Licensor provides 
# the Work (and each Contributor provides its Contributions) on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, 
# without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, 
# or FITNESS FOR A PARTICULAR PURPOSE. You are solely responsible for determining the appropriateness 
# of using or redistributing the Work and assume any risks associated with Your exercise of permissions 
# under this License. 

   # Licensed under the Apache License, Version 2.0 (the "License");
   # you may not use this file except in compliance with the License.
   # You may obtain a copy of the License at

   #     http://www.apache.org/licenses/LICENSE-2.0
   
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



# sample calls at bash terminal
#>> ./PATHOFTHISFILE/multiJobSubmit.py
#>> python multiJobSubmit.py
#>> nohup ./multiJobSubmit.py > LABEL.out 2> LABEL.err < /dev/null &


##################################################################################################


import subprocess, os, time, math

#net data size.  to do: write interface python program for translating csv to matlab info or to intermediate hdf5 (then write simple matlab routine for reading in hdf5 files).  latter option useful in plotting in python later
#Ntraj=25
#
nproc=4 #specify the number of processors free to use in computational environment (recommended to set <= to number of cores available)

fcase='2017.03.06_selected'
outFileNameBase='./TestBatch_' + fcase  # absolute or relative path (former may be preferred in some aps)
readFolderNameBase='../JLSdataDec20/' + fcase  # absolute or relative path (former may be preferred in some apps)
os.system('mkdir -p ' + outFileNameBase) #this should also be in the code, but just make at os level in case MATLAB mkdir fails
outFileNameMATLAB='\'' + outFileNameBase + '\''  #add quotes via escape characters (works in python 2) to string for MATLAB's string parsing.  
readFolderNameMATLAB='\'' + readFolderNameBase + '\'' 

p=subprocess.Popen('ls ' + readFolderNameBase + '/*.csv',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
stdout,stderr=p.communicate()
Ntraj = len(stdout.split())
print(stdout.split())
print('Number of  csv files found in folder: ',Ntraj)

#CRUDE division of labor scheme
dirlist=range(1,Ntraj+1) #i still prefer indices starting at 1...
totalJobs=len(dirlist)
jobsPerRun=int(math.ceil(totalJobs/nproc)) #ceiling ensures all cases hit
caseind2=range(0,totalJobs+jobsPerRun,jobsPerRun) #adds padding so last job isn't empty
# caseind.pop() #remove the last spurious item

def divyup(dirlist, nproc):
    """
    take N data files and divide up in "nproc" sublists based on input list specified in "dirlist".
    """
    # myinput = range(Ntraj)
    myinput = dirlist #just pass list in directly vs. assuming integer indexing
    input_size = len(myinput)
    slice_size = input_size / nproc
    remain = input_size % nproc
    result = []
    iterator = iter(myinput)
    for i in range(nproc):
        result.append([])
        for j in range(slice_size):
            result[i].append(iterator.next())
        if remain:
            result[i].append(iterator.next())
            remain -= 1
    return result

caseind = divyup(dirlist,nproc)

print('Division of labor across N=',nproc, 'processors')
print('*'*72)
for ix,i in enumerate(caseind):
    print('Range assigned to processor ',ix, ': ', i[0],i[-1])
    print('Load: ',i[-1] - i[0])

print

SOURCEPATH="./"
#setup file location info
matlabpath='cd ' + SOURCEPATH + ';' #this is the location of the matlab batch script called in "basecmd" [the codes eventually leaves this main dir and goes into the specified data dir]
basecmd='nohup matlab -nosplash -nodisplay -r '  #provide options for matlab command line call.

#NOTE on linking dependencies like lightspeed.  code links to absolute paths of 3rd party libraries in the install folder.  do not move these files otherwise 
#script will crash.



#setup matlab command passed in terminal;  test first by running script with code below "EXE PORTION" commented out
matlabarglist='\"'+ 'runstuffUrsaBatch(99,[${INDLOW}:${INDHI}],' + outFileNameMATLAB +',' + readFolderNameMATLAB + ');\"'  #use of quotes and escape characters needed in MATLAB 2015.

nohupcmd='> ${INDLOW}.out 2> ${INDLOW}.err < /dev/null &'
for i in caseind:
# for i in caseind[10:-10]:  #if re-running subset

    ub=i[-1]
    lb=i[0]
    print lb,ub

    os.environ["INDLOW"]=str(int(lb)) 
    os.environ["INDHI"]=str(int(ub))
    print matlabarglist
    cmd=matlabpath + basecmd + matlabarglist+ nohupcmd
    print('script set to run following at terminal from CWD within subprocess:')
    print(cmd)
    
    #EXE PORTION############################################
    #os.system(cmd)
    #print('testing ENV variable set')
    #p=subprocess.Popen('echo $INDLOW',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    #stdout,stderr=p.communicate()
    #print(stdout)

    print('Submitting jobs (check stderr output below...empty if worked)')
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    stdout,stderr=p.communicate()
    print(stderr)


