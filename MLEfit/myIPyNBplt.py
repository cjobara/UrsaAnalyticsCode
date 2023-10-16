
# Written by Chris Calderon 2014 (Chris.Calderon@UrsaAnalytics.com)
#
#
# Copyright 2014 Ursa Analytics, Inc.

   # Licensed under the Apache License, Version 2.0 (the "License");
   # you may not use this file except in compliance with the License.
   # You may obtain a copy of the License at

   #     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
myIPyNBplt:

Module containing some functions allowing flexible matplotlib plotting in IPyNB [hides some messy details in matplotlib func calls]
Also contains some I/O routines commonly used when messing with data.

"""

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import h5py
import numpy as np
import scipy as sp

def basicPltSetup(axisfontsize='22',legendfontsize='18',pltFont='Times New Roman'):
    """
    Simple wrapper for using custom font (Helvetica) and desired axis tick / font sizes (returns a dictionary)

    For distribution, probably want to remove fontname since users might have a hard time installing custom fonts 
    """
    pltPars={}
    axis_font = {'fontname':pltFont, 'size':axisfontsize}
    title_font = {'fontname':pltFont, 'size':axisfontsize, 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    
    # font_path = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/' + pltFont +'.ttf'
    # font_prop = font_manager.FontProperties(fname=font_path, size=axisfontsize) #used in custom legend below 
    # legendfont_prop = font_manager.FontProperties(fname=font_path, size=legendfontsize) #used in custom legend below
    
    fig = plt.figure() 
    ax  = fig.add_subplot(111,adjustable='box')

    #pack things up for revising plot later in dict (see if possible to adjust before plot made)
    pltPars['axis_font'] = axis_font
    # pltPars['title_font'] = title_font
    # pltPars['font_prop'] = font_prop
    # pltPars['leg_prop'] = legendfont_prop #e.g. of use: lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,ncol=1,numpoints=1,prop=pltPars['leg_prop'])
    pltPars['axisfontsize'] = axisfontsize

    return fig,ax,pltPars

def sharedXPltSetup(axisfontsize='22',legendfontsize='18',pltFont='Helvetica'):
    """
    Simple wrapper for setting up shared x-axis with custom font (Helvetica) and desired axis tick / font sizes (returns a dictionary)
    After adding data to plots, likely need to rerun
    fig.subplots_adjust(hspace=0) 
    in ipythonNB


    For distribution, probably want to remove fontname since users might have a hard time installing custom fonts 
    """
    
    pltPars={}
    axis_font = {'fontname':pltFont, 'size':axisfontsize}
    title_font = {'fontname':pltFont, 'size':axisfontsize, 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    
    font_path = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/' + pltFont +'.ttf'
    # font_prop = font_manager.FontProperties(fname=font_path, size=axisfontsize) #used in custom legend below 
    # legendfont_prop = font_manager.FontProperties(fname=font_path, size=legendfontsize) #used in custom legend below

    fig = plt.figure() 
    fig, (axTop, ax) = plt.subplots(2, sharex=True)
    fig.subplots_adjust(hspace=0) #cannot think of situation where this isn't desirable in shared x-axis setting

    #pack things up for revising plot later in dict (see if possible to adjust before plot made)
    pltPars['axis_font'] = axis_font
    pltPars['title_font'] = title_font
    # pltPars['font_prop'] = font_prop
    # pltPars['leg_prop'] = legendfont_prop #e.g. of use: lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,ncol=1,numpoints=1,prop=pltPars['leg_prop'])
    pltPars['axisfontsize'] = axisfontsize

    return fig,ax,axTop,pltPars

def xyLabels(ax,pltPars,xlabel='',ylabel='',title='',pltFont='Helvetica'): #don't use Helvetica in deployed versions (likely a pain to setup)
    """
    set the axis number font size and type.  option, put  x and y labels [for pubs, leave empty and use latex overlay ] 
    (requires pltPar dict from basicPltSetup) 
    """
    ax.set_xlabel(xlabel, **pltPars['axis_font'])
    ax.set_ylabel(ylabel, **pltPars['axis_font'])
    ax.set_title(title, **pltPars['axis_font'])
    
    #set font size
    for item in ([ax.xaxis.label, ax.yaxis.label ] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(pltPars['axisfontsize'])
    
    # Set the tick font type
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(pltFont)


def wouth5(floc,DATAMAT): #write a simple hdf5 file into "floc";  assumes numpy array passed as DATAMAT.  handling other datatypes with h5py is fairly easy
    """
    write an hdf5 file dumping everything into generic name useful for playing with other programs "/dataset0."  

    floc: file location  (can use "~" in path) 
    DATAMAT: numpy matrix or vector


    """
    if floc[0]=='~':  #replace standard unix home shortcut with explicit path
                floc='/Users/calderoc' + floc[1:] #purdue login: calderoc / princeton: ccaldero; since grad school i use either calderoc or ccalderoN as unix logins (keep symlink to both home folders in *nix style systems)
    f=h5py.File(floc,'w')#default name for h5import dataset when no arg given (so i use this in my scripts as well)
    dset = f.create_dataset('dataset0', data=DATAMAT)
    f.close()
    return 0

#to get above into R use:
#>> library(hdf5);hdf5load(floc);DATAMAT=dataset0
#to get above into MATLAB, use
#>> DATAMAT=h5read(floc); %this is a simple wrapper script i wrote in ~/mfiles


def loadh5(floc,dsetname='/dataset0'): #load a numeric hdf5 file. "floc" is the path of ASCII file to read.
    """
    read an hdf5 file assuming everything of interest is stored in dsetname "/dataset0."  h5py can be leveraged for more complex HDF5 storage schemes
  
    """

    if floc[0]=='~':  #replace standard unix home shortcut with explicit path
        floc='/Users/calderoc' + floc[1:]
    f=h5py.File(floc,'r')
    # print 'Read file: ' + floc
    # print 'Pulling Dataset Named: ' + dsetname
    g=f[dsetname]  #default name for h5import dataset when no arg given (so i use this in my scripts as well)
    #f.close()  #closing before return causes problems.  look into possible issue associated with not formally closing a file upon exit of a function 
    #apparently after you go out of scope, "auto close" occurs avoiding memory leaks.
    return np.array(g)  #return numpy array (get shape with g.shape)

def wout(floc,DATAMAT): #write out an ASCII file to "floc" given a 2D numpy array "DATAMAT"
    """
    write an ASCII file.  

    floc: file location  (can use "~" in path) 
    DATAMAT: numpy matrix or vector
    """
    if floc[0]=='~':  #replace standard unix home shortcut with explicit path
                floc='/Users/calderoc' + floc[1:]
    np.savetxt(floc, DATAMAT)
    return 0



