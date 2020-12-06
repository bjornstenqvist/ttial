#!/usr/bin/env python
# coding: utf-8

# ### A notebook for running and visualizing results from
# # The Tea Is A Lie
# ### - TTIAL: a Mesh Analysis Software for Steady State Diffusion 
# 
# #### available at Zenodo
# #### DOI: XYZ
# 
# By Björn Stenqvist, Div. of Physical Chemistry, Lund University, Sweden
#  
# The following Python packages are needed to run the notebook:<br> numpy<br> matplotlib

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np, os
import matplotlib, matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl

fontSizeLabel=10
fontSizeLegend=8
fontSizeAll=8
plt.rcParams.update({'font.size': fontSizeAll, 'figure.figsize': [3.513475, 2.713475],'xtick.labelsize':fontSizeAll,'ytick.labelsize':fontSizeAll})
rc('font',**{'family':'serif','sans-serif':['Times']}) # set font
rc('text', usetex=True) # use tex
color_c='black' # color for reference
color_mu='black' # color for this work
insetColor='black' # color of lines in insets
cmaps = 'plasma'
R=8.31446261815324 # gas constant
T=298.15 # temperature
RT=R*T
mum=1e-6 # micrometer

def findValue(filename,key):
    with open (filename, "r") as hfile:
      sp = hfile.read()

    lines = sp.split("\n")
    for line in lines:
        words = line.split(" ")
        if words[0] == key:
            return words[1]
        
def getDimensions(filename):
    height = float(findValue(filename,'height'))
    width = float(findValue(filename,'width'))
    return height, width


# ## Generate input

# In[ ]:


get_ipython().run_cell_magic('writefile', 'input.txt', 'model_nbr     0     # model number, 0 gives Brick and Mortar\nd             5.0   # width of bricks\ns             0.2   # horizontal spacing between bricks\ng             0.2   # vertical spacing between bricks\nt             1.0   # brick thickness\nN             5     # nbr of layers of bricks\nomega         1.0   # offset ratio, negative gives random\nc_out         1.0   # concentration on one side of system (zero on the other)\nS_out         1.0   # solubility outside system\nS_mv          1e1   # mortar solubility (verticle)\nS_mh          1e1   # mortar solubility (horizontal)\nS_bv          1e-1  # brick solubility (verticle)\nS_bh          1e-1  # brick solubility (horizontal)\nD_mv          1.0   # mortar diffusion coefficient (verticle)\nD_mh          1.0   # mortar diffusion coefficient (horizontal)\nD_bv          1.0   # brick diffusion coefficient (verticle)\nD_bh          1.0   # brick diffusion coefficient (horizontal)\nseed          1     # seed, negative gives random\nheight        5.4   # height of system, only relevant when loading an external mesh\nwidth         10.1  # width of system, only relevant when loading an external mesh\nload_external false # load external mesh, true/false\nNc            26    # number of columns of nodes in mesh\nNr            29    # number of rows of nodes in mesh\nS_mtv         5.0   # top solubility (verticle, only for model_nbr 1)\nS_mth         5.0   # top solubility (horizontal, only for model_nbr 1)\nz_break       1.3   # z-break (only for model_nbr 1)')


# ## Run software

# In[ ]:


get_ipython().system(' ../ttial')


# ## Visualize results
# 
# ### Concentration

# In[ ]:


def plotConcentrations():
    height, width = getDimensions('output.txt')
    fig, ((ax1),(ax2)) = plt.subplots(nrows=1, ncols=2)

    conc_ver = np.loadtxt('conc_ver_matrix.txt')
    M,N = np.shape(conc_ver)
    im = ax1.matshow(conc_ver,cmap=cmaps)
    ax1.set_xticks([0,N/2,N-1])
    ax1.set_xticklabels(['0',str(width/2),str(width)])
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yticks([0,M/2,M-1])
    ax1.set_yticklabels(['0',str(height/2),str(height)])
    ax1.tick_params(right=False,left=False,top=False,bottom=False)
    ax1.set_title('conc ver matrix', fontsize=12)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)

    conc_hor = np.loadtxt('conc_hor_matrix.txt')
    M,N = np.shape(conc_hor)
    im2 = ax2.matshow(conc_hor,cmap=cmaps)
    ax2.set_xticks([0,N/2,N-1])
    ax2.set_xticklabels(['0',str(width/2),str(width)])
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_yticklabels({})
    ax2.tick_params(right=False,left=False,top=False,bottom=False)
    ax2.set_title('conc hor matrix', fontsize=12)
    ax2.set_xlabel('x', fontsize=12)

    maxV = np.max([np.max(conc_ver),np.max(conc_hor)])
    minV = np.min([np.min(conc_ver),np.min(conc_hor)])
    if minV < 0:
        print('Warning! Negative concentration')
    minV = 0
    cbticks = [minV,0.5*(minV+maxV),maxV]
    cblabels = [str(minV),str(0.5*(minV+maxV)),str(maxV)]
    cb_ax = fig.add_axes([0.95, 0.25, 0.04, 0.5])
    cbar = mpl.colorbar.ColorbarBase(cb_ax, ticks=cbticks, boundaries=np.linspace(minV,maxV,1000), cmap=cmaps,norm=mpl.colors.Normalize(vmin=minV, vmax=maxV))
    cbar.ax.set_yticklabels(cblabels, fontsize=12)
    plt.savefig('concentrations.pdf',dpi=300,bbox_inches='tight')
plotConcentrations()


# ### Flux

# In[ ]:


def plotFlux():
    height, width = getDimensions('output.txt')
    fig, ((ax1),(ax2)) = plt.subplots(nrows=1, ncols=2)

    j_ver = np.loadtxt('j_ver_matrix.txt')
    M,N = np.shape(j_ver)
    im = ax1.matshow(j_ver,cmap=cmaps)
    ax1.set_xticks([0,N/2,N-1])
    ax1.set_xticklabels(['0',str(width/2),str(width)])
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yticks([0,M/2,M-1])
    ax1.set_yticklabels(['0',str(height/2),str(height)])
    ax1.tick_params(right=False,left=False,top=False,bottom=False)
    ax1.set_title('j ver matrix', fontsize=12)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)

    j_hor = np.loadtxt('j_hor_matrix.txt')
    M,N = np.shape(j_hor)
    im2 = ax2.matshow(j_hor,cmap=cmaps)
    ax2.set_xticks([0,N/2,N-1])
    ax2.set_xticklabels(['0',str(width/2),str(width)])
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_yticklabels({})
    ax2.tick_params(right=False,left=False,top=False,bottom=False)
    ax2.set_title('j hor matrix', fontsize=12)
    ax2.set_xlabel('x', fontsize=12)

    maxV = np.max([np.max(j_ver),np.max(j_hor)])
    minV = np.min([np.min(j_ver),np.min(j_hor)])
    cbticks = [minV,0.5*(minV+maxV),maxV]
    cblabels = [str(minV),str(0.5*(minV+maxV)),str(maxV)]
    cb_ax = fig.add_axes([0.95, 0.25, 0.04, 0.5])
    cbar = mpl.colorbar.ColorbarBase(cb_ax, ticks=cbticks, boundaries=np.linspace(minV,maxV,1000), cmap=cmaps,norm=mpl.colors.Normalize(vmin=minV, vmax=maxV))
    cbar.ax.set_yticklabels(cblabels, fontsize=12)
    plt.savefig('flux_matrix.pdf',dpi=300,bbox_inches='tight')
    
    fig, (ax) = plt.subplots(nrows=1, ncols=1)
    j_ver = np.loadtxt('j_ver_vector.txt')
    j_hor = np.loadtxt('j_hor_vector.txt')
    ax.plot(j_ver[:,0],j_ver[:,1], lw=1, color='blue',ls='-',label=r'Vertical')
    ax.plot(j_hor[:,0],j_hor[:,1], lw=1, color='red',ls='-',label=r'Horizontal')
    ax.plot(np.linspace(0,height,2),[0,0], lw=1, color='black',ls='-')
    ax.set_xlim((0, height))
    ax.set_xlabel(r'Depth [m]', fontsize=12)
    ax.set_ylabel(r'Flux [kg/m$^3$]', fontsize=12)
    ax.legend(loc='center right', frameon=False,fontsize=12)
    plt.savefig('flux_vector.pdf',dpi=300,bbox_inches='tight')
plotFlux()


# ### Absolute activity

# In[ ]:


def plotAbsoluteActivity():
    height, width = getDimensions('output.txt')
    fig, (ax1) = plt.subplots(nrows=1, ncols=1)

    absactivity = np.loadtxt('V_matrix.txt')
    M,N = np.shape(absactivity)
    im = ax1.matshow(absactivity,cmap=cmaps)
    ax1.set_xticks([0,N/2,N-1])
    ax1.set_xticklabels(['0',str(width/2),str(width)])
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yticks([0,M/2,M-1])
    ax1.set_yticklabels(['0',str(height/2),str(height)])
    ax1.tick_params(right=False,left=False,top=False,bottom=False)
    ax1.set_title('V matrix', fontsize=12)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)

    maxV = np.max(absactivity)
    minV = np.min(absactivity)
    cbticks = [minV,0.5*(minV+maxV),maxV]
    cblabels = [str(minV),str(0.5*(minV+maxV)),str(maxV)]
    cb_ax = fig.add_axes([0.8, 0.1, 0.04, 0.75])
    cbar = mpl.colorbar.ColorbarBase(cb_ax, ticks=cbticks, boundaries=np.linspace(minV,maxV,1000), cmap=cmaps,norm=mpl.colors.Normalize(vmin=minV, vmax=maxV))
    cbar.ax.set_yticklabels(cblabels, fontsize=12)
    plt.savefig('solubility.pdf',dpi=300,bbox_inches='tight')
plotAbsoluteActivity()


# ### Diffussion coefficient

# In[ ]:


def plotDiffusionCoefficients():
    height, width = getDimensions('output.txt')
    fig, ((ax1),(ax2)) = plt.subplots(nrows=1, ncols=2)

    diffcoeff_ver = np.loadtxt('Dv_matrix.txt')
    M,N = np.shape(diffcoeff_ver)
    im = ax1.matshow(diffcoeff_ver,cmap=cmaps)
    ax1.set_xticks([0,N/2,N-1])
    ax1.set_xticklabels(['0',str(width/2),str(width)])
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yticks([0,M/2,M-1])
    ax1.set_yticklabels(['0',str(height/2),str(height)])
    ax1.tick_params(right=False,left=False,top=False,bottom=False)
    ax1.set_title('Dv matrix', fontsize=12)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)

    diffcoeff_hor = np.loadtxt('Dh_matrix.txt')
    M,N = np.shape(diffcoeff_hor)
    im2 = ax2.matshow(diffcoeff_hor,cmap=cmaps)
    ax2.set_xticks([0,N/2,N-1])
    ax2.set_xticklabels(['0',str(width/2),str(width)])
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_yticklabels({})
    ax2.tick_params(right=False,left=False,top=False,bottom=False)
    ax2.set_title('Dh matrix', fontsize=12)
    ax2.set_xlabel('x', fontsize=12)

    maxV = np.max([np.max(diffcoeff_ver),np.max(diffcoeff_hor)])
    minV = np.min([np.min(diffcoeff_ver),np.min(diffcoeff_hor)])
    if minV < 0:
        print('Warning! Negative diffusion coefficients')
    cbticks = [minV,0.5*(minV+maxV),maxV]
    cblabels = [str(minV),str(0.5*(minV+maxV)),str(maxV)]
    cb_ax = fig.add_axes([0.95, 0.25, 0.04, 0.5])
    cbar = mpl.colorbar.ColorbarBase(cb_ax, ticks=cbticks, boundaries=np.linspace(minV,maxV,1000), cmap=cmaps,norm=mpl.colors.Normalize(vmin=minV, vmax=maxV))
    cbar.ax.set_yticklabels(cblabels, fontsize=12)
    plt.savefig('diffusion_coefficients.pdf',dpi=300,bbox_inches='tight')
plotDiffusionCoefficients()


# ### Solubility

# In[ ]:


def plotSolubility():
    height, width = getDimensions('output.txt')
    fig, ((ax1),(ax2)) = plt.subplots(nrows=1, ncols=2)

    s_ver = np.loadtxt('sv_matrix.txt')
    M,N = np.shape(s_ver)
    im = ax1.matshow(s_ver,cmap=cmaps)
    ax1.set_xticks([0,N/2,N-1])
    ax1.set_xticklabels(['0',str(width/2),str(width)])
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yticks([0,M/2,M-1])
    ax1.set_yticklabels(['0',str(height/2),str(height)])
    ax1.tick_params(right=False,left=False,top=False,bottom=False)
    ax1.set_title('sv matrix', fontsize=12)
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)

    s_hor = np.loadtxt('sh_matrix.txt')
    M,N = np.shape(s_hor)
    im2 = ax2.matshow(s_hor,cmap=cmaps)
    ax2.set_xticks([0,N/2,N-1])
    ax2.set_xticklabels(['0',str(width/2),str(width)])
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_yticklabels({})
    ax2.tick_params(right=False,left=False,top=False,bottom=False)
    ax2.set_title('sh matrix', fontsize=12)
    ax2.set_xlabel('x', fontsize=12)

    maxV = np.max([np.max(s_ver),np.max(s_hor)])
    minV = np.min([np.min(s_ver),np.min(s_hor)])
    cbticks = [minV,0.5*(minV+maxV),maxV]
    cblabels = [str(minV),str(0.5*(minV+maxV)),str(maxV)]
    cb_ax = fig.add_axes([0.95, 0.25, 0.04, 0.5])
    cbar = mpl.colorbar.ColorbarBase(cb_ax, ticks=cbticks, boundaries=np.linspace(minV,maxV,1000), cmap=cmaps,norm=mpl.colors.Normalize(vmin=minV, vmax=maxV))
    cbar.ax.set_yticklabels(cblabels, fontsize=12)
    plt.savefig('solubility.pdf',dpi=300,bbox_inches='tight')
plotSolubility()


# In[ ]:




