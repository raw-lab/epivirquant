#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from skimage.filters import threshold_otsu
from skimage.color import label2rgb
from skimage import measure, restoration, exposure
from pypher.pypher import zero_pad

def corr(corrCount, XB, shortname, PSF, nLR_iter, szMetric, px2nm, outDir):
    print("Starting processing for: " + shortname)
    startCorr = time.time();
    XB_avgPx = np.mean(XB)                # Get mean pixel value in current working image 
    currOut = os.path.join(outDir, "Step-3_Corr", f"CORR{corrCount}_{shortname}")         
    if os.path.exists(currOut) == False:  # Check if current out directory exists         
       os.mkdir(currOut);                 # Create current out directory
    dname = os.path.join(outDir, "Step-3_Corr", f"CORR{corrCount}_{shortname}")
    plt.imshow(XB,cmap='gray')
    plt.title("XB-"+shortname)
    plt.savefig(os.path.join(dname, f"XB-{shortname}.png"), dpi=150)
    plt.close('all')
    L = np.size(XB,0)             # Get length of current DAPI image
    dim = np.size(XB,1)           # Get dim of current DAPI image
    pad_L = int(L*(16/L))         # Pad size of y-axis by experimentally acquired ratio
    pad_dim = int(dim*(16/dim))   # Pad size of x-axis by experimentally acquired ratio
    XB = zero_pad(XB,(L+pad_L,dim+pad_dim),position='center')  # Pad for discontinuities
    XB = restoration.richardson_lucy(XB, PSF, nLR_iter)  # Perform LR algorithm
    XB = XB[pad_L:L,pad_dim:dim]                       # Truncate to avoid edge effects
    XB = zero_pad(XB,(L,dim),position='center')        # Pad back to original dimensions
    XB[XB == 0] = XB_avgPx                             # Set px = 0 to mean px
    plt.figure("RL-deconvblind")
    plt.imshow(XB,cmap='gray'), plt.title("XB-"+shortname+": deconvoluted")
    plt.savefig(os.path.join(dname, f"XB-{shortname}_Deconvoluted.png"), dpi=150)
    plt.close('all')
    threshold = threshold_otsu(XB)    # Get current working image Otsu threshold
    otsuMask = XB > threshold         # Binarize image using Otsu threshold
    labels = measure.label(otsuMask)  # Acquire labels from binary objects
    rprops = measure.regionprops(labels,intensity_image=XB)  # Acquire object properties
    dVec = []              # Initialize list to populate with object diameters
    iVec = []              # Initialize list to populate with object mean intensities
    eVec = []              # Initialize list to populate with object linear eccentricites 
    oIdx = []              # Initialize list to populate with object idx to omit
    majorVec = []          # Initialize list to populate with semi-major axis lengths
    oCount = 0             # Initialize counter
    for prop in rprops:    # Loop over number of object properties
        if szMetric == 1:  # If user-defined size metric option is equal to 1
            currMajor = prop.axis_major_length    # Current object semi-major axis length 
            currMajorSc = currMajor*px2nm         # Scale semi-major axis for omit check
            majorVec.append(currMajor)            # Store current semi-major axis length
            eqD = prop.equivalent_diameter_area   # Get current equivalent diameter area
            scEqD = eqD*px2nm                     # Scale diameter by nm-to-px ratio
            dVec.append(scEqD)                    # Store current mean axis length
        if szMetric == 2:                         # If size metric option set to 2
            currMajor = prop.axis_major_length    # Current object semi-major axis length
            currMajorSc = currMajor*px2nm         # Scale semi-major axis for omit check
            majorVec.append(currMajor)            # Store current semi-major axis length
            currMinor = prop.axis_minor_length    # Current object semi-minor axis length
            currDia = (currMajor + currMinor)/2   # Compute avg. of semi-major/minor axes
            currDia = currDia*px2nm               # Scale diameter by nm-to-px ratio
            dVec.append(currDia)                  # Store current scaled mean axis length
        iVec.append(prop.intensity_mean)          # Store current object mean intensity
        ecc = prop.eccentricity                   # Get objects' linear eccentricity
        eVec.append(ecc)                          # Store current linear eccentricity
        oCount = oCount + 1                       # Update counter
    currAvg = round(np.mean(dVec),2)              # Get current image avg. object size
    numObj = len(dVec)                            # Get current FITC image # of objects  
    plt.imshow(otsuMask,cmap='gray'), plt.title("XB-"+shortname+": Otsu mask")
    plt.savefig(os.path.join(dname, f"XB-{shortname}_OtsuMask.png"), dpi=150)
    labelsColored = label2rgb(labels,colors=["blue"],bg_label=0)  # Color objects
    labelsColored = exposure.adjust_gamma(labelsColored,gamma=0.1,gain=1)  # Scale int.
    plt.imshow(labelsColored), plt.title("XB-"+shortname+": Labels")
    plt.savefig(os.path.join(dname, f"XB-{shortname}_Labels.png"), dpi=150)
    plt.close('all')

    tempInt = iVec  # Get current list of object intensities
    intensity =  np.array(tempInt)*1000     # Scale intensities for convenience 
    tempSize = dVec     # Get working image list of object sizes
    size = np.array(tempSize)   # Convert list to array
    [m,b] = np.polyfit(tempSize,tempInt,1)    # Fit curve to intensity vs. size data
    x = np.linspace(np.min(size),np.max(size),np.size(size))  # Create vector of x-ticks
    y = m*x + b     # Compute line given slope & y-intercept
    plt.plot(x,y,'k--',markersize=14)
    plt.axvline(currAvg,color='b',linestyle =":")
    plt.text(currAvg+(currAvg*0.01),np.min(tempInt),"\u03BC = "+str(currAvg)+" nm",
             size=12,fontweight='bold',color='blue')
    plt.scatter(tempSize,tempInt,c=tempInt,cmap='hot',edgecolors='black')
    plt.ylabel("Mean intensity (arb. unit)"), plt.xlabel("Object diameter (nm)")
    plt.title("XB-"+str(corrCount)+" Object mean intensity vs. diameter")
    plt.savefig(os.path.join(dname, f"XB-{corrCount}_IvsD.png"), dpi=150)
    plt.close('all')    # Close all open figures
    plt.scatter(tempSize,eVec,c=eVec,cmap='hot',edgecolors='black')
    plt.ylabel("Linear eccentricity")                                                     
    plt.xlabel("Object diameter (nm)")
    plt.title(f"XB-{corrCount} Linear eccentricty vs. diameter")
    plt.savefig(os.path.join(dname, f"XB-{corrCount}_EvsD.png"), dpi=150)
    plt.close('all')    # Close all open figures

    print("Processing ended for: " + shortname)
    return shortname, corrCount, currAvg, iVec, dVec, m, b, numObj, round(time.time()-startCorr, 4)

def get_corrolation(outDir, PSF, nLR_iter, szMetric, px2nm, XB0, XBsnames, nCorr, sphereSize, cpus):
    avgVec = []    # Initialize list to populate with object size averages
    intVec = []    # Initialize list to populate with object intensities
    szVec = []     # Initialize list to populate with object sizes
    mVec = []      # Initialize list to populate with intensity vs diameter slopes
    bVec = []      # Initialize list to populate with intensity vs diameter y-intercepts
    objVec = []    # Initialize list to populate with number of objects identified per image

    print(" "); print('Step 3 begin:');     # Begin step 3
    print("Performing LR algorithm for", nLR_iter, "iterations...")
    start = time.time();    # Get current UTC time in seconds

    if os.path.exists(os.path.join(outDir, "Step-3_Corr")) == False:    # Check if Step 3 out dir exists  
        os.mkdir(os.path.join(outDir, "Step-3_Corr"))   # Create Step 3 output directory

    from joblib import Parallel, delayed
    results = Parallel(n_jobs=cpus)(delayed(corr)(ii, XB0[ii], XBsnames[ii], PSF, nLR_iter, szMetric, px2nm, outDir) for ii in range(len(XB0)))

    for shortname, corrCount, currAvg, iVec, dVec, m, b, numObj, runtime in results:
        print("Evaluated image", corrCount+1, "/", nCorr, ":", shortname, "in", runtime, "seconds");
        avgVec.append(currAvg)  # Store current image avg. object size
        intVec.append(iVec)
        szVec.append(dVec)  # Store current list of object sizes
        mVec.append(m)
        bVec.append(b)
        objVec.append(numObj)

    # Generate PDF
    x = np.sort(np.array(sum(szVec,[])))
    mu = round(np.mean(x),2); sigma = round(np.std(x),2)
    std1 = mu - sigma; std2 = mu + sigma
    y = stats.norm.pdf(x,mu,sigma)
    plt.plot(x,y,'-o',color='k',markersize=3)
    [xmin,xmax,ymin,ymax] = plt.axis()
    n = np.max(y)/ymax
    plt.axvline(mu,ymax=n,color='b',linestyle =":",lw=2)
    plt.axvline(std1,color='r',linestyle =":",lw=2)
    plt.axvline(std2,color='r',linestyle =":",lw=2)
    plt.grid()
    plt.xlim(0, 500)
    ymin = min(y); ymax = max(y)
    plt.fill_betweenx((ymin,ymax),std1,std2,color='gray',alpha=0.4)
    plt.xlabel("Object diameter (nm)"), plt.ylabel("P(d)")
    plt.title("XB-all: object size probability density")
    plt.legend(["Diameters (d)","\u03BC = "+str(mu)+" nm",
                "\u03C3 = \u00B1 "+str(sigma)+" nm"],loc='best',fontsize=12)
    plt.savefig(os.path.join(outDir, "Step-3_Corr", "XB_SizeDistribution.png"), dpi=150)
    plt.close('all')    # Close all open figures

    #Plot Histogram
    n, bins, patches = plt.hist(x,bins=21,edgecolor='#e0e0e0',linewidth=0.4,alpha=0.8,zorder=3)
    n = n.astype('int') # it MUST be integer
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.viridis(n[i]/max(n))) 
    plt.grid(axis='y',zorder=0)
    plt.xlim(0, 500)
    plt.axvline(mu, color='0',linestyle ="--", label=r'$\mu_{obj}$' + f" = {mu:.2f}", zorder=5)
    plt.xlabel("Object diameter (nm)"), plt.ylabel("Frequency")
    plt.title("XB-all: object size histogram")
    plt.savefig(os.path.join(outDir, "Step-3_Corr", "XB_SizeHistogram.png"), dpi=150)
    plt.close('all')    # Close all open figures
    
    #Final Correction Factor
    XB_mu = np.mean(avgVec)    # Get mean size of all objects over all DAPI images
    CORR = sphereSize/XB_mu    # Compute fluorescence size correction coefficient
    print("XB correction factor: "+str(round(CORR,4)))
    logCorr(XBsnames,objVec,avgVec,CORR,XB_mu,nCorr, sigma, outDir)     # Log correction info
    stop = time.time() - start;     # Get correction run time
    print("Step 3 end: "+str(round(stop,4))+" seconds"); print(" ")
    print('o-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-o');
    return CORR


def logCorr(XBsnames, objVec, avgVec, CORR, XB_mu, nCorr, sigma, outDir):
    with open(os.path.join(outDir, "Step-3_Corr", "corrLog.txt"), 'w+') as corrLog:
        corrLog.write("|>==============================================<| \n")
        corrLog.write("|>================== CORR Log ==================<| \n")
        for iii in range(len(XBsnames)):
            corrLog.write(" "+XBsnames[iii]+" \n")
            corrLog.write(" Number of objects: "+str(objVec[iii])+" \n")
            corrLog.write(" Mean object size:  "+str(avgVec[iii])+" nm \n")
            corrLog.write("|................................................| \n")
            corrLog.write("|>==============================================<| \n")
        corrLog.write(" Mean size over "+str(len(XBsnames))+" images:      "+
                      str(round(XB_mu,2))+" nm \n")
        corrLog.write(" Measurment error: " + str(sigma) + " \n")
        corrLog.write(" Objects identified per image: "
                      +str(round(np.sum(objVec)/nCorr,2))+" obj/img \n")
        corrLog.write(" Total number of objects: " + str(np.sum(objVec)) + " \n")
        corrLog.write(" Final correction factor:      "+str(round(CORR,6))+" \n")
        corrLog.write("|>==============================================<| \n")
        corrLog.write("|................................................| \n")