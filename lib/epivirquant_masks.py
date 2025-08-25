import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from skimage.filters import threshold_otsu
from skimage.color import label2rgb
from skimage import measure
from skimage import restoration
from skimage import exposure
from pypher.pypher import zero_pad

fsep = '/'

def generate_masks(outDir, XG0, XGsnames, nF, PSF, nLR_iter, szMetric, px2nm, CORR, SM_constraint):
    print(" "); print('Step 4 begin:');
    start = time.time();
    if os.path.exists(outDir+fsep+"Step-4_genMasks") == False:
        os.mkdir(outDir+fsep+"Step-4_genMasks")
    count = 1
    avgVec = []     # object size averages
    intVec = []     # object intensities
    szVec = []      # object sizes
    mVec = []       # intensity vs diameter slopes
    bVec = []       # intensity vs diameter y-intercepts
    objVec = []     # objects identified
    xyVec = []
    omit = []       # omission indices
    train0 = []     # binary masks for NN training/testing
    test0 = []      # deconvoluted FITC images
    for i in range(len(XG0)):
        XG = XG0[i]
        XG_avgPx = np.mean(XG)  # mean pixel value in current working image
        shortname = XGsnames[i];
        currOut = outDir+fsep+"Step-4_genMasks"+fsep+"genMask"+str(count)+"_"+shortname
        if os.path.exists(currOut) == False:
            os.mkdir(currOut);
        dname = outDir+fsep+"Step-4_genMasks"+fsep+"genMask"+str(count)+"_"+shortname;
        print("Evaluating image "+str(count)+"/"+str(nF)+": "+shortname)
        plt.imshow(XG,cmap='gray')
        plt.title("XG-"+shortname)
        plt.savefig(dname+fsep+"XG-"+shortname+".png",dpi=150)
        plt.close('all')
        print("Performing LR algorithm for "+str(nLR_iter)+" iterations...")
        L = np.size(XG,0)       # length of current FITC image
        dim = np.size(XG,1)     # dim of current FITC image
        pad_L = int(L*(32/L))        # Pad size of y-axis 
        pad_dim = int(dim*(24/dim))  # Pad size of x-axis
        XG = zero_pad(XG,(L+pad_L,dim+pad_dim),position='center')  # Pad for discontinuities
        XG = restoration.richardson_lucy(XG,PSF,nLR_iter)
        XG = XG[pad_L:L,pad_dim:dim]    # Truncate to avoid edge effects
        XG = zero_pad(XG,(L,dim),position='center') # Pad back to original dimensions
        XG[XG == 0] = XG_avgPx  # Set px = 0 to mean px
        train0.append(XG)
        plt.figure("RL-deconvblind")
        plt.imshow(XG,cmap='gray'), plt.title("XG-"+shortname+": deconvoluted")
        plt.savefig(dname+fsep+"XG-"+shortname+"_Deconvoluted.png",dpi=150)
        plt.close('all')

        print("Scanning FITC image: "+shortname)
        threshold = threshold_otsu(XG)    # Get Otsu threshold
        otsuMask = XG > threshold         # Apply Otsu threshold
        labels = measure.label(otsuMask)  # Acquire labels
        plt.imshow(otsuMask,cmap='gray'), plt.title("XG-"+shortname+": Otsu mask")
        plt.savefig(dname+fsep+"XG-"+shortname+"_OtsuMask.png",dpi=150)
        rprops = measure.regionprops(labels,intensity_image=XG)  # Acquire object properties 
        dVec = []   # Object diameters
        iVec = []   # Object mean intensities
        eVec = []   # Object linear eccentricites
        oIdx = []   # Object idx to omit
        majorVec = []   # Semi-major axis lengths
        coords = [] # Object point coordinates
        oCount = 0
        for prop in rprops:
            if szMetric == 1:
                currMajor = prop.axis_major_length
                currMajorSc = currMajor*px2nm*CORR  # Scale semi-major axis
                majorVec.append(currMajor)
                eqD = prop.equivalent_diameter_area # Current equivalent diameter area
                scEqD = eqD*px2nm   # Scale diameter by nm-to-px ratio
                scEqD = scEqD*CORR  # Scale diameter correction coefficient
                dVec.append(scEqD)     
            if szMetric == 2:
                currMajor = prop.axis_major_length
                currMajorSc = currMajor*px2nm*CORR  # Scale semi-major axis
                majorVec.append(currMajor)
                currMinor = prop.axis_minor_length  # Semi-minor axis length 
                currDia = (currMajor + currMinor)/2 # Compute avg. of semi-major/minor axes
                currDia = currDia*px2nm # Scale diameter by nm-to-px ratio
                currDia = currDia*CORR  # Scale diameter correction coefficient
                dVec.append(currDia)
            coords.append(prop.coords)
            iVec.append(prop.intensity_mean)
            ecc = prop.eccentricity # Linear eccentricity
            eVec.append(ecc)
            if currMajorSc == 0 or ecc > 0.9999:
                oIdx.append(oCount) # Store index to omit
            if currMajorSc > SM_constraint:
                oIdx.append(oCount)
            oCount = oCount + 1    
            oIdx = [*set(oIdx)]; oIdx.sort(reverse=True)  # Remove duplicates and sort elements
        # Remove false positives from all lists
        for j in oIdx:
            del dVec[j]
            del iVec[j]
            del eVec[j]
            del coords[j]
        currAvg = round(np.mean(dVec),2)
        avgVec.append(currAvg)  # Avg. object size  
        szVec.append(dVec)      # Object sizes
        intVec.append(iVec)     # Object intensities
        numObj = len(dVec)      # Number of objects in image
        objVec.append(numObj)
        omit.append(oIdx)       # Current omission index
        xyVec.append(coords)
        print("Number of objects: "+str(numObj))                                              
        print(shortname+" mean size: "+str(currAvg)+" nm")
        labelsColored = label2rgb(labels,colors=["green"],bg_label=0)  # Color objects
        labelsColored = exposure.adjust_gamma(labelsColored,gamma=0.1,gain=1)
        plt.imshow(labelsColored), plt.title("XG-"+shortname+": Labels")
        plt.savefig(dname+fsep+"XG-"+shortname+"_Labels.png",dpi=150)
        plt.close('all')    
        # i vs d Plots 
        tempInt = intVec[i]
        intensity =  np.array(tempInt)*1000 # Scale intensities 
        tempSize = szVec[i]
        size = np.array(tempSize)
        [m,b] = np.polyfit(tempSize,tempInt,1)  # Fit curve to intensity vs. size data
        mVec.append(m)
        bVec.append(b)
        x = np.linspace(np.min(size),np.max(size),np.size(size))  # Create vector of x-ticks
        y = m*x + b # Compute line given slope & y-intercept

        plt.plot(x,y,'k--',markersize=14)
        plt.axvline(currAvg,color='b',linestyle =":")
        plt.text(currAvg+(currAvg*0.01),np.min(tempInt),"\u03BC = "+str(currAvg)+" nm",
                size=12,fontweight='bold',color='blue')
        plt.scatter(tempSize,tempInt,c=tempInt,cmap='hot',edgecolors='black')
        plt.ylabel("Mean intensity (arb. unit)"), plt.xlabel("Object diameter (nm)")
        plt.title("XG-"+str(count)+" Object mean intensity vs. diameter")
        plt.savefig(dname+fsep+"XG-"+str(count)+"_IvsD.png",dpi=150)
        plt.close('all')

        plt.scatter(tempSize,eVec,c=eVec,cmap='hot',edgecolors='black')
        plt.ylabel("Linear eccentricity"), plt.xlabel("Object diameter (nm)")
        plt.title("XG-"+str(count)+" Linear eccentricty vs. diameter")
        plt.savefig(dname+fsep+"XG-"+str(count)+"_EvsD.png",dpi=150)
        plt.close('all')
        count = count + 1
    
    # Generate PDF 
    x = np.sort(np.array(sum(szVec,[])))
    mu = round(np.mean(x),2); sigma = round(np.std(x),2)
    std1 = mu - sigma; std2 = mu + sigma
    y = stats.norm.pdf(x,mu,sigma)
    plt.plot(x,y,'-o',color='k',markersize=3)
    [xmin,xmax,ymin,ymax] = plt.axis(); n = np.max(y)/ymax
    plt.axvline(mu,ymax=n,color='b',linestyle =":",lw=2)
    plt.axvline(std1,color='r',linestyle =":",lw=2)
    plt.axvline(std2,color='r',linestyle =":",lw=2)
    plt.xlabel("Object diameter (nm)"), plt.ylabel("P(d)")
    plt.title("XG-all: object size probability density")
    plt.legend(["Diameters (d)","\u03BC = "+str(mu)+" nm",
                "\u03C3 = \u00B1 "+str(sigma)+" nm"],loc='best',fontsize=12)
    plt.savefig(outDir+fsep+"Step-4_genMasks"+fsep+"XG_SizeDistribution.png",dpi=150)
    plt.close('all')
    # Partition Counts 
    XG_mu = np.mean(avgVec)  # Get mean size of all objects over all DAPI images
    print("Mean size over "+str(len(XGsnames))+" images: "+str(round(XG_mu,2))+" nm")
    nSz1 = 0; nSz2 = 0; nSz3 = 0; nSz4 = 0; nSz5 = 0
    for jj in range(len(szVec)):
        currSizes = np.array(szVec[jj])
        for kk in range(len(currSizes)):
            if currSizes[kk] < 100:
                nSz1 = nSz1 + 1
            if 100 <= currSizes[kk] <= 220:
                nSz2 = nSz2 + 1
            if 220 < currSizes[kk] < 500:
                nSz3 = nSz3 + 1
            if 500 <= currSizes[kk] <= 1200:
                nSz4 = nSz4 + 1
            if 1200 < currSizes[kk] < 3000:
                nSz5 = nSz5 + 1
    print("Total # of objects identified: "+str(np.sum(objVec)))
    print("Size breakdown:")
    print("    < 100 nm:           "+str(nSz1)+" objects")
    print("    100 nm to 220 nm:   "+str(nSz2)+" objects")
    print("    220 nm to 500 nm:   "+str(nSz3)+" objects")
    print("    500 nm to 1200 nm:  "+str(nSz4)+" objects")
    print("    1200 nm to 3000 nm: "+str(nSz5)+" objects")
    
    logCount0(XGsnames,objVec,avgVec,XG_mu,nF,outDir,fsep,nSz1,nSz2,nSz3,nSz4,nSz5)
    exportTSV(XGsnames, szVec, xyVec, outDir, fsep)
    stop = time.time() - start;
    print("Step 4 end: "+str(round(stop,4))+" seconds"); print(" ")
    print('o-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-o')


# Log Initial Count
def logCount0(XGsnames,objVec,avgVec,XG_mu,nF,outDir,fsep,nSz1,nSz2,nSz3,nSz4,nSz5):
    with open(outDir+fsep+"Step-4_genMasks"+fsep+"countLog.txt",'w+') as countLog:
        countLog.write("|>===============================================<| \n")
        countLog.write("|>================== Count Log ==================<| \n")
        for jjj in range(len(XGsnames)):
            countLog.write(" "+XGsnames[jjj]+" \n")
            countLog.write(" Number of objects: "+str(objVec[jjj])+" \n")
            countLog.write(" Mean object size:  "+str(avgVec[jjj])+" nm \n")
            countLog.write("|.................................................| \n")
            countLog.write("|>===============================================<| \n")
        countLog.write(" Mean size over "+str(len(XGsnames))+" images: "+
                       str(round(XG_mu,2))+" nm \n")
        countLog.write(" Total # of objects identified: "+str(np.sum(objVec))+" \n")
        countLog.write(" Size breakdown: \n")
        countLog.write("     < 100 nm:           "+str(nSz1)+" objects \n")         
        countLog.write("     100 nm to 220 nm:   "+str(nSz2)+" objects \n")
        countLog.write("     220 nm to 500 nm:   "+str(nSz3)+" objects \n")
        countLog.write("     500 nm to 1200 nm:  "+str(nSz4)+" objects \n")               
        countLog.write("     1200 nm to 3000 nm: "+str(nSz5)+" objects \n")
        countLog.write("|>===============================================<| \n")
        countLog.write("|.................................................| \n")
# (filename, size, x coord, y coord)
def exportTSV(XGsnames,szVec, xyVec, outDir,fsep):
    with open(outDir+fsep+"Step-4_genMasks"+fsep+"sizeCoords.tsv",'w+') as countLog:
        for i, name in enumerate(XGsnames):
            for j in range(len(szVec[i])):
                countLog.write(name+"\t"+str(szVec[i][j])+"\t"+str(xyVec[i][j][0][0])+"\t"+str(xyVec[i][j][0][1])+"\n")
