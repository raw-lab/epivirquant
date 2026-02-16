#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.spatial.distance import cdist
from skimage import measure
from skimage.filters import threshold_otsu

def getVP(outDir, snameVP, XBVP, dConstraint, pad, px2nm):
    fsep = os.sep
    if not os.path.exists(outDir + fsep + "Step-1_VP"):
        os.mkdir(outDir + fsep + "Step-1_VP")
    
    plt.figure("XB-" + snameVP)
    plt.imshow(XBVP, cmap='gray')
    plt.title("XB-" + snameVP)
    plt.savefig(f"{outDir}{fsep}Step-1_VP{fsep}XB-{snameVP}", dpi=150)

    print('o-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-o\n')
    print("Step 1 begin:\nScanning calibration image for VP candidates...")
    
    start = time.time()
    threshold = threshold_otsu(XBVP)
    otsuMask = XBVP > threshold
    plt.figure("XB-" + snameVP + ": Otsu mask")
    plt.imshow(otsuMask, cmap='gray')
    plt.title("XB-" + snameVP + ": Otsu mask")
    plt.savefig(f"{outDir}{fsep}Step-1_VP{fsep}XB-{snameVP} - Otsu mask", dpi=150)
    
    otsuMask_eroded = ndimage.binary_erosion(otsuMask, iterations=1)
    
    plt.figure("XB-" + snameVP + ": Otsu mask eroded")
    plt.imshow(otsuMask_eroded, cmap='gray')
    plt.title("XB-" + snameVP + ": Otsu mask eroded")
    plt.savefig(f"{outDir}{fsep}Step-1_VP{fsep}XB-{snameVP} - Otsu mask eroded", dpi=150)
    
    labels = measure.label(otsuMask_eroded)
    rprops = measure.regionprops(labels, intensity_image=XBVP)

    # calculates the centers and bounding box of each vlp
    bBox = np.array([prop.bbox for prop in rprops])
    cVec = np.array([prop.centroid for prop in rprops])
    # Creates a matrix of the distances between every point
    dVec = np.array(cdist(cVec, cVec))
    # List to ensure we don't go over the same pair of points twice
    vpsVisited: tuple[int, int] = []
    optBox = []
    minDist = 999999999999
    numVPs = 0
    # Find potential VP candidates
    for i, currPoint in enumerate(dVec):
        for j, dist in enumerate(currPoint):
            if 0 < dist < dConstraint and (i,j) not in vpsVisited and (j,i) not in vpsVisited:
                vpsVisited.append((i,j))
                # get the bounding box of both, then crop the image to fit both points
                point1 = bBox[i]; point2 = bBox[j]
                yMin = min(point1[0], point2[0]) - pad
                yMax = max(point1[2], point2[2]) + pad
                xMin = min(point1[1], point2[1]) - pad
                xMax = max(point1[3], point2[3]) + pad
                # if any of the points are off the image due to padding, skip it
                if yMin < 0 or yMax < 0 or xMin < 0 or xMax < 0:
                    continue
                VP = XBVP[yMin:yMax, xMin:xMax]
                # Double check the potential VP actually has 2 objects
                threshold = threshold_otsu(VP)
                otsuMask = VP > threshold
                otsuMask_eroded = ndimage.binary_erosion(otsuMask, iterations=1)
                labels = measure.label(otsuMask_eroded)
                rprops = measure.regionprops(labels, intensity_image=VP)
                # the opt box will be the pair that is the closest together
                if len(rprops) == 2 and dist < minDist:
                    optBox = VP
                    minDist = dist
                    numVPs += 1
   
    print(f"Number of VP candidates identified: {numVPs}")
    if numVPs == 0:
        sys.exit(f"Epivirquant ERROR: zero VPs identified within {dConstraint * px2nm} nm. Increase distance constraint or select new calibration image.")
    else:
        print(f"VP saved to: {outDir}{fsep}Step-1_VP{fsep}")

    plt.figure("Initial Optimization Box")
    plt.imshow(optBox, cmap='gray')
    #plt.plot([x2, y2], [x1, y1],color="blue",marker="o",markerfacecolor="r")
    plt.title("VP_Final")
    plt.savefig(f"{outDir}{fsep}Step-1_VP{fsep}XB-{snameVP}_VP_Final.png", dpi=150)
    plt.close('all')
    
    print(f"Closest proximity at VP: {round(minDist * px2nm, 2)} nm")
    
    stop = time.time() - start
    print(f"\nStep 1 end: {round(stop, 4)} seconds")
    print('\no-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-o')
    return optBox