#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import math
import numpy as np
from skimage import io, img_as_float
import argparse

from epivirquant_lib import (epivirquant_pairing, epivirquant_decon, 
                             epivirquant_calibration, epivirquant_masks)
__version__ = '0.1'
__date__ = '04-07-2024'
__authors__ = 'Richard Allen White III, Sadie M. Hollenack, Jose Luis Figueroa III'

def main():
    args = parser.parse_args()

    if os.path.exists(args.outDir):
        print("Deleting old output folder: " + args.outDir)
        shutil.rmtree(args.outDir)
    os.mkdir(args.outDir)

    # Step 0: Load Images
    XB0, XG0, XBVP, snameVP, nCorr, XBsnames, XGsnames, nF = load_images(args.dapi, args.fitc, args.calibration)

    # Step 1: Get Optimizaion Box
    optBox = epivirquant_pairing.getVP(args.outDir, snameVP, XBVP, args.pad)

    # Step 2: Blind Deconvolution 
    PSF = epivirquant_decon.decon(args.outDir, optBox, args.a, args.b, args.sig, args.r, args.tau, args.s, args.v, args.nMLE_iter, "gam")
    
    # Step 3: Get the Correction Coefficent and Quantify the DAPI images
    px2nm = args.scaleMetric / args.scaleLength 
    CORR = epivirquant_calibration.get_corrolation(args.outDir, PSF, args.nLR_iter, args.szMetric, px2nm, XB0, XBsnames, nCorr, args.sphereSize)
    
    # Step 4: Quantify the FITC images
    epivirquant_masks.generate_masks(args.outDir, XG0, XGsnames, nF, PSF, args.nLR_iter, args.szMetric, px2nm, CORR, args.SM_constraint)

    print("EpiVirQuant complete.")

parser = argparse.ArgumentParser(add_help=False)
parser.set_defaults()
required = parser.add_argument_group('Input file(s) required')
required.add_argument('--dapi', action='append', default=[], help='Path to DAPI images directory')
required.add_argument('--fitc', action='append', default=[], help='Path to FITC images directory')
required.add_argument('--calibration', action='append', default=[], help='Path to calibration DAPI image')

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--genFigs', type=bool, default=True, help='Toggle the generation of figures on (True) or off (False). [True]')
optional.add_argument('--sphereSize', type=float, default=175, help="Diameter of microspheres in nanometers.")
optional.add_argument('--scaleLength', type=float, default=585, help='Length of scale bar for imaging equipment in pixels (px). [585]')
optional.add_argument('--scaleMetric', type=float, default=20e3, help="Represented length of scale bar in nm.")
optional.add_argument('--pad', type=int, default=14, help="Set padding around VP centroids to expand bounding box. [14od]")
optional.add_argument('--dConstraint', type=int, default=30, help="Set user-defined px distance constraint for VP candidates. [30]")
optional.add_argument('--fSize', type=int, default=31, help="Set dimensions of hybrid point-spread function, fSize-by-fSize. [31]")
optional.add_argument('--psfMethod', type=str, default='gam', choices=['gam', 'hyb', 'gau'], help="Set formula for point-spread function creation. The available options are: 'gam' (gamma sinc fn); 'hyb' (hybrid sinc fn); and 'gau' (gaussian). [gam]")
optional.add_argument('--nMLE_iter', type=int, default=10, help="Set maximum likelihood estimation (MLE) number of iterations. [10]")
optional.add_argument('--a', type=int, default=1, help="Parameter to control amplitude for gaussian component of GL PSF. [1]")
optional.add_argument('--b', type=int, default=1, help="Parameter to control amplitude for sinc component of hybrid PSF. [1]")
optional.add_argument('--sig', type=float, default=1, help="Parameter to control std of gaussian component of hybrid PSF. [1]")
optional.add_argument('--r', type=float, default=2.718281828, help="Parameter to control width for sinc component of hybrid PSF. [2.718281828]")
optional.add_argument('--tau', type=float, default=0.5, help="Parameter to control periodicity for sinc component of GL PSF. [0.5]")
optional.add_argument('--v', type=float, default=2.25, help="Parameter to control GL-PSFi (initial PSF) vertical stretch. [2.25]")
optional.add_argument('--s', type=float, default=0.0, help="Parameter to control GL-PSFi (initial PSF) vertical shift. [0.0]")
optional.add_argument('--tauVec', type=list, default=np.arange(1/(100*math.pi), 11/(100*math.pi), 1/(100*math.pi)), help="Tau parameters [1/(100*math.pi), 11/(100*math.pi), 1/(100*math.pi)]")
optional.add_argument('--vVec', type=list, default=np.arange(1/math.pi, 1.1*math.pi, 1/math.pi), help="Vector with v parameters [1/math.pi, 1.1*math.pi, 1/math.pi]")
#>----------------------------------| Size Correction |----------------------------------<#
optional.add_argument('--nLR_iter', type=int, default=80, help="Set number of iterations for Lucy-Richardson algorithm. [80]")
optional.add_argument('--szMetric', type=int, default=1, help="""Set metric to determine identified object sizes. The available
options are:
   1) Equivalent diameter area: the diameter of a circle having
      the same area as the identified object area (default).
   2) The average of the identified objects' semi-major and
      semi-minor axes. [1]""")
optional.add_argument('--SM_constraint', type=int, default=8000, help="Set threshold of semi-major axis to determine false positives. [8000]")
optional.add_argument('--outDir', type=str, default="EpiVirQuant_Output", help="Set EpiVirQuant output directory name. [EpiVirQuant_Output]")
optional.add_argument('--version', '-v', action='version',
                    version=f'Cyclops: \n version: {__version__} {__date__}',
                    help='show the version number and exit')
optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

def load_images(path_DAPI, path_FITC, path_CALIB):
    fsep = os.sep
    print("Initializing parameters and dirs...")
    fpath1 = path_DAPI[0]
    fpath2 = path_FITC[0]

    fnamesXB = os.listdir(fpath1)
    fnamesXG = os.listdir(fpath2)

    fExtension = os.path.splitext(fnamesXB[0])[-1]

    nCorr = len(fnamesXB)
    nF = len(fnamesXG)
    print("Loading images from DAPI and FITC dirs...")
    XBfnames = [f"{fpath1}{fsep}{name}" for name in fnamesXB]
    XBsnames = [name[:-len(fExtension)].replace(" ", "").replace("blue beads ", "BB-").replace(" + ", "+") for name in fnamesXB]
    XB0 = [img_as_float(io.imread(fname)) for fname in XBfnames]
    XGfnames = [f"{fpath2}{fsep}{name}" for name in fnamesXG]
    XGsnames = [name[:-len(fExtension)].replace(" ", "").replace("blue beads ", "BB-").replace(" + ", "+") for name in fnamesXG]
    XG0 = [img_as_float(io.imread(fname)) for fname in XGfnames]
    fpath3 = path_CALIB[0]
  
    currVP = io.imread(fpath3)
    XBVP = img_as_float(currVP)
    snameVP = fpath3.split(fsep)[-1].replace(" ", "").replace("bluebeads ", "BB-")[:-len(fExtension)]
    print(f"Number of DAPI images: {nCorr}")
    print(f"Number of FITC images: {nF}")
    print("Initialization complete.\n")
    return XB0, XG0, XBVP, snameVP, nCorr, XBsnames, XGsnames, nF

if __name__ == "__main__":
    main()