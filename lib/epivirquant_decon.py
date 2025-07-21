#!/usr/bin/env/ python
# -*- coding: utf-8 -*-

__version__ = '0.1'
__date__ = '04-07-2024'
__authors__ = 'Richard Allen White III & Jose Luis Figueroa III'

import os
import numpy as np
import time
import matplotlib.pyplot as plt
from pypher.pypher import psf2otf
from scipy.ndimage import gaussian_filter
import math
from skimage import img_as_ubyte
from skimage.measure import shannon_entropy
from skimage.feature import graycomatrix, graycoprops
eps = np.finfo(float).eps 

def create_PSF(fSize, a, b, sigma, r, tau, v, s, psfMethod):
    x = np.linspace(-fSize, fSize, fSize)
    y = np.linspace(-fSize, fSize, fSize)
    x, y = np.meshgrid(x, y)
    PSFi = np.zeros((fSize, fSize))
    R = np.sqrt(np.square(x) + np.square(y))
    if psfMethod == "gam":
        for i in range(fSize):
            for j in range(fSize):
                #equation 6
                PSFi[i, j] = a * np.exp(v / ((math.gamma(1 + tau * R[i, j])) * (math.gamma(1 - tau * R[i, j])))) + s
        sumPSF = PSFi.sum()
        if sumPSF != 0:
            PSFi /= sumPSF
    return PSFi, x, y

def get_MLE(optBox, fSize, PSFi, nMLE_iter):
    L, dim = optBox.shape
    J = [optBox,optBox,0,np.zeros((L*dim,2))]
    P = [PSFi,PSFi,0,np.zeros((fSize*fSize,2))]   
    fw = np.fft.fft2(np.ones((L, dim)))
    nMLE_iter = 10  # Set desired number of MLE iterations (default 10)
    for k in range(nMLE_iter):
        J_k, P_k = update_MLE(J, P, k)
        OTF = psf2otf(P_k, (L, dim))
        plusNoise = np.real(np.fft.ifft2(OTF * np.fft.fft2(J_k)))
        plusNoise[plusNoise <= 0] = np.finfo(plusNoise.dtype).eps
        currEst = J[0] / (plusNoise + np.finfo(plusNoise.dtype).eps)
        Psi_k = np.fft.fft2(currEst)
        J[2] = J[1]
        OTF = psf2otf(P[1], (L, dim))

        scale = np.real(np.fft.ifft2(np.conj(OTF) * fw)) + np.sqrt(eps)
        OTFtemp = np.fft.ifft2(np.conj(OTF) * Psi_k)
        J[1] = np.maximum((J_k * np.real(OTFtemp) / scale), 0)
        J_k_flat = np.ndarray.flatten(J_k, order='F')
        J_flat = np.ndarray.flatten(J[1], order='F')
        J4 = J[3]
        Jtemp = J_flat - J_k_flat
        J[3] = np.vstack((Jtemp, J4[:, 0])).T
        P[2] = P[1]
        Jfreq = np.fft.fft2(J[2])
        OTF_k = np.conj(Jfreq) * fw
        scale = otf2psf(OTF_k, fSize) + np.sqrt(eps)
        Jfreq_k = otf2psf(np.conj(Jfreq) * Psi_k, fSize)+ np.sqrt(eps)
        P[1] = np.maximum((P_k * Jfreq_k) / scale, 0)
        P[1] = P[1] / np.sum(P[1])
        P_flat = np.ndarray.flatten(P[1], order='F')
        P_k_flat = np.ndarray.flatten(P_k, order='F')
        P4 = P[3]
        Ptemp = P_flat - P_k_flat
        P[3] = np.vstack((Ptemp, P4[:, 0])).T
        if k == (nMLE_iter - 1):
            PSFr = P[1]
            Xdec = J[1]
    return PSFr, Xdec

def update_MLE(J, P, k):
    if k != 0:
        alpha_k_J = (J[3][0:, 0].T @ J[3][0:, 1]) / (J[3][0:, 1].T @ J[3][0:, 1] + eps)
        alpha_k_J = np.maximum(np.minimum(alpha_k_J, 0), 0)
        J_k = np.maximum(J[1] + alpha_k_J * (J[1] - J[2]), 0)
        alpha_k_P = (P[3][:, 0].T @ P[3][:, 1]) / (P[3][:, 1].T @ P[3][:, 1] + eps)
        alpha_k_P = np.maximum(np.minimum(alpha_k_P, 0), 0)
        P_k = np.maximum(P[1] + alpha_k_P * (P[1] - P[2]), 0)
        P_k = P_k / np.sum(P_k)
    else:
        J_k = np.maximum(J[1], 0)
        P_k = np.maximum(P[1], 0)
        P_k = P_k / (np.sum(P_k) + eps)
    return J_k, P_k

def otf2psf(OTF, fSize):
    if np.sum(OTF) == 0:
        PSF = np.zeros_like(OTF)
    else:
        PSF = np.fft.ifft2(OTF)
        if np.max(abs(PSF.imag)) / np.max(abs(PSF)) <= np.finfo(PSF.dtype).eps:
            PSF = PSF.real
        PSF = np.roll(PSF, int(np.floor(fSize / 2)), axis=0)
        PSF = np.roll(PSF, int(np.floor(fSize / 2)), axis=1)
        PSF = PSF[0:fSize, 0:fSize]
    return PSF

def opt_metrics(Xdec, entVec, energyVec):
    Xdec_uint = img_as_ubyte(Xdec)
    S = shannon_entropy(Xdec_uint, base=2)
    entVec.append(S)
    GLCM = graycomatrix(Xdec_uint, [2], [0], symmetric=True, normed=True)
    energy = graycoprops(GLCM, "energy")[0, 0]
    energyVec.append(energy)
    return entVec, energyVec

def getTextures(Xdec,currOutDir):
    L = np.size(Xdec,0)          # Get length of current Xdec
    dim = np.size(Xdec,1)        # Get dim of current Xdec
    [x,y] = np.mgrid[0:L,0:dim]  # Create grid for surface plot
    f = plt.figure(); ax = f.add_subplot(projection='3d')
    ax.plot_surface(x,y,Xdec,rstride=1,cstride=1,cmap="turbo",linewidth=2)
    ax.view_init(9.25,0)         # Rotate surface
    plt.xlabel("Y (px)",fontsize=7), plt.ylabel("X (px)",fontsize=7)
    plt.tick_params(axis='both',labelsize=8)
    ax.set_zlabel("Intensity (arb. unit)",fontsize=7)
    plt.title("Final Virgilian Pair intensity textures"), plt.show()
    plt.savefig(currOutDir+os.sep+"VP_final_textures.png",dpi=300)
    plt.close('all')             # Close all open figures

def genPSFOTF(PSF,currOutDir):
    L = np.size(PSF,0)           # Get length of current OTF
    dim = np.size(PSF,1)         # Get dim of current OTF
    [x,y] = np.mgrid[0:L,0:dim]  # Create grid for surface plot
    f1 = plt.figure(); ax = f1.add_subplot(projection='3d')
    ax.plot_surface(x,y,PSF,rstride=1,cstride=1,cmap="coolwarm",linewidth=2)
    plt.ylabel("Y"); plt.xlabel("X"); ax.set_zlabel("h(x,y)")
    plt.title("Final point-spread function"), plt.show()
    plt.savefig(currOutDir+os.sep+"PSF_final.png",dpi=150)
    plt.close('all')             # Close all open figures
    OTF = psf2otf(PSF,(L,dim))   # Convert PSF to OTF
    f2 = plt.figure(); ax = f2.add_subplot(projection='3d')
    ax.plot_surface(x,y,np.real(OTF),rstride=1,cstride=1,cmap="RdYlBu",linewidth=2)
    plt.ylabel("Y"); plt.xlabel("X"); ax.set_zlabel("Z")
    plt.title("Final optical transfer function"), plt.show()
    plt.savefig(currOutDir+os.sep+"OTF_final.png",dpi=150)                                  
    plt.close('all')             # Close all open figures

def logDecon(outDir,Count,optStop,minS,maxE,minSE,fSize,tau,v):
    with open(outDir+os.sep+"Step-2_Decon"+os.sep+"deconLog.txt",'w+') as deconLog:
        deconLog.write("|>==============================================<| \n")
        deconLog.write("|>================= CyDecon Log ================<| \n")
        deconLog.write("Number of sweeps: "+str(Count)+" \n")
        deconLog.write("Sweep rate        "+str(Count/optStop)+" sweeps/second \n")
        deconLog.write("Optimal PSF selected: \n")
        deconLog.write("  Minimum entropy: "+str(minS)+" bits \n")
        deconLog.write("  Max GLCM energy: "+str(maxE)+" \n")
        deconLog.write("  min(S - E):      "+str(minSE)+" \n")
        deconLog.write("Final PSF parameters: \n")
        deconLog.write("  Filter size: "+str(fSize)+"-by-"+str(fSize)+" \n")
        deconLog.write("  tau:         "+str(tau)+" \n")
        deconLog.write("  v:           "+str(v)+" \n")
        deconLog.write("|>==============================================<| \n")
        deconLog.write("|................................................| \n")  

def decon(outDir, optBox, a, b, sig, r, tau, s, v, nMLE_iter, psfMethod):
    out_dir = outDir + os.sep + "Step-2_Decon"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print("\nStep 2 begin:")
    print("Loading gamma sinc PSF...")

    start_time = time.time()
    L, dim = optBox.shape
    f_max = min(L, dim)

    if f_max % 2 == 0:
        f_max -= 1

    filters = np.arange(3, f_max+1, 2)
    ent_vec = []
    energy_vec = []
    count = 0
    psf_params = []
    b_vec = []
    psf_vec = []
    #formula 7
    tau_vec = np.arange(0, 1/(10*math.pi), 1/(100*math.pi))
    v_vec = np.arange(0, 1*math.pi, 1/10*math.pi)
    opt_count = 1
    curr_out_dir = os.path.join(out_dir, "optimization")

    if not os.path.exists(curr_out_dir):
        os.mkdir(curr_out_dir)

    opt_start = time.time()
    TEMPCOUNT = 0
    tau = tau_vec[0]
    v = v_vec[0]
    for f_size in filters:
        tau = round(tau, 10)
        v = round(v, 10)
        print(f"fSize: {f_size}; tau: {tau}; v: {v}")
        PSFi, X, Y = create_PSF(f_size, a, b, sig, r, tau, v, s, psfMethod)
        PSFr, Xdec = get_MLE(optBox, f_size, PSFi, nMLE_iter)
        ent_vec, energy_vec = opt_metrics(Xdec, ent_vec, energy_vec)
        count += 1
        psf_params.append([f_size, tau, v, s])
        TEMPCOUNT+=1
        for tau in tau_vec:
            tau = round(tau, 10)
            v = round(v, 10)
            print(f"fSize: {f_size}; tau: {tau}; v: {v}")
            PSFi, X, Y = create_PSF(f_size, a, b, sig, r, tau, v, s, psfMethod)
            PSFr, Xdec = get_MLE(optBox, f_size, PSFi, nMLE_iter)
            ent_vec, energy_vec = opt_metrics(Xdec, ent_vec, energy_vec)
            count += 1
            psf_params.append([f_size, tau, v, s])

            for v in v_vec:
                tau = round(tau, 10)
                v = round(v, 10)
                print(f"fSize: {f_size}; tau: {tau}; v: {v}")
                PSFi, X, Y = create_PSF(f_size, a, b, sig, r, tau, v, s, psfMethod)
                PSFr, Xdec = get_MLE(optBox, f_size, PSFi, nMLE_iter)
                ent_vec, energy_vec = opt_metrics(Xdec, ent_vec, energy_vec)
                count += 1
                psf_params.append([f_size, tau, v, s])

    opt_stop = time.time() - opt_start
    

    print("|.............................................|")
    print(f"Number of sweeps: {count}")
    print(f"Sweep rate: {round(count / opt_stop, 4)} sweeps/second")
    print("Parameter sweeps complete.")

    minSidx = np.argmin(ent_vec)
    minS = ent_vec[minSidx]
    maxEidx = np.argmax(energy_vec)
    maxE = energy_vec[maxEidx]

    print("Optimal PSF selected:")
    print(f"  Minimum entropy: {minS} bits")
    print(f"  GLCM energy: {maxE}")

    opt_params = psf_params[minSidx]
    f_size = opt_params[0]
    tau = opt_params[1]
    v = opt_params[2]

    print("Final PSF parameters:")
    print(f"  Filter size: {f_size}-by-{f_size}")
    print(f"  tau: {tau}")
    print(f"  v: {v}")

    PSF, _, _ = create_PSF(f_size, a, b, sig, r, tau, v, s, psfMethod)
    PSF, _ = get_MLE(optBox, f_size, PSF, nMLE_iter)

    print("Final PSF constructed.")

    plt.scatter(np.arange(0, count), ent_vec, c=ent_vec, cmap='hot', edgecolors='black')
    plt.ylabel("Shannon entropy (bits)")
    plt.xlabel("Iteration")
    plt.autoscale()
    plt.xlim([0, count + 1])
    plt.title("Shannon entropy vs. iteration")
    plt.savefig(os.path.join(curr_out_dir, "EntVsIter.png"), dpi=150)
    plt.close('all')

    plt.scatter(np.arange(0, count), energy_vec, c=energy_vec, cmap='hot', edgecolors='black')
    plt.ylabel("Energy")
    plt.xlabel("Iteration")
    plt.autoscale()
    plt.xlim([0, count + 1])
    plt.title("GLCM energy vs. iteration")
    plt.savefig(os.path.join(curr_out_dir, "EnergyVsIter.png"), dpi=150)
    plt.close('all')

    getTextures(Xdec, curr_out_dir)
    genPSFOTF(PSF, curr_out_dir)
    #getDist(optBox, curr_out_dir, "B-measure_initial.png", os.sep)
    #getDist(Xdec, curr_out_dir, "B-measure_final.png", os.sep)

    fig, axes = plt.subplots(1, 2)
    axes[0].imshow(optBox, cmap="gray")
    axes[0].set_xlabel("X (px)")
    axes[0].set_ylabel("Y (px)")
    axes[0].set_title("Initial optimization box")
    axes[1].imshow(Xdec, cmap="gray")
    axes[1].set_xlabel("X (px)")
    axes[1].set_ylabel("Y (px)")
    axes[1].set_title("Final optimization box")
    plt.subplots_adjust(left=0.08, right=0.99, bottom=0.01, top=0.99)
    plt.savefig(outDir + os.sep + "Step-2_Decon" + os.sep + "InitialVsFinal_optBox.png", dpi=300)
    plt.close('all')
        #def logDecon(outDir,Count,optStop,minS,maxE,minSE,fSize,tau,v):
    logDecon(outDir, count, opt_stop, minS, maxE, f_size, tau, v)
    stop = time.time() - start_time
    print(f"Step 2 end: {round(stop, 4)} seconds")
    print('o-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-o')

