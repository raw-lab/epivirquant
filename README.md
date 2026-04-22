# Welcome to EpiVirQuant

## About 

EpiVirQuant directly counts and measures the size of viral-like particles, bacteria, archaea in epifluorescence microscopy images commonly used in viral ecology, viral estimation and enumeration, and microbial ecology. Quantification and sizing is based on data-driven blind deconvolution with a tunable point-spread function. 

Running EpiVirQuant
-------------
[Install conda if you have not](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

To run using the data provided in ```./data```, run the ```runDemo.sh``` script in your terminal. This script will create a conda environment with all of the dependencies. 

If you want to create the environment without quantifying the demo data, follow these steps </br>

Create the mamba environment 
```bash
mamba create -y -n EpiVirQuant -c conda-forge joblib configargparse matplotlib numpy scipy scikit-image pypher stardist
```
Build the packages:
```bash
pip install .
```

Run Epivirquant
```bash
epivirquant.py --dapi "{PATH TO DAPI IMAGES}" --fitc "{PATH TO FITC IMAGES}" --calibration "{PATH TO CALIBRATION IMAGE}"
```


## 📄 License

Creative Commons Attribution-NonCommercial (CC BY-NC 4.0) — See LICENSE file

**Citations:**

If you are publishing results obtained using EpiVirQuant, please cite: <br />  
- Figueroa III JL, Hollenack SM, Bellanger M, Fulghum B, Visscher PT, White III RA. 2026.  <br />
Resolving and quantifying viral-like particles via blind deconvolution. [BMC Methods](https://doi.org/10.1186/s44330-026-00060-z). 3:10.  <br />

Pre-Print EpiVirQuant <br />
- Figueroa III JL, Bellanger M, Fulghum B, Visscher PT, White III RA. 2024  <br />
Resolving and Quantifying Viral-Like Particles via Blind Deconvolution  <br />
[BioRxiv](https://doi.org/10.1101/2024.04.21.590467)  <br />

---

## Contributing to EpiVirQuant

We welcome contributions of other experts expanding features in EpiVirQuant. Please contact us via support. 

---

## 📞 Support

- **Issues:** [open an issue](https://github.com/raw-lab/epivirquant/issues).  
- **Email:** [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br /> 
             [Jose Luis Figueroa III](mailto:jlfiguer@charlotte.edu) <br />
---

**Made with ❤️ to prevent people counting VLPs by eye for hours - for the viral ecology community**
