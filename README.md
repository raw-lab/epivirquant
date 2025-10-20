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

Citing EpiVirQuant
-------------
If you are publishing results obtained using EpiVirQuant

BioRxiv
Figueroa III JL, Hollenack SM, Bellanger M, Fulghum B, Visscher PT, White III RA. 2024  <br />
Resolving and Quantifying Viral-Like Particles via Blind Deconvolution  <br />
[BioRxiv](https://doi.org/10.1101/2024.04.21.590467)  <br />



Contact
-------

The informatics point-of-contact for this project is [Dr. Richard Allen White III](https://github.com/raw-lab).  
If you have any questions or feedback, please feel free to get in touch by email.  
[Dr. Richard Allen White III](mailto:rwhit101@charlotte.edu)<br /> 
[Jose Luis Figueroa III](mailto:jlfiguer@charlotte.edu) <br />
[Sadie Marie Hollenack](mailto:bhollena@charlotte.edu) <br />
[Bryan Fulghum](mailto:bfulghu2@charlotte.edu)  <br />
[Madeline Bellanger](mailto:mbellang@charlotte.edu) <br />

Or [open an issue](https://github.com/raw-lab/epivirquant/issues). 
