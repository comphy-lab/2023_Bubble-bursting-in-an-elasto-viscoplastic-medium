# Bursting bubble in an elasto-viscoplastic medium

This repository contains the codes used for simulating bubble bursting in an elasto-viscoplastic medium. 
The article is currently under preparation.

## Prerequisite

Basilisk C can be installed from [here](http://basilisk.fr). The installation steps are outlined [here](http://basilisk.fr/src/INSTALL). \
In case of compatibility issues, please feel free to contact me: [argb@mech.kth.se](mailto:argb@mech.kth.se). \
For post-processing codes, Python 3.X is required. 

## Compilation

The code can be compiled as: 
```bash
qcc -O2 -Wall -fopenmp burst_evp_new.c -o burst_evp_new -lm
```

## Running the code

The arguments to the code are: $\mathcal{J}$, $De$. For $\mathcal{J}=0.1,De=1.0$,
```bash 
./burst_evp_new 0.1 1.0 > out.log 2>&1 &
```
