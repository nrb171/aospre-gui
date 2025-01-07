# Installing AOSPRE-GUI.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running AOSPRE-GUI](instructions.md)

# Prerequisites

## MATLAB or MATLAB Runtime (Required)
AOSPRE-GUI is built using MATLAB. While installing MATLAB with a license is easier, you can run MATLAB programs without a license using the MATLAB Runtime. So, you must have access to one of the following on your system: 
- [MATLAB](https://www.mathworks.com/products/matlab.html)
- [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html)

Assistance installing MATLAB is outside the scope of this document. If you encounter issues, contact your system administrator.

## AOSPRE (Recommended)
AOSPRE-GUI is a graphical user interface for setting up simulations using the AOSPRE package. However, for especially large simulations, it may be preferable to use the AOSPRE-GUI software on a local machine and then transfer the generated namelists, scanning tables, and scripts to a different (presumably more powerful) machine with AOSPRE installed. 

You can download the AOSPRE package from the [AOSPRE GitHub repository](https://github.com/NCAR/AOSPRE). Follow the installation instructions located here: [AOSPRE Installation](https://github.com/NCAR/AOSPRE/blob/main/docs/building.md).

## Python (Optional)
There is a standalone Python script that can be used to rename WRF output files ([wrfoutToAOSPRE.py](../helpers/wrfoutToAOSPRE.py)). This is not necessary as there is a MATLAB function in the GUI that can perform the same action, but it is included as a template in case your simulation file does not follow the default wrfout naming conventions. 

Any version of Python 3 will be sufficient.

# Installation
1. Open a terminal and navigate to the directory where you want to install AOSPRE-GUI.
2. Clone the AOSPRE-GUI repository using the following command:
```bash
git clone https://github.com/nrb171/aospre-gui/
```
3. Navigate to the `aospre-gui` directory:
```bash
cd aospre-gui
```
And you're ready to move onto the next section ([instructions](instructions.md))!