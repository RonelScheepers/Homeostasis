# Homeostasis
This repository contains the MATLAB and SINGULAR code for the mathematical analysis and associated figures in the publication "Robust Homeostasis of Cellular CHolesterol via Antithetic Control", submitted for review to the journal Fronties in Cell and Developmental Biology.

Authors: Ronel Scheepers and Robyn P. Araujo, May 2023, Queensland University of Technology, Brisbane, Australia.

THe mass-action analysis was performed by implementing the universal algorithm developed by Araujo and Liotta (Universal Structures for Adaptation in Biochemical Reaction Networks, 2023) via the open-source software Singular. 
The following instructions describe how this software can be access or downloaded:
The open–source software Singular can either be accessed via Jupyterlab, an extensible environment for interactive and reproducible computing that is based on the Jupyter Notebook and architecture, or via running Singular from the computer’s Terminal. In the first instance, follow the steps below to run Singular as a web–based interface on a local computer:
Download Anaconda from https://www.anaconda.com/products/distribution,
Open Anaconda Navigator from applications or through search,
Open JupyterLab from Navigator,
Install Jupyter-kernel-singular on terminal by copying code on  https://anaconda.org/conda-forge/jupyter-kernel-singular,

In the launcher, select Singular under the Notebook tab. The active file is now ready for the code to be entered.
If the code does not run, install Python by following the instructions on the Singular homepage http://www.singular.uni-kl.de/, under the Graphical Interface tab,
Finally, check that Python is in PATH by opening Terminal and typing the command: echo $PATH.

To run Singular from the Terminal, the application can be downloaded from
https://www.singular.uni-kl.de/index.php/singular-download.html

The numerical analysis in Section 4.3 was performed on data generated using the files labelled as StabilityOriginalParameters.mlx and StabilityRandomParameters.mlx respectively while the MATLAB script Chol_Input_Tiled.m generated the Figures 4a - 4d.


The code listings for analysing RPA in the cholesterol concentration of the plasma membrane (Cp), the free cholesterol (Cf) and esterified cholesterol (E) is provided in the SingularRPA_Cp.pdf, Singular_Cf.pdf and Singular_E.pdf files respectively.

For any further questions, please contact ronel.scheepers@hdr.qut.edu.au
