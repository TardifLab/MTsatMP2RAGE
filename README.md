# MTsatMP2RAGE

**MATLAB code for calculating T1 maps as in the manuscript, and running simmulations.**

## Requirements:

- Requires the qMRlab toolbox (http://qmrlab.org/) which can be downloaded here: https://github.com/qMRLab/qMRLab
- Uses the following simulation code: https://github.com/TardifLab/OptimizeIHMTimaging
- The simulations for B1 correction are based on the following toolbox and manuscript, but the above link contains the most updated functions, so use those! https://github.com/TardifLab/MTsatB1correction
- Assumes an MP2RAGE is used for T1 mapping, and thus for the noise calculations. So we need the MP2RAGE code from: https://github.com/JosePMarques/MP2RAGE-related-scripts
- Contains some misc code such as colormaps, and the necessary MP2RAGE code for correcting M0 https://github.com/christopherrowley/NeuroImagingMatlab
  - Specificalling in /QuantitativeFitting/MP2RAGE-related-scripts-master/Rowley_custom
- The polyFitN toolbox: https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
- And a symbolic toolbox: https://www.mathworks.com/matlabcentral/fileexchange/9577-symbolic-polynomial-manipulation 
---

## Data

Sample data for calculating MTsat with MP2RAGE is provided in /sampleData

- The data is stored as individual matlab matrices, and they have been skull stripped and cropped to reduce data size
- They have been labeled as "input*" or "output*", with the output files to be used as a reference.

Folder b1/ contains the necessary simulation outputs to b1 correct the MTsat data with the model-based approach

## Code Overview:

Code is broken up into sections depending on what you are looking to do.

- **MTeffect_in_T1mapping/** the main scripts for running the simulations, as well as the resulting values were saved.

- **ProcessingScripts/** sample scripts for calculating maps

- **simNewSingleProtocol/** simpified code for simulating the B1+ correction factor for different acquisition parameters

- **b1/** contains the calculated simulation results for B1 correction factors used in this study. These should work for your data if you use the same acquisition parameters in the MT-weighted image (including sat pulse time + flip angle)

The inversion recovery data for the phantom was fit using the qMRlab GUI. I was not able to adapt the code to get it to run outside of that.

**If used, please reference the following publication:**

Rowley CD, Nelson MC, Campbell JSW, Leppert IR, Pike GB, Tardif CL. Fast magnetization transfer saturation imaging of the brain using MP2RAGE T1 mapping. https://doi.org/10.1002/mrm.30143

If you used the simulation code, please reference the following:

Karakuzu, Agah, et al. "qMRLab: Quantitative MRI analysis, under one umbrella." Journal of Open Source Software 5.53 (2020).

Rowley CD, Campbell JSW, Leppert IR, Nelson MC, Pike GB, Tardif CL. Optimization of acquisition parameters for cortical inhomogeneous magnetization transfer (ihMT) imaging using a rapid gradient echo readout. Magnetic Resonance in Medicine. 2023; https://onlinelibrary.wiley.com/doi/10.1002/mrm.29754

If you use the MP2RAGE calculation code, please also cite: 

Marques, José P., et al. "MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field." Neuroimage 49.2 (2010): 1271-1281.

**For additional help with this code, contact rowleycd@mcmaster.ca**
