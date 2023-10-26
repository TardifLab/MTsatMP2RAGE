# MNI Sample Protocols

**The scripts are prefilled with the imaging parameters provided in the sample protocols from the Montreal Neurological Institute.**

## Notes:
These scripts make the following assumptions:
- You already have loaded your images into matlab
- All images are in register and have the same matrix size (eg, 256x256x192 ); including the B1 map
- Preprocessing steps have been done (eg, Gibbs unringing and denoising)

If you are not sure how to load your images into matlab, or about other preprocessing, please see ProcessingScripts/Sample_code_calculate_all_maps.m 

# B1 correction:
- Sequence simulations for correcting MTsat have already been run for these protocols.
- Figures from the simulations and fitting are in MNI_sampleProtocols/b1Correction/simOutputs/Figures
- The fitValues files are located in MNI_sampleProtocols/b1Correction/simOutputs
- hMRI spoiling simulation results have been hard coded into the matlab scripts. But the results are in MNI_sampleProtocols/b1Correction/simOutputs/hMRI
- Sample input data in minc format can be found in MNI_sampleProtocols/sampleData ... coming soon
- Reference output data can be found in MNI_sampleProtocols/b1Correction/simOutputs/processing
