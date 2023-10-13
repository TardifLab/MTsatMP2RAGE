# MNI Sample Protocols

**The scripts are prefilled with the imaging parameters provided in the sample protocols from the Montreal Neurological Institute.**

## Notes:
These scripts make the following assumptions:
- You already have loaded your images into matlab
- All images are in register and have the same matrix size (eg, 256x256x192 ); including the B1 map
- Preprocessing steps have been done (eg, Gibbs unringing and denoising)

If you are not sure how to load your images into matlab, or about other preprocessing, please see ProcessingScripts/Sample_code_calculate_all_maps.m 
