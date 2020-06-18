# Linear Fit to Produce T1-maps from Variable Flip Angle method
Note: make sure you have the Image Processing Toolbox installed (https://www.mathworks.com/products/image.html).

This repo contains two matlab scripts:
# vfat1maplinear.m 
This will load your SPGR images (as dicoms) in the folder specified (currently hard-coded, sorry) and the flip angles specified (also hard-coded to 2,5,10,15 degrees) and create a VFA T1 map using a linear fit. It will save the T1 map as a dicom.


# vfat1maplinear_b1correction.m 
This version includes one extra step - B1 correction to improve the T1 estimation. B1 inhomogeneities impact the actual flip angle that mangetization in a voxel rotates due to on-resonance RF pulses. It's often represented as a normalized correction factor of the nominal flip angle (nominal means the flip angle you actually intended & set at the scanner). i.e.:

\alpha_actual = B1_correction factor.* \alpha_nominal

An important note about the B1 maps used in this script, with many thanks to Emil Ljungberg and Mathieu Boudreau  KCL for some key information and resources:
The B1 maps this script is written to import were produced by 2db1map. 2db1map is built on the Bloch Siegert method by Sacolick et al. For more info about the Bloch Siegert method, check out the resources in this issue https://github.com/qMRLab/qMRLab/issues/394 from Mathieu Boudreau.

The 2db1map would produce B1 maps in absolute units of mT which is calculated from two images obtained by the acquisition (with and without bloch siegert prep pulse). The implementation used to produce the B1 data used here gives you a flip angle scaling map (an image with voxel value = 1 where the nominal flip angle is achieved and lower than one where the nominal flip angle is lower). This means you have an image with decimal values, which is what you want as a scaling factor (B1_correction) but which is a problem since DICOM only supports integer values from 0-2^15. So GE has scaled the image by 10x the flip angle of the 2db1map acquisition (*in degrees*). Key info: it is multiplied by the flip angle you see from dicominfo of the B1 map (not your VFA images). This was 20 degrees in this case.
