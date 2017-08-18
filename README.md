# Quaternion-Sparse-Coding
Sparse coding and dictionary learning based on quaternion algebra

## Short descriptions
_Quaternion dictionary training in KQSVD_:
* The main functions include main_DictTraining.m and main_Reconstruction.m.
- main_DictTraining.m: Calling function for dictionary training step in K-QSVD algorithm. The trained disctionary is generated in folder of "training result".
- main_Reconstruction.m: Calling function for color image reconstruction using K-QSVD algorithm.
* The other called functions: Type "help XXX" will show the help information of this function. 
* Subfolder of Dataset classifying: For Storage of the input training image and reconstruction image. Here the training images are “flower” while the reconstructed images are “animal”, which verifies the robustness of our method to the variety of patches’ source. 
* Subfolder of Training result: for Storage of the trained dictionary.


_Color Image denoising using KQSVD_:
* main_Denoising.m: Calling function for Color Image De-noising using K-QSVD.
* The other called functions: Type "help XXX" will show the help information of this function.
* Subfolder of original: For Storage of the input image.
* Subfolder of Training result: For Storage of the trained dictionary without noise. 

## Reference and citations
Please cite following papers for utilizing the codes:
* Yu, Licheng, et al. "[Quaternion-based sparse representation of color image.](http://ieeexplore.ieee.org/abstract/document/6607436/)" Multimedia and Expo (ICME), 2013 IEEE International Conference on. IEEE, 2013.
* Xu, Yi, et al. "[Vector sparse representation of color image using quaternion matrix analysis.](http://ieeexplore.ieee.org/abstract/document/7024169/)" IEEE Transactions on Image Processing 24.4 (2015): 1315-1329.
