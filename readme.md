# DirichletSteklovLandmarkMatching
This repository contains the code for the paper ["Non-Isometric Shape Matching via Functional Maps on Landmark-Adapted Bases"]() by Mikhail Panine, Maxime Kirgo and Maks Ovsjanikov.

In this paper, we propose a principled approach for non-isometric landmark-preserving non-rigid shape matching. Our method is based on the functional maps framework, but rather than promoting isometries we focus instead on near-conformal maps that preserve landmarks exactly. We achieve this, first, by introducing a novel landmark-adapted basis using an intrinsic Dirichlet-Steklov eigenproblem. Second, we establish the functional decomposition of conformal maps expressed in this basis. Finally, we formulate a conformally-invariant energy that promotes high-quality landmark-preserving maps, and show how it can be solved via a variant of the recently proposed ZoomOut method that we extend to our setting. Our method is descriptor-free, efficient and robust to significant mesh variability. We evaluate our approach on a range of benchmark datasets and demonstrate state-of-the-art performance on non-isometric benchmarks and near state of the art performance on isometric ones.

<p align="center">
  <img align="center"  src="/figures/teaser.png", width=800>
</p>

The main steps involved in our method to map a source shape to a target shape is described in Sect. 3 of the paper and illustrated below. The source shape is highlighted in orange and the target shape in blue.

<p align="center">
  <img align="center"  src="/figures/pipeline.png", width=800>
</p>


Main Functions
--------------
```matlab
Steklov_settings = compute_steklov_settings(num_landmarks, NN_type,InitialGuess,DS_num_eigs,radii_factor,weight_Orthonormality,weight_Proper,weight_Bijectivity,num_LB_eigs,ZO_start,ZO_step);

% Input:
%	num_landmarks: number of landmarks to use
%	NN_type: whether to use our fast approximation or the principled definition ('principled','fast'), 'fast' is our preferred option
%	InitialGuess: type of initial guess ('normal_derivatives','trivial','landmark_harmonics','conformal_energy'), 'normal_derivatives' is our preferred method
%	DS_num_eigs: number of Dirichlet-Steklov eigen-functions
%	radii_factor: size of radius around each landmark
%	weight_Orthonormality: orthonormality weight
%	weight_Proper: properness weight
%	weight_Bijectivity: bijectivity weight
%	num_LB_eigs: number of Laplace-Beltrami eigenfunctions
%	ZO_start: ZoomOut refinement start
%	ZO_step: number of ZoomOut refinement steps
%
% Output:
%   Steklov_settings: the Settings object used to parameterize our method
```

```matlab
[Src_refined,Tar_refined,fullp2pTarSrc, fullp2pSrcTar,fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = compute_steklov(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings);
% Input:
%	Src: the source shape structure (see shape_matching.m)
%	Src_landmarks: the landmarks on the source shape
%	Tar: the target shape structure (see shape_matching.m)
%	Tar_landmarks: the landmarks on the target shape
%	Steklov_settings: the Settings object used to parameterize our method (see function above)
% Output:
%	Src_refined: the source shape structured with circular boundaries at the landmarks
%	Tar_refined: the target shape structured with circular boundaries at the landmarks
%	fullp2pTarSrc: the p2p map from the original (i.e. before introducing circular boundaries) target shape to the original source shape, before ZoomOut refinement
%	fullp2pSrcTar: the p2p map from the original (i.e. before introducing circular boundaries) source shape to the original target shape, before ZoomOut refinement
%	fullp2pTarSrc_ZO: the p2p map from the original (i.e. before introducing circular boundaries) target shape to the original source shape, after ZoomOut refinement
%	fullp2pSrcTar_ZO: the p2p map from the original (i.e. before introducing circular boundaries) source shape to the original target shape, after ZoomOut refinement
```

Comments
--------
- The script ```shape_matching.m``` shows how to perform shape matching using our method.
- The scripts ```reproduce_fig_12_step[1-3].mat``` allow to reproduce Fig. 12 left of our paper:
	- First, clone the Faust remeshed dataset shapes (```vtx_5k``` folder) to the folder ```./data/FAUST/vtx_5k/```
	- ```reproduce_fig_12_step1.mat``` computes the remeshed shapes, wrapped into .mat files;
	- ```reproduce_fig_12_step2.mat``` computes the geodesic errors on all shapes;
	- ```reproduce_fig_12_step3.mat``` plots the figure using the computed data.


Acknowledgments
----------------
- This work was funded by EDF R&D and the ANRT as part of the CIFRE grant 2019/0433.
- The ```+MESH``` module was developped by [Jing Ren](https://github.com/llorz).
- The ```gputoolbox``` module was developped by [Alec Jacobson](https://github.com/alecjacobson/gptoolbox).

Citation
--------
If this code contributes to academic work, please cite:
```bib
@inproceedings{panine2022dirichletsteklov,
  title={Non-Isometric Shape Matching via Functional Maps on Landmark-Adapted Bases},
  author={Panine, Mikhail and Kirgo, Maxime and Ovsjanikov, Maks},
  booktitle={Computer Graphics Forum},
  year={2022},
  organization={Wiley Online Library}
}
```

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/). For any commercial uses or derivatives, please contact us (mikhail.panine@usi.ch, maximekirgo@gmail.com, maks@lix.polytechnique.fr).
