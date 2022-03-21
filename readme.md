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
y = f(x);

% Input:
%	x: input
%
% Output:
%   y: output
```

Comments
--------
- The script ```shape_matching.m``` shows how to perform shape matching using our method.


Acknowledgments
----------------
- This work was funded by EDF R&D and the ANRT as part of the CIFRE grant 2019/0433.
- The ```+MESH``` module was developped by [Jing Ren](https://github.com/llorz).

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
