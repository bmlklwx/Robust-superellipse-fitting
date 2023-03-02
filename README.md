# A Robust Superellipse Fitting Algorithm

A 2D version of the EMS algorithm: [CVPR 2022] Robust and Accurate Superquadric Recovery: a Probabilistic Approach.
> [**Robust and Accurate Superquadric Recovery: a Probabilistic Approach**](https://arxiv.org/abs/2111.14517 "ArXiv version of the paper.")  
> Weixiao Liu, Yuwei Wu, [Sipu Ruan](https://ruansp.github.io/), [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)

Our [original work](https://github.com/bmlklwx/EMS-superquadric_fitting) is to fit [superquadrics](https://en.wikipedia.org/wiki/Superellipsoid) (3D generalization of superellipse) to point clouds.
This is a simple variant to the original paper to solve [superellipse](https://en.wikipedia.org/wiki/Superellipse) (also know as Lamé curve) fitting problem in 2D cases.
The demo (test_script.m) shows the fitting results to randomly generated superellipse-shaped point clouds, with large amount of noise and outliers.
This repo also contains MATLAB functions to sample points almost uniformly on the side of superellipse, and to draw superellipse.

<img src="/figures/demo1.jpg" alt="superquadrics1" width="250"/><img src="/figures/demo2.jpg" alt="superquadrics2" width="250"/><img src="/figures/demo3.jpg" alt="superquadrics3" width="250"/>

<img src="/figures/demo4.jpg" alt="superquadrics4" width="250"/><img src="/figures/demo5.jpg" alt="superquadrics5" width="250"/><img src="/figures/demo6.jpg" alt="superquadrics6" width="250"/>

For visitors interested in more complex 3D superquadrics fitting, please visit this [repository](https://github.com/bmlklwx/EMS-superquadric_fitting).

If you find this repo useful, please cite

> W. Liu, Y. Wu, S. Ruan and G. S. Chirikjian, "Robust and Accurate Superquadric Recovery: a Probabilistic Approach," <br />
> 2022 IEEE/CVF Conference on Computer Vision and  Pattern Recognition (CVPR), New Orleans, LA, USA, 2022, pp. 2666-2675, <br />
> doi: 10.1109/CVPR52688.2022.00270.
