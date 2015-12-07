# ShapeID

--------------
## Interpolation between two meshes ##

1. Change paths and names of the meshes in ReparametrizationDemo.m, and run the script. You should see a folder "meshes/" that contains subfolder "reparametrized".

2. Change paths and names of the meshes in PoissonMeshInterpDemo.m accordingly, and run the script. The resulting interpolated meshes can be found under subfolder "reparametrized".

--------------
## Building Correspondences by Composition along Small Hops ##

This experiment is divided into four modules. The original input is a collection of triangular meshes; the goal is to build correspondences maps between any pair of meshes.

1. Run any pairwise registration algorithm to obtain an affinity matrix for the entire collection. This module takes a collection of meshes as input, and for each pair outputs a (dis-)similarity score.

2. Embed the affinity matrix into $\mathbb{R}^2$ using MDS or diffusion maps. This module takes an affinity matrix and outputs embedding coordinates for each mesh.
