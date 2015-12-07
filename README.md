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

3. On the plot that renders the embedded dataset, write a module that takes two inputs -- a list of coordinates $\left\{\left(x_i, y_i\right)\right\}$, and the coordinate (x,y) of an arbitrary new point -- then outputs a list of indces with corresponding weights.

4. The last module takes as input all meshes in the collection and the list of indices with weights, and "synthesize" a new mesh.

Once all four modules are available, we can build correspondences between two far-apart meshes in four steps: 1) join the pair with a straight line; 2) sample on this line equi-distantly, and synthesize an artificial mesh at each point; 3) compute correspondences between each consecutie pair of meshes on the straight line; 4) "knit-together" the pairwise correspondences on consecitive small hops to build a global one.
