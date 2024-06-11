# Boundary Element Method (BEM) for 3D problem

Implementation of BEM for Earth's gravity potential problem by solving the Laplace equation $\nabla^{2}\ u = 0$, where the boundary is represented by the Earth's surface. This program solves the discretized boundary integral equations to find the potential on the boundary. First it sets up the linear system of equations which is then solved using the BiCGStab linear solver. Some sections are parallelized using the OpenMP library to make the computations faster. The validity of the achieved results was subsequently checked visually using the *3D KMaDG Vizualizátor SlimDX* app developed by Róbert Špir.

![MOP_3D_viz](https://github.com/Tomas-Homola/MOP_linBazFun3D_seriovyKod/assets/61438447/29a26798-012c-46c9-b135-6097fb2d8824)
