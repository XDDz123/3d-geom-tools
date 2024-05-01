# 3D Geometry Tools
This project implements an assortment of algorithms for processing 3D geometry.

## Iterative Closest Point (ICP)
The ICP algorithm could be described as follows: </br>
1. Select a number of points for matching.
2. Match every point to their cloest counterpart in the other point set.
3. Reject pairs based on distance or normal angles.
4. Minimize the error function such that after rotations and translations distance between two scans are minimal.
  
### Point-to-point ICP
Point-to-point ICP minimizes the error function: $\large E = \displaystyle \frac{1}{N_p} \displaystyle \sum_{i=1}^{N_p} |Rp_i + t - q_i|^2$. </br> 
This can be solved in closed form by the following: </br> 
1. Find the centers of mass of both scans: $\large \bar{p} = \displaystyle \sum_{i=1}^{N_p} p$ and $\large \bar{q} = \displaystyle \sum_{i=1}^{N_q} q$.  </br>
2. Subtract centroids from scans: $\large\tilde{p} = p_i - \bar{p}$ and $\large\tilde{q} = q_i - \bar{q}$.
3. Perform SVD on matrix $\large A = \displaystyle \sum_{i=1}^{N_p} \tilde{p}_i\tilde{q}_i^T$ and obtain matrices $\large U$, $\large V^t$. </br>
4. If $\large rank(A) = 3$, then $\large R = UV^T$ and $\large t = \bar{q} - R\bar{p}$. </br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/72bc49c0-980f-4392-933e-a86ea50f0422" width="490" height="420"> 

### Point-to-plane ICP
Point-to-point ICP minimizes the error function: $\large E = \displaystyle \frac{1}{N_p} \displaystyle \sum_{i=1}^{N_p} \[(Rp_i + t - q_i) \cdot n_i\]^2$. </br> 
This can be computed by solving the equation $\large Ax = b$ in closed form, where:  </br>
```math
\large
A = \begin{bmatrix}~
p_1 \times n_1 & n_1 \\
p_2 \times n_2 & n_2 \\
\vdots & \vdots
\end{bmatrix}
~~~~
x = \begin{bmatrix}~
r_x \\
r_y \\
r_z \\
t_x \\
t_y \\
t_z 
\end{bmatrix}
~~~~
b = \begin{bmatrix}~
-(p_1 - q_1) \cdot n_1 \\
-(p_2 - q_2) \cdot n_2 \\
\vdots
\end{bmatrix}
```
Then solve the least squares problem: $\large x = (A^TA)^{-1}A^Tb$.  </br></br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/eb727abb-0016-4e22-ad06-8c9f2ad00ae2" width="490" height="420"> 

## Uniform Laplace
### Mean Curvature
The uniform discretization is given by the following equation:</br>
```math
\large
\Delta_{uni} f(v_i) := \dfrac{1}{\mid\mathcal{N}_1(v_i)\mid}\sum\limits_{v_j \in \mathcal{N}_1(v_i)} (x_j - x_i) \approx -2Hn
```
The laplace operator $\large L$ can be computed by:</br>
```math
\large
L_{ij} = -1 \text{~if~} i \neq j, j \in \mathcal{N}_1(v_i) \\
```
```math
\large
L_{ij} = \mid\mathcal{N}_1(v_i)\mid \text{~if~} i = j
```
```math
\large
L_{ij} = 0 \text{~otherwise}
```
This means that we can construct the laplace operator through the
neighborhood information at each vertex. Specifically, the
number of neighbors of each vertex is recorded on the diagonal and the
indices of their corresponding neighbors are set to 1.

The mean curvature can then be computed by applying the laplace operator
on the mesh vertices:
```math
\large
H = \dfrac{\|\Delta_s\text{x}\|}{2}
```
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/f4dd1870-8b7d-4f9d-982b-64a037420a75" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/d1569844-d08a-4624-9c73-a08163c5ee89" width="50%" height="50%"> 

### Gaussian Curvature
Gaussian curvature can be computed from the angle deficit at each vertex. </br>

```math
\large
\text{angle deficit} = 2\pi - \sum\limits_{j} \theta_j
```
where $j$ is the angle at the current vertex in each of its connected triangles.</br>
In a perfectly flat region we would expect the angles to add up to $2\pi$. </br>

To obtain the gaussian curvature, we would need to normalize using Area $A$:</br>
```math
\large
K = (2\pi - \sum\limits_{j} \theta_j) / A
```

The form of area chosen for this implementation was barycentric cells,
where edge mid points and triangle barycenters are connected to form an area. This method was chosen for its simplicity, as the resulting area
is simply $1/3$ of the triangle areas. </br></br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/1e9435b0-f579-4be6-bef9-ef7e0494048d" width="50%" height="50%"> 

## Non-uniform Laplace (Discrete Laplace-Beltrami)
Laplace-Beltrami with the cotangent discretization is given by the following equation: </br>
```math
\large \Delta_{S} f(v_i) := \dfrac{1}{2A_i}\sum\limits_{v_j \in \mathcal{N}_1(v_i)} (cot\alpha_{ij}+cot\beta_{ij}) (f(v_j) - f(v_i))
```
In matrix form we can define the discrete laplace operator $\large L$ by: </br>
```math
\large L = M^{-1}C
```
```math
\large C_{i_j} = (cot\alpha_{ij} + cot\beta_{ij})/2 \text{~if~} i\neq j, j \in \mathcal{N}_1(v_i)
```
```math
\large C_{i_j} = -\sum_{v_j \in \mathcal{N}_1(v_i)} ((cot\alpha_{ij} + cot\beta_{ij})/2) \text{~if~} i=j
```
```math
\large C_{i_j} = 0 \text{~otherwise}
```
```math
\large M^{-1} = \text{diag}(...,\frac{1}{A_i},...)
```
where $\large \alpha_{ij}$ and $\large \beta_{ij}$ are the angles opposite to each edge connected to the current vertex.

At each neighbor of each vertex, we locate the connected faces that contain both vertices (e.g.~the current edge) and record the angles $\large \alpha_{ij}$ and $\large \beta_{ij}$ (at the other
vertex that is not the current pair), which eventually formulates the $\large C$ matrix. In addition, each connected face was accumulated in a vector to construct the $\large M$ matrix. </br></br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/6aaaabf2-7ee1-400b-a06e-87900c503f2a" width="50%" height="50%"> 

## Modal Analysis
Based on equation 4 in the paper [Spectral Geometry Processing with Manifold Harmonics](https://doi.org/10.1111/j.1467-8659.2008.01122.x) by Vallet and Lévy. </br>
The basis vectors can be computed as follows: 
```math 
\large \Delta \phi_i = \lambda_i \phi_i
```
```math 
\large M^{-1} C \phi_i = \lambda_i \phi_i
```
```math
\large M^{-1} C M^{-\frac{1}{2}} M^{\frac{1}{2}} \phi_i = \lambda_i \phi_i
```
```math
 \large M^{-\frac{1}{2}} M^{-\frac{1}{2}} C M^{-\frac{1}{2}} M^{\frac{1}{2}} \phi_i = \lambda_i \phi_i
```
```math
 \large M^{-\frac{1}{2}} C M^{-\frac{1}{2}} M^{\frac{1}{2}} \phi_i = \lambda_i M^{\frac{1}{2}} \phi_i
```
```math
 \large D = M^{-\frac{1}{2}} C M^{\frac{1}{2}}
```
```math
 \large \alpha_i = M^{\frac{1}{2}} \phi_i
```
```math
 \large D \alpha_i = \lambda_i \alpha_i
```
```math
 \large \phi_i = M^{-1/2} \alpha_i
```
We first find the $\large k$ smallest eigen vectors of $\large D$.
Then obtain the basis vectors $\large \phi_i$ by mapping them into canonical basis (multiplying by $\large M^{-1/2}$).

The reconstruction for each dimension (of the vertices) can be computed by:
```math
\large x := [x_1, ..., x_n]
```
```math
\large x \leftarrow \sum\limits_{i=1}^{k}(x^T\phi_i)\phi_i
```
## Mesh Smoothing
### Explicit Laplacian Mesh Smoothing
Explicit smoothing is computed with the following equation:
```math
\large P^{(t+1)} = (I + \lambda L) P^{(t)}
```
where the Laplace-Beltrami operator is applied to the given mesh in
small steps as defined by $\large \lambda$. <\br> 
This process aims to slowly smooth out the given mesh/surface, however it is conditionally stable depending on the chosen magnitude of $\large \lambda$.

### Implicit Laplacian Mesh Smoothing
In contrast to explicit smoothing, implicit smoothing is unconditionally stable. <\br>
Implicit smoothing is computed with the following equations:
```math
\large (I - \lambda L) P^{(t+1)} = P^{(t)}
```
Since $\large L = M^{-1} C$ is no longer symmetric after the normalization
from \(M\), we must first symmetrize
```math
\large (M - \lambda C) P^{(t+1)} = M P^{(t)}
```
the resulting sparse system is symmetric positive definite and can be solved using methods such as iterative conjugate gradients.
