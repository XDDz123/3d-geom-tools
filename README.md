# 3D Geometry Tools
This project implements an assortment of algorithms for processing 3D geometry:
* [Iterative Closest Point](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#iterative-closest-point-icp)
* [Uniform Laplace](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#uniform-laplace)
* [Non-uniform Laplace (Discrete Laplace-Beltrami)](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#non-uniform-laplace-discrete-laplace-beltrami)
* [Modal Analysis](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#modal-analysis)
* [Mesh Smoothing](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#mesh-smoothing)
* [Geodesics](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#geodesics-in-heat)
  * [Path finding](https://github.com/XDDz123/3d-geom-tools?tab=readme-ov-file#path-finding)

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
L_{ij} = \begin{cases}~
-1 & \text{~if~} i \neq j, j \in \mathcal{N}_1(v_i) \\
\mid\mathcal{N}_1(v_i)\mid & \text{~if~} i = j \\
~0 & \text{~otherwise}
\end{cases}
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
\large C_{i_j} = \begin{cases}
~(cot\alpha_{ij} + cot\beta_{ij})/2 & \text{~if~} i\neq j, j \in \mathcal{N}_1(v_i) \\
-\sum_{v_j \in \mathcal{N}_1(v_i)} ((cot\alpha_{ij} + cot\beta_{ij})/2) & \text{~if~} i=j \\
~0 & \text{~otherwise}
\end{cases}
```
```math
\large M^{-1} = \text{diag}(...,\frac{1}{A_i},...)
```
where $\large \alpha_{ij}$ and $\large \beta_{ij}$ are the angles opposite to each edge connected to the current vertex.

At each neighbor of each vertex, we locate the connected faces that contain both vertices (e.g. the current edge) and record the angles $\large \alpha_{ij}$ and $\large \beta_{ij}$ (at the other
vertex that is not the current pair), which eventually formulates the $\large C$ matrix. In addition, each connected face was accumulated in a vector to construct the $\large M$ matrix. </br></br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/6aaaabf2-7ee1-400b-a06e-87900c503f2a" width="50%" height="50%"> 

## Modal Analysis
Based on equation 4 in the paper [Spectral Geometry Processing with Manifold Harmonics](https://doi.org/10.1111/j.1467-8659.2008.01122.x) by Vallet and Lévy. </br>
The basis vectors can be computed from the following derivation: 
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
#### Sample
The following reconstructs the model using $k$ eigen vectors. </br></br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/1d530eeb-f341-4b60-a5ca-fb23ac25f0e4" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/263e1107-2ff9-4b3d-b32f-f380fbc687a0" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/6792d1aa-0a55-4e90-8e01-72bda70ee3ce" width="50%" height="50%"> 

## Mesh Smoothing
### Explicit Laplacian Mesh Smoothing
Explicit smoothing is computed with the following equation:
```math
\large P^{(t+1)} = (I + \lambda L) P^{(t)}
```
where the Laplace-Beltrami operator is applied to the given mesh in
small steps as defined by $\large \lambda$. </br> 
This process aims to slowly smooth out the given mesh/surface, however it is conditionally stable depending on the chosen magnitude of $\large \lambda$.

### Implicit Laplacian Mesh Smoothing
In contrast to explicit smoothing, implicit smoothing is unconditionally stable.
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

#### Output
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/5d29c7df-efcd-4eca-bfe1-a4446d83546c" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/b1394aee-87a0-44ba-8ac3-3451b807fec1" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/9171ba8f-4569-40dc-b432-be27a024ec87" width="50%" height="50%"> 

## Geodesics in Heat
This section implements techniques described in the paper [Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow](https://ddg.math.uni-goettingen.de/pub/GeodesicsInHeat.pdf) by Crane, Weischedel, and Wardetzky. </br>
In general terms, the aim is to compute the distance between vertices on a surface mesh or manifold by calculating the geodesic distance based on the physical mechanics of heat flow.  </br></br>
The outline of the heat method can described as follows: </br>
1. Standard geodesic distance is computed by solving the [eikonal equation](https://en.wikipedia.org/wiki/Eikonal_equation) $\large |\triangledown\phi| = 1$, where $\phi$ is the geodesic distance. </br>
   However, most processes involved in computing this measure are often non-linear.</br></br>
2. In order to solve for $\large \phi$, we take advantage between the relationship between heat distribution and distance using Varadhan’s formula:</br>
 ```math
 \large \displaystyle{\phi(x, y) = \lim_{t \to 0} \sqrt{-4t\text{log}k_{t,x}(y)}} \text{, where } k_{t,x} \text{ is the heat kernel.}
 ```
3. Reconstructions of $\large k_{t,x}$ are found to be extremely sensitive to numerical errors. </br> The paper cirumvents this problem with the following procedure: </br>
   * Given a reconstruction of the heat kernel $\large u$
   * Compute its gradient $\large -\triangledown u$
   * Then normalize such that the output gradient field $\large X = -\displaystyle{\frac{\triangledown u }{\|\triangledown u\|}}$ is only dependent on the direction of the gradient and robust to errors in magnitude. </br></br>
4. Lastly, $\large \phi$ can be computed by minimizing the equation $\large \|\triangledown \phi - X\|^2$. This minimization problem can be solved linearly with $\large \Delta \phi = \triangledown X$.
### Implementation
The geodesic distance on manifolds can be computed with the following steps:
1. Compute the cotangent Laplacian matrix $\large L_C$ and the weight (area) matrix $\large A$. </br>
* The discretization of the Laplacian can be computed from:
```math
\large 
(Lu)_i=\frac{1}{2A_i}\displaystyle{\sum_j(cota_{ij}+cot\beta_{ij})(u_j-u_i)}
```
```math
\text{Where } A_i = \frac{\text{Area of one ring neighbors at i}}{3}, \alpha_{ij} ~\&~ \beta_{ij} \text{~are the angles opposite to edge}_{ij}
```
* In matrix form we can define the discrete Laplacian operator with
```math
\large 
L_{C_{ij}} = \begin{cases}
(cot\alpha_{ij} + cot\beta_{ij})/2, & \text{if  \(i\neq j, j \in \mathcal{N_1}(v_i)\)} \\
-\sum_{v_j \in \mathcal{N}_1(v_i)} ((cot\alpha_{ij} + cot\beta_{ij})/2), & \text{if \(i=j\)} \\
0, & \text{otherwise}
\end{cases}
```
* For Dirichlet boundary conditions, boundary vertices have $\large L_{C_{ij}} = 0$ for $\large j \in {\mathcal{N_1}}(v_i)$. </br>
2. Compute the optimal time $\large t = mh^2$, where $\large m=1$ and $\large h$ is the mean distance between all vertices.</br>
3. Compute the heat flow $\large u$ with $\large (A - tL_C)u = \delta_\gamma$, where the Kronecker delta $\large \delta_\gamma$ is a vector where heat source vertices are set to 1 and other vertices to 0.</br>
4. Compute the gradient $\triangledown u$ with $\large \triangledown u = \frac{1}{2A_f}\displaystyle{\sum_i u_i(N \times e_i)}$, where $\large A_f$ is the face area, $\large N_i$ is the face normal, and $\large e_i$ is the vector (edge) opposite of the current vertex orientated counter-clockwise.</br>
5. Compute the normalized vector field $\large X = \triangledown u / \| \triangledown u\|$.</br>
6. Compute the divergences of $\large X$ with $\large \triangledown\cdot X = \frac{1}{2}\displaystyle{\sum_j cot\theta_1(e_1 \cdot X_j) + cot\theta_2 (e_2 \cdot X_j)}$, where $\large \theta_1$ and $\large \theta_2$ are the 2 remaining vertex angles of the current face $\large j$, $\large e_1$ and $e\large _2$ are the edges opposite to $\large \theta_1$ and $\large \theta_2$ respectively.</br>
7. Compute the geodesics with $\large L_c\phi=\triangledown\cdot X$ by solving the linear system.
  
### Results
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/5e2f1da8-8a97-4eb8-ae3f-781c5d1565c2" width="50%" height="50%">
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/26e47ece-b692-4af9-b2a1-84416978af33" width="50%" height="50%"> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/af5a7f11-998a-4fa3-a096-d86e5504cdc3" width="50%" height="50%"> 

### Path finding
A path finding algorithm was implemented to compute the shortest path between two vertices on a manifold. </br> </br> 
Given a start vertex and an end vertex, the geodesic distance at the end is computed. </br> 
Starting at the start vertex, the algorithm moves the current vertex to the neighbor that is closest to the end vertex until current vertex reaches the designated destination. </br> </br> 
In the samples below, the blue line marks the shortest path. </br>  </br> 
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/2f30d69a-ffbd-4f8b-8962-0b0dc16d030e" width="50%" height="50%">
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/6027d04d-592b-4f83-b8d0-1fffb6aab5a7" width="50%" height="50%">
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/a282f65d-0c26-4988-a5b8-cdd68c67b960" width="50%" height="50%">

### Geodesics under noise
Gaussian noise was added to the mesh vertices, where the $\sigma$ of the gaussian distribution is based on the bounding box size of the mesh in each dimension scaled by a factor $k$. </br>
As expected, the algorithm is robust. Similar to the results of the paper, the computed geodesic distances are reasonable even when the amount of noise is relatively large.  </br> </br>
<img src="https://github.com/XDDz123/3d-geom-tools/assets/20507222/9c7aae53-7d8a-492c-ad31-028a1dae1a1d" width="50%" height="50%">
