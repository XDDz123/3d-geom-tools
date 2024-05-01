# 3D Geometry Tools
This project implements an assortment of algorithms for processing 3D geometry.

## ICP
The ICP algorithm could be described as follows: </br>
1. Select a number of points for matching
2. Match every point to their cloest counterpart in the other point set.
3. Reject pairs based on distance or normal angles
4. Minimize the error function such that after rotations and translations, the distance between the two scans are minimal
  
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

