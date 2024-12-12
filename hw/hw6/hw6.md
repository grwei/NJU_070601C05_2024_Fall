---
export_on_save:
  html: true
html:
  embed_local_images: true
  embed_svg: true
  offline: false
  toc: true
print_background: true
---

# Homework 6

- **Course:** Numerical Solutions to PDEs - FALL 2024
- **Instructor:** Zhou, Bowen ([周博闻](https://as.nju.edu.cn/54/79/c11339a218233/page.htm))
- **Due date:** Nov. 22, 2024
- **Submitted date:** Dec. DD, 2024
- **Problem set:** [PS6.pdf](https://box.nju.edu.cn/d/439906db314e411489a3/files/?p=%2FProblemSets%2FPS6.pdf)
- **Course website:** <https://grwei.github.io/NJU_070601C05_2024_Fall/>

> &ensp; &emsp; Describe the setup and each step in your solutions with words and clearly label your final answers. Use Matlab for plotting and programming and include your code as an appendix to your problem set.

## Table of Contents {ignore=true}

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [Homework 6](#homework-6)
  - [Problem 1](#problem-1)
  - [Problem 2](#problem-2)
  - [Acknowledgement](#acknowledgement)
  - [Contact Information](#contact-information)
  - [Appendix](#appendix)
    - [Matlab code for problem 1](#matlab-code-for-problem-1)

<!-- /code_chunk_output -->

## Problem 1

> &ensp; &emsp; Poisson equations appear often in fluid mechanics, e.g. in potential flow calculations, which include groundwater flow and very high Reynolds number situations where viscous effects are neglected. The Poisson equation also appears when solving for the pressure in incompressible flow Navier-Stokes calculations, and solving it actually accounts for much of the computational cost (as we will see later). Hence, there is a lot of interest in developing fast methods to solve linear systems of equations that arise from elliptic equations.
> &ensp; &emsp; Here, you will test out the basic ideas of iterative methods by applying the point Jacobi, Gauss-Seidel, and successive over-relaxation algorithms to the solution of the steady-state Laplace equation describing steady flow in square open channel. Consider gravity-driven flow in a square open channel, with 3 side walls and the top open to the atmosphere. The full flow equations are given by
>
> $$\frac{\partial u_i}{\partial t} + \frac{\partial u_i u_j}{\partial x_j} = - \frac{1}{\rho} \frac{\partial P}{\partial x_i} + \nu \frac{\partial^2 u_i}{\partial x_j x_j} - g \delta_{i3},$$
>
> $$\frac{\partial u_i}{\partial x_i} = 0.$$
>
> &ensp; &emsp; Assuming steady-state and fully-developed flow (no variation in x), the equations reduce to a balance between bottom friction and gravity:
>
> $$-g \sin{\theta} = \nu \left( \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 y}{\partial z^2} \right),$$
>
> where $u$ is the horizontal velocity, $g$ is gravitational acceleration, $\nu$ is the viscosity, and $x$ is streamwise, $y$ is spanwise, and $z$ is the vertical coordinate. Boundary conditions are no slip on the three side walls ($y = −L, \, y = L$ and $z = 0$) and free slip on the free surface ($z = H$), where $L = 1 \, \text{m}, H = 1 \, \text{m}$. The bottom slope angle is $\theta = 0.08$ and $\nu = 10^{−6} \, \mathrm{m}^2 / \mathrm{s}$.
>
> &ensp; &emsp; You will solve for the flow in a two-dimensional cross section as shown in the figure below. A computational grid with $\Delta y = \Delta z = 0.1 \, \mathrm{m}$ should be used; then repeat the calculations with spacings of $0.05 \, \mathrm{m}$ and discuss the differences as you answer each question below.
>
> 1. Write a program to compute the steady-state solution using the Jacobi iteration method. Discretize the equation using a second-order scheme. State your equation for $u$ and demonstrate how boundary conditions will be applied. Use one-sided differences for the no-flux boundary condition. Initialize your array with zeros. Continue the iterations until the maximum difference between successive iterations is less than $0.001$.
> 2. You can monitor the progress of the solution by plotting the value of the solution at the center of the domain: $(y, z) = (0, 0.5)$. How many iterations are required until the solution at this point steadily varies by no more than $0.001$ between iterations?
> 3. Plot a contour map of the velocity field in the channel cross section. Also plot vertical profiles ($u(z)$) at the channel midpoint and one-fourth of the way from the edge.
> 4. Repeat the above using the Gauss-Seidel iteration method and then with successive over-relaxation (SOR). Find the optimal value of $\omega$ (to within $0.05$) by trial and error. Show a plot that compares the convergence of the methods, and discuss the performance.

&ensp; &emsp; 本节讨论二维矩形区域上的 Poisson 方程

$$
\begin{equation*}
  \tag{1.1}
  \left\{
  \begin{aligned}
    & u_{yy} + u_{zz} = - f_0, \quad y \in [0, 2L], \, z \in [0, H], \\
    & \left. u \right|_{y = 0, \, 2L} = \left. u \right|_{z = 0} = \left. u_z \right|_{z = H} = 0
  \end{aligned}
  \right.
\end{equation*}
$$

的数值解, 式中

$$
\begin{equation*}
  \tag{1.2}
  f_0 := (g / \nu) \sin{\theta}
\end{equation*}
$$

为常数.

&ensp; &emsp; 问题 (1.1) 的解析解可用分离变量法求出. 注意到

$$
\begin{equation*}
  \tag{1.3}
  w :=u + y(y - L) f_0 / 2
\end{equation*}
$$

## Problem 2

> &ensp; &emsp; Consider the 1D boundary value problem
>
> $$
> \begin{align*}
>   & \frac{\mathrm{d}^2 \phi}{\mathrm{d}x^2} = \sin{(k \pi x)}, \quad 0 \le x \le 1, \\
>   & \phi(0) = \phi(1) = 0,
> \end{align*}
> $$
>
> use centered difference scheme on a uniform mesh
>
> $$
> \begin{align*}
>   \frac{\phi_{j+1} - 2\phi_j + \phi_{j-1}}{\Delta^2} & = \sin{(k \pi j \Delta)}, \quad j = 1, 2, \cdots, N - 1, \\
>   \phi_0 &= \phi_N = 0,
> \end{align*}
> $$
>
> start with an initial guess of $\phi^{[0]} = 0$, and iterate with Gauss-Seidel scheme
>
> $$
> \begin{equation*}
>   \phi_j^{[k+1]} = \frac{1}{2} \left[ \phi_{j+1}^{[k]} + \phi_{j-1}^{[k+1]} - \Delta^2 \sin{(k \pi j \Delta)} \right].
> \end{equation*}
> $$
>
> 1. Use $N = 64$ grid points, for each $k = 1,2,4,8,16$, plot the convergence curve, with the number of iterations on the $x$-axis, and the maximum of the absolute value of the residual $\max{(|r_j|)}$ on the $y$ axis. The residual
> $$ r_j^{[k]} := \sin{(k \pi x_j)} - \left. \left( \phi_{j+1}^{[k]} - 2 \phi_{j}^{[k]} + \phi_{j-1}^{[k]} \right) \right/ \Delta^2.$$
> 2. Now consider a slightly more complicated right-hand side for the boundary value problem
> 3. Use
>

&ensp; &emsp; text

## Acknowledgement

&ensp; &emsp; I am grateful to ...

## Contact Information

- **Author:** Guorui Wei (危国锐)
- **E-mail:** [313017602@qq.com](mailto:313017602@qq.com)
- **Website:** <https://github.com/grwei>

## Appendix

### Matlab code for problem 1

```matlab {.line-numbers}
%
```

徐老师, 目前
