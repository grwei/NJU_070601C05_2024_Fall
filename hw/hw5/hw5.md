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

# Homework 5

- **Course:** Numerical Solutions to PDEs - FALL 2024
- **Instructor:** Zhou, Bowen ([周博闻](https://as.nju.edu.cn/54/79/c11339a218233/page.htm))
- **Due date:** Nov. 15, 2024
- **Submitted date:** Dec. 4, 2024
- **Problem set:** [PS5.pdf](https://box.nju.edu.cn/d/439906db314e411489a3/files/?p=%2FProblemSets%2FPS5.pdf)
- **Course website:** <https://grwei.github.io/NJU_070601C05_2024_Fall/>

> &ensp; &emsp; Describe the setup and each step in your solutions with words and clearly label your final answers. Use Matlab for plotting and programming and include your code as an appendix to your problem set.

## Table of Contents {ignore=true}

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [Homework 5](#homework-5)
  - [Table of Contents {ignore=true}](#table-of-contents-ignoretrue)
  - [Problem 1](#problem-1)
    - [Problem 1.1](#problem-11)
    - [Problem 1.2](#problem-12)
    - [Problem 1.3](#problem-13)
    - [Problem 1.4](#problem-14)
    - [Problem 1.5](#problem-15)
    - [Problem 1.6](#problem-16)
    - [Problem 1.7](#problem-17)
    - [Problem 1.8](#problem-18)
  - [Acknowledgement](#acknowledgement)
  - [Contact Information](#contact-information)
  - [Appendix](#appendix)
    - [Solid wall B.C. 对 Fick 扩散的约束 (ChatGPT mini 4o)](#solid-wall-bc-对-fick-扩散的约束-chatgpt-mini-4o)
      - [Physical Interpretation](#physical-interpretation)
      - [Mathematical Interpretation](#mathematical-interpretation)
      - [Summary](#summary)
    - [Matlab code for problem 1](#matlab-code-for-problem-1)
      - [hw5\_1.m](#hw5_1m)
      - [DoSag2D.m](#dosag2dm)
      - [bSrc.m](#bsrcm)
      - [getSolution.m](#getsolutionm)
      - [prepareParams.m](#prepareparamsm)
      - [reset.m](#resetm)
      - [solve2D.m](#solve2dm)
      - [solve2DSteady.m](#solve2dsteadym)

<!-- /code_chunk_output -->

## Problem 1

> &ensp; &emsp; **Streeter–Phelps Model.** The Streeter–Phelps equation is used in the study of water pollution as a water quality modeling tool. The model describes how dissolved oxygen (DO) decreases in a river or stream along a certain distance by degradation of biochemical oxygen demand (BOD). The equation was derived by H. W. Streeter, a sanitary engineer, and Earle B. Phelps, a consultant for the U.S. Public Health Service, in 1925, based on field data from the Ohio River. The equation is also known as the DO sag equation.
>
> &ensp; &emsp; In this problem set, we consider the Streeter-Phelps model in a 2D river. The governing equations for BOD and DO includes advection and diffusion process, with additional source and sink terms.
>
> $$
> \begin{equation*}
> \tag{1}
> \frac{\partial B}{\partial t} + u_i \frac{\partial B}{\partial x_i} = D_i \frac{\partial^2 B}{\partial x_i^2} - K_{\text{r}} B + \delta(x - x_0, y - y_0) Q_{\text{B}},
> \end{equation*}
> $$
>
> $$
> \begin{equation*}
> \tag{2}
> \frac{\partial O}{\partial t} + u_i \frac{\partial O}{\partial x_i} = D_i \frac{\partial^2 O}{\partial x_i^2} + K_{\text{a}} (O_{\text{sat}} - O) - K_{\text{r}} B,
> \end{equation*}
> $$
>
> where $B$ and $O$ are BOD and DO concentrations with units $[\mathrm{M} \mathrm{L}^{−3}]$. $K_{\text{r}}$ $[\mathrm{T}^{−1}]$ is the reaction rate between BOD and DO. $K_{\text{a}}$ $[\mathrm{T}^{−1}]$ is the mass transfer coefficient for oxygen. The reaeration rate is given by $K_{\text{a}} (O_{\text{sat}} −O)$, where $O_{\text{sat}}$ is the saturation oxygen content in water exposed to air, it is mostly affected by water temperature. The reaeration term describes the flux of oxygen across the air-water interface into the water. $Q_{\text{B}}$ $[\mathrm{M} \mathrm{T}^{−1}]$ is the rate of BOD discharge into the river. $\delta(x −x_0, y −y_0)$ is a delta function indicating a point source. $i = 1, 2$ for the $x$ and $y$ directions. Repeated index suggest summation over each value of the repeated index, e.g.,
>
> $$
> \begin{equation*}
>   u_i \frac{\partial B}{\partial x_i} \equiv u_1 \frac{\partial B}{\partial x_1} + u_2 \frac{\partial B}{\partial x_2} = u \frac{\partial B}{\partial x} + v \frac{\partial B}{\partial y}.
> \end{equation*}
> $$
>
> &ensp; &emsp; **1.** Discretize the above two equations using the forward in time, centered in space (FTCS) scheme. Set the proper lateral boundary conditions for the problem (HINT: use solid wall boundaries in the cross flow $y$ direction, set Neumann (zero-gradient) conditions for the outflow boundary, and Dirichlet conditions for the inflow boundary. Try using the ghost point technique for setting boundary conditions).
>
> &ensp; &emsp; **2.** Perform von Neumann stability analysis for the discretized equations. Determine the proper proper Courant number $\text{Cr} = u \Delta t / \Delta x, \, r = D \Delta t / {\Delta x}^2$, and Peclet number $\text{Pe} = u \Delta x / D$ values for a stable solution (HINT: you may combine the two equations, by ignoring the BOD discharge term $\delta(x −x_0, y − y_0) Q_{\text{B}}$ in Eq. (1), and the $K_{\text{a}} O_{\text{sat}} − K_{\text{r}} B$ terms in Eq. (2). You may further simply the analysis by assuming $D_x = D_y$).
>
> &ensp; &emsp; **3.** Perform modified equation analysis, and discuss the truncation errors.
>
> &ensp; &emsp; **4.** Write a Matlab script to solve the problem. Use the following assumptions:
>
> - The length and width of the river is $L_x = 150 \, \text{m}, \, L_y = 30 \, \text{m}$.
> - The flow is assumed to be in the $x$-direction, i.e. $u = 0.4 \, \text{m} \, \text{s}^{−1}, v = 0$.
> - The diffusion coefficients $D_x = D_y = 0.5 \, \text{m}^2 \, \text{s}^{−1}$.
> - The reaction rate $K_{\text{r}} = 0.01 \, \text{s}^{−1}$.
> - The mass transfer coefficient for oxygen $K_{\text{a}} = 0.02 \, \text{s}^{−1}$.
> - The saturation oxygen content $O_{\text{sat}} = 8 \, \text{mg} \; \text{L}^{−1}$.
> - The BOD discharge rate $Q_{\text{B}} = 70 \, \text{mg} \; \text{s}^{−1}$.
> - Place the point source $Q_{\text{B}}$ at $x = 10 \, \text{m}, \, y = 15 \, \text{m}$.
> - Use $\Delta x = \Delta y = 1 \, \text{m}$ for the homework. (You may use finer spacings).
>
> &ensp; &emsp; **5.** Plot a $x, y$ contour of the steady state (i.e. $\partial / \partial t = 0$) solution for both BOD and DO. (HINT: use Matlab tool: `contour` or `contourf`). Plot the DO sag curve and BOD curve along the center line $y =15 \, \text{m}$, i.e. plot $B(x, y = 15 \, \text{m})$ and $O(x, y = 15 \, \text{m})$.
>
> &ensp; &emsp; **6.** According to Environmental Protection Agency (EPA) requirement, the minimum steady state DO in a river must be above $6 \, \text{mg} \; \text{L}^{−1}$. Does the current discharge satisfy the EPA requirement? If not, propose alternative discharge methods, test with your model, and present your results. Note that you must maintain the total discharge rate at $70 \, \text{mg} \; \text{s}^{−1}$.
>
> &ensp; &emsp; **7.** Suggest 3 or more ways to make your code more realisitic.
>
> &ensp; &emsp;  **8.** Suggest 3 or more ways to speed up your code.

&ensp; &emsp; 本章讨论二维 Streeter-Phelps 模型 (1)(2) 的数值解. 用 forward in time, centered in space (FTCS) 方案, 将 Eq. (1) 离散为

$$
\begin{equation*}
    \tag{3a}
    \begin{aligned}
        & \quad \; \frac{1}{\tau} \left( B^{n+1}_{i,j} - B^{n}_{i,j} \right) + \frac{u^{n}_{i, j}}{2 h_x} \left( B^{n}_{i+1,j} - B^{n}_{i-1,j} \right) + \frac{v^n_{i,j}}{2 h_y} \left( B^{n}_{i,j+1} - B^{n}_{i,j-1} \right) \\
        & = \frac{D_x}{h^2_x} \left( B^{n}_{i+1,j} - 2 B^{n}_{i,j} + B^{n}_{i-1,j} \right) + \frac{D_y}{h^2_y} \left( B^{n}_{i,j+1} - 2 B^{n}_{i,j} + B^{n}_{i,j-1} \right) - K_{\text{r}} B^n_{i,j} + \frac{Q_\text{B}}{h_x h_y} \mathrm{I}_{i,j}(x_0, y_0),
    \end{aligned}
    \end{equation*}
$$

where

$$
n = 0, 1, \cdots, N, \quad i = 0, 1, \cdots I, \quad j = 0, 1, \cdots, J,
$$

示性函数

$$
\begin{equation*}
    \tag{4}
    \mathrm{I}_{i,j}(x_0, y_0) := \left\{
        \begin{aligned}
            & 1, \quad 源点 (x_0, y_0) \text{ 落在 $(i, j)$ 格点所代表的那个有限体积元}, \\
            & 0, \quad \text{otherwise}.
        \end{aligned}
    \right.
\end{equation*}
$$

可以设置多个源点. 将 Eq. (2) 离散为

$$
\begin{equation*}
    \tag{3b}
    \begin{aligned}
        & \quad \; \frac{1}{\tau} \left( O^{n+1}_{i,j} - O^{n}_{i,j} \right) + \frac{u^{n}_{i, j}}{2 h_x} \left( O^{n}_{i+1,j} - O^{n}_{i-1,j} \right) + \frac{v^n_{i,j}}{2 h_y} \left( O^{n}_{i,j+1} - O^{n}_{i,j-1} \right) \\
        & = \frac{D_x}{h^2_x} \left( O^{n}_{i+1,j} - 2 O^{n}_{i,j} + O^{n}_{i-1,j} \right) + \frac{D_y}{h^2_y} \left( O^{n}_{i,j+1} - 2 O^{n}_{i,j} + O^{n}_{i,j-1} \right) + K_{\text{a}} \left( O_{\text{sat}} - O^n_{i,j} \right) - K_{\text{r}} B^n_{i,j}.
    \end{aligned}
    \end{equation*}
$$

格式 (3) 中的速度场 $u^n, v^n$ 被假定已由数值模式的动力框架给出.

### Problem 1.1

&ensp; &emsp; 我们要在二维有界区域

$$
\begin{equation*}
    \tag{5}
    t \ge 0, \quad (x, y) \in \Omega := [0, L_x] \times [0, L_y]
\end{equation*}
$$

上求解模型 (1)(2), 为此, 须为模型设置合适的初始条件和边界条件, 这里取为

$$
\begin{equation*}
    \tag{6}
    \left\{
    \begin{aligned}
        & (\text{I.C.}) \quad
        \begin{aligned}
            % \tag{6}
            \left. F \right|_{t = 0} = F^0(x, y), \quad (x,y) \in \Omega,
        \end{aligned} \\
        & (\text{B.C.}) \quad
        \left\{
        \begin{aligned}
            % \tag{7}
            & (\text{Dirichlet}) & \quad \left. F \right|_{x = 0} = F_0(y, t), \\
            & (\text{Neumann}) & \quad \left. F_x \right|_{x = L_x} = F_1(y, t), \\
            & (\text{solid wall}) & \quad \left( \boldsymbol{n} \cdot \boldsymbol{u} \right)_{y = 0, \, L_y} = 0, \\
            & ({\color{blue}{\text{solid wall 2}}}) & \quad \left[ \boldsymbol{n} \cdot (- \kappa \nabla F) \right]_{y = 0, \, L_y} = 0.
        \end{aligned}
        \right. \\
    \end{aligned}
    \right.
\end{equation*}
$$

Eq. (6) 中, $F$ 指代 $B$ 和 $O$.

&ensp; &emsp; 相应地, 为格式 (3) 引入由

$$
\begin{align*}
    \tag{7a}
    & \frac{F^n_{-1, j} + F^n_{1, j}}{2} = F_0(t_n, y_j), \\
    \tag{7b}
    & \frac{F^n_{1, j} - F^n_{-1, j}}{2 h_x} = F_1(t_n, y_j), \\
    \tag{8a}
    & \color{gray}{F^n_{i, -1} = 3 F^n_{i, 0} - 3F^n_{i, 1} + F^n_{i, 2}}, \\
    \tag{8b}
    & \color{gray}{F^n_{i, J+1} = 3 F^n_{i, J} - 3 F^n_{i, J-1} + F^n_{i, J-2}},
\end{align*}
$$

决定的 halo points $F^n_{-1, j}, F^n_{I + 1, j}, F^n_{i, -1}, F^n_{i, J + 1}$. <span style="color: gray">其中, Eq. (8) 的设计考虑是使格式 (3) 在网格边界 $j = 0, J$ 处能以 $\mathcal{O}(h_y)$ 精度离散模型 (1)(2) 在边界 $y = 0, L_y$ 处的 Fick 扩散项. Halo points $F^n_{i, -1}, F^n_{i, J+1}$ 的设计 (8) 仅考虑 Fick 扩散, 是因为 (6) 中 B.C. 第二项 (y 边界 solid wall) 保证 $v^n_{i, 0} = v^n_{i, J} = 0$, 所以无论 $F^n_{i, -1}, F^n_{i, J+1}$ 如何选取, 都不会破坏格式 (3) 对模型 (1)(2) 在边界 $y = 0, L_y$ 处 $v$-平流 (*exactly* zero) 的离散</span>. <b style="color: blue;">一个重要问题是, 在物理上, solid wall 边界对 Fick 扩散的约束是怎样的? 若是梯度扩散通量为零, 即 $\color{black}{ \boldsymbol{n} \cdot (- \kappa \nabla F) = 0}$, 则与齐次 Neumann 边界 $\color{black}{\partial F / \partial n = 0}$ 等效?</b> 发现 ChatGPT [支持](#solid-wall-bc-对-fick-扩散的约束-chatgpt-mini-4o)此观点, 且有符合常识的物理解释, 所以我将 (8) 改为

$$
\begin{equation*}
    \tag{$8'$}
    F^n_{i, -1} = F^n_{i, 1}, \quad F^n_{i, J+1} = F^n_{i, J-1}.
\end{equation*}
$$

<span style="color: blue;">但, 若我坚持用 (8), 则会引入多大的误差呢?</span> **毕竟, Eq. (8), 尽管不蕴含 (8'), 似乎也是在合理前提下** ($\partial^2 F / \partial n^2$ 连续到边界? 这个条件过于苛刻吗? 零通量边界会引起二阶导突变甚至不存在?) **导出的?**

&ensp; &emsp; 引入 halo points (7)(8') 后, 在速度场 $u^n_{i,j}, v^n_{i,j}$ 已获得的前提下 (速度场的求解方案 is beyond the scope of this article), 只需根据模型 (1)(2) 的初始条件 (6) 制作初始场

$$
\begin{equation*}
    \tag{9}
    B^0_{i,j} = B^0(x_i, y_j), \quad O^0_{i,j} = O^0(x_i, y_j), \quad x_i = i h_x, \, y_j = j h_y,
\end{equation*}
$$

就可驱动格式 (3) 向前积分. 以后, 简单起见, 将速度场取为常向量, 即认为 $u^n \equiv u, \, v^n \equiv v$, where $u, v$ 是已知的常数.

### Problem 1.2

&ensp; &emsp; 本节对格式 (3) 作 von Neumann 稳定性分析. 为简化问题, I ignored the BOD discharge term $\delta(x −x_0, y − y_0) Q_{\text{B}}$ in Eq. (1), and the $K_{\text{a}} O_{\text{sat}} − K_{\text{r}} B$ terms in Eq. (2). 简化后的格式 (3a) 和 (3b) 有相同形式的增长因子,

$$
\begin{equation*}
    \tag{10}
    \begin{aligned}
        G & = (1 - 2 r_x - 2 r_y - K \tau) \\
        & + 2 r_x \cos{(\xi_x h_x)} - \mathrm{i} \, \mathrm{Cr}_x \, \sin{(\xi_x h_x)} \\
        & + 2 r_y \cos{(\xi_y h_y)} - \mathrm{i} \, \mathrm{Cr}_y \, \sin{(\xi_y h_y)},
    \end{aligned}
\end{equation*}
$$

其中 $K$ 指代 $K_{\mathrm{r}}$ 和 $K_{\mathrm{a}}$,

$$
\begin{equation*}
    \mathrm{Cr}_i := u_i \tau / h_i^2, \; r_i := D \tau / h_i^2, \; P_{\mathrm{e}, i} := u_i h / D,
\end{equation*}
$$

下标 $i = 1, 2, 3$ 分别表示 $x, y, z$ 方向. 要使格式稳定, 就是要选取 $\tau, h_x, h_y$, 使加强的 von Neumann 条件 ([Hw 3](../hw3/hw3.md#problem-1) 的 Eq. (1.22))

$$
\begin{equation*}
    \tag{11}
    |G| \le 1, \quad \forall |\xi_x h_x| \le \pi, \, |\xi_y h_y| \le \pi
\end{equation*}
$$

被满足.

&ensp; &emsp; 设 $A$ 是定义了加法运算的集合, $A_1, A_2 \subset A$. 定义集合间的加法

$$
\begin{equation*}
    \tag{12}
    A_1 + A_2 := \left\{ \left. a_1 + a_2 \; \right| \; a_i \in A_i \right\}.
\end{equation*}
$$

记平面直角坐标系上的椭圆

$$
\begin{equation*}
    \tag{13}
    C(a, b) := \left\{(x, y) \; \left| \; \frac{x^2}{a^2} + \frac{y^2}{b^2} = 1 \right. \right\}, \quad a, b > 0.
\end{equation*}
$$

当 $(\xi_x h_x, \xi_y h_y)$ 取遍 $[-\pi, \pi)^2$, 由 (10) 知 $G = x + \mathrm{i} y$ 在复平面上的轨迹是

$$
\begin{equation*}
    \tag{14}
    z_0 + C(2 r_x, |\mathrm{Cr}_x|) + C(2 r_y, |\mathrm{Cr}_y|),
\end{equation*}
$$

其中

$$
\begin{equation*}
    \tag{15}
    z_0 := 1 - 2 r_x - 2 r_y - K \tau
\end{equation*}
$$

在实轴上. 由 (14) 看出, 命题 (11) 的一个充分条件是

$$
\begin{equation*}
    \tag{16}
    |z_0| + \max{\{2 r_x, |\mathrm{Cr}_x|\}} + \max{\{2 r_y, |\mathrm{Cr}_y|\}} \le 1,
\end{equation*}
$$

一个必要条件是

$$
\begin{equation*}
    \tag{17}
    |z_0| + 2 r_x + 2 r_y \le 1.
\end{equation*}
$$

故在

$$
\begin{equation*}
    \tag{18}
    2 r_x \ge |\mathrm{Cr}_x|, \, 2 r_y \ge |\mathrm{Cr}_y|
\end{equation*}
$$

的大前提下, 命题 (11) 的充分条件 (16) 即是必要条件 (17), 从而是充要的. 对于 (18) 以外的情形, (16) 是充分但不必要的. 简单起见, 我们用命题 (16) 的、从而也是命题 (11) 的一个充分**不必要**条件,

$$
\begin{equation*}
    \tag{19}
    |\text{Cr}_y| \le |\text{Cr}_x| \le 2 r_x = 2 r_y \le (1 - K \tau) / 2,
\end{equation*}
$$

来决定网格参数 $h_x, h_y, \tau$. 特别地, 当 $v = 0, \, D_x = D_y = D$, 条件 (19) 等价于

$$
\begin{equation*}
    h_x = h_y =: h \le {2 D / |u|}, \quad \tau \le \left(4 D / h^2 + \max{\{ K_{\mathrm{a}}, K_{\mathrm{r}} \}} \right)^{-1}.
\end{equation*}
$$

&ensp; &emsp; 一般地, 要写出 (11) 的充要条件, 或许可先用 Lagrange multipliers ([2019](https://abel.math.harvard.edu/~knill/teaching/math21a2019/18-lagrange.pdf), [2022](https://abel.math.harvard.edu/~knill/teaching/math21a2022/handouts/lecture29.pdf)) 方法写出 $\max_{|\xi h| \le \pi}|G|$ 的表达式, 再由 (11) 获得 $G$ 的其余参数应满足的约束. 计算似乎相当繁琐, 所以留给感兴趣的读者. 我用 Mathematica ([MMA]((https://reference.wolfram.com/language/ref/ForAll.html))) 软件算得 (11) 的一个充要条件是

$$
\begin{equation*}
    \tag{20}
    \text{MMA 似乎算不出来? 是我给的命令不好吗?}
\end{equation*}
$$

```wolfram{.line-numbers}
G = 1 - 2 rx - 2 ry - k tau + 2 rx cx - I crx sx + 2 ry cy - I crx sy
ForAll[{cx, sx, cy, sy}, cx^2 + sx^2 == 1 && cy^2 + sy^2 == 1, Abs[G] <= 1]
```

若作简化

$$
\begin{equation*}
    \tag{21}
    \text{Cr}_y = 0, \; r_x = r_y = r, \; K\tau = 0,
\end{equation*}
$$

也就是将命题 (21) 添加到大前提, 则用 MMA 显式化命题 (11)

```wolfram{.line-numbers}
Resolve[ForAll[{cx, sx, cy}, cx^2 + sx^2 == 1 && -1 <= cy <= 1, (crx*sx)^2 + (2*cx*r + 2*cy*r - 2*r - 2*r + 1)^2 <= 1], Reals]
FullSimplify[%, r >= 0]
```

将得到形式很复杂的复合命题

$$
\begin{equation*}
    \tag{22}
    \text{crx}^2+4 (r-1) r\leq 0\land 4 r\leq 1\land \left(\text{crx}=2 r\lor r=0\lor \text{crx}+2 r=0\lor \left(\left(\left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right)\leq 0\lor \left(r^2 (2 r-1) \left(4 r^2-\text{crx}^2\right)\leq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\right)\lor 6 r<1\lor 2 r^2\geq r\lor r^2\leq 0\right)\land \left(\text{crx}\neq 0\lor \text{crx}^2 r (4 r-1)\geq 0\lor \text{crx}^2 r (4 r-1) \left(\text{crx}^4-4 \text{crx}^2 r (r+1)+4 r^2\right)\geq 0\lor \left(r^2 \left(\text{crx}^{10}-4 \text{crx}^8 (r (9 r-4)+1)+8 \text{crx}^6 r^2 \left(62 r^2-48 r+11\right)+16 \text{crx}^4 r^2 (4 r (r (3 r (3 r-4)+10)-3)+1)-16 \text{crx}^2 r^4 (8 r (5 r-2)+3)+64 r^6\right)\geq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\right)\lor \left(r^2 \left(\text{crx}^{10}+4 \text{crx}^8 r (7 r-4)-8 \text{crx}^6 r^2 (2 r (r+8)-7)+16 \text{crx}^4 r^2 (4 r (r (5 r (5 r-4)+11)-3)+1)-16 \text{crx}^2 r^4 (8 r (5 r-2)+3)+64 r^6\right)\leq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^4+4 \text{crx}^2 r (3 r-2)+8 r^2 \left(8 r^2-4 r+1\right)\neq 0\right)\land \left(\text{crx}^{10}+4 \text{crx}^8 r (7 r-4)+16 \text{crx}^4 r^2 (4 r (r (5 r (5 r-4)+11)-3)+1)+64 r^6\neq 8 \text{crx}^6 r^2 (2 r (r+8)-7)+16 \text{crx}^2 r^4 (8 r (5 r-2)+3)\lor \text{crx}^2 r^2+4 r^4\leq 0\lor \left(\text{crx}^2+4 r^2\neq 0\land \text{crx}^4+\text{crx}^2 (4 (2-5 r) r-2)+4 r^2\neq 0\right)\right)\right)\lor \left(\text{crx}^4+4 \text{crx}^2 r (7 r-3)+4 r^2\neq 0\land \text{crx}^2 r (4 r-1) \left(\text{crx}^4+4 \text{crx}^2 r (7 r-3)+4 r^2\right)\leq 0\right)\lor \left(\text{crx}^2 r (4 r-1) \left(\text{crx}^4+\text{crx}^2 (4 (2-5 r) r-2)+4 r^2\right)\leq 0\land \left(r^2\leq 0\lor \text{crx}^4+\text{crx}^2 (4 (2-5 r) r-2)+4 r^2\neq 0\right)\right)\right)\land \left(\text{crx}=0\lor \left(\text{crx}^2=2 r\land r>0\right)\lor \left(\left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\land \text{crx}^8 r^2+16 \text{crx}^2 r^6\leq 4 \text{crx}^4 \left(\text{crx}^2+1\right) r^4\right)\lor r=0\lor \left(\text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\neq 0\land \text{crx}^2 r^2 \left(\text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\right)\geq 0\land \text{crx}^2 r (6 r-1)\leq 0\right)\lor \text{crx}^2 r (2 r-1)\geq 0\lor \text{crx}^2 \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right)\geq 0\right)\land \left(2 r^2\geq r\lor \left(\text{crx}^2-2 r\right)^2\leq 0\lor r (2 r-1) \left(\text{crx}^2-2 r\right)>0\right)\land \left(\left(r>0\land \left(r (2 r-1) \left(\text{crx}^2-2 r\right)<0\lor \left(\text{crx}^2 (10 r-3)+2 r (2 r (8 r-3)+1)>0\land 6 r<1\right)\right)\right)\lor \text{crx}^2+4 r^2=0\lor r^2 (2 r-1) \left(16 r^4-\text{crx}^4\right)\leq 0\lor \text{crx}^8 r^2+64 r^8\leq 4 \text{crx}^4 \left(4 r^6+r^4\right)\right)\right)\right)\land \left(\text{crx}=2 r\lor r=0\lor \text{crx}+2 r=0\lor \left(\left(\text{crx}\neq 0\lor \text{crx}^2 r (4 r-1)\geq 0\lor \text{crx}^2 r (4 r-1) \left(\text{crx}^4-4 \text{crx}^2 r (r+1)+4 r^2\right)\geq 0\lor \left(r^2 \left(\text{crx}^{10}-4 \text{crx}^8 (r (9 r-4)+1)+8 \text{crx}^6 r^2 \left(62 r^2-48 r+11\right)+16 \text{crx}^4 r^2 (4 r (r (3 r (3 r-4)+10)-3)+1)-16 \text{crx}^2 r^4 (8 r (5 r-2)+3)+64 r^6\right)\geq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\right)\lor \left(r^2 \left(\text{crx}^{10}+4 \text{crx}^8 r (7 r-4)-8 \text{crx}^6 r^2 (2 r (r+8)-7)+16 \text{crx}^4 r^2 (4 r (r (5 r (5 r-4)+11)-3)+1)-16 \text{crx}^2 r^4 (8 r (5 r-2)+3)+64 r^6\right)\leq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^4+4 \text{crx}^2 r (3 r-2)+8 r^2 \left(8 r^2-4 r+1\right)\neq 0\right)\land \left(\text{crx}^{10}+4 \text{crx}^8 r (7 r-4)+16 \text{crx}^4 r^2 (4 r (r (5 r (5 r-4)+11)-3)+1)+64 r^6\neq 8 \text{crx}^6 r^2 (2 r (r+8)-7)+16 \text{crx}^2 r^4 (8 r (5 r-2)+3)\lor \text{crx}^2 r^2+4 r^4\leq 0\lor \left(\text{crx}^2+4 r^2\neq 0\land \text{crx}^4+\text{crx}^2 (4 (2-5 r) r-2)+4 r^2\neq 0\right)\right)\right)\lor \left(\text{crx}^4+4 \text{crx}^2 r (7 r-3)+4 r^2\neq 0\land \text{crx}^2 r (4 r-1) \left(\text{crx}^4+4 \text{crx}^2 r (7 r-3)+4 r^2\right)\leq 0\right)\lor \text{crx}^2 r (4 r-1) \left(\text{crx}^4+\text{crx}^2 (4 (2-5 r) r-2)+4 r^2\right)\geq 0\right)\land \left(\text{crx}=0\lor \left(\text{crx}^2=2 r\land r>0\right)\lor \left(\left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\land \text{crx}^8 r^2+16 \text{crx}^2 r^6\leq 4 \text{crx}^4 \left(\text{crx}^2+1\right) r^4\right)\lor \left(\text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\neq 0\land \text{crx}^2 r^2 \left(\text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\right)\geq 0\land \text{crx}^2 r (6 r-1)\leq 0\right)\lor \text{crx}^2 r (2 r-1)\geq 0\lor \text{crx}^2 \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right)\geq 0\lor \text{crx}^2 r^2>0\right)\land \left(6 r^2\leq r\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\leq 0\lor r \left(\text{crx}^2 (10 r-3)+2 r (2 r (8 r-3)+1)\right)>0\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+8 r^2 (6 r (3 r-1)+1)<0\right)\land \left(\left(r>0\land \left(r (2 r-1) \left(\text{crx}^2-2 r\right)<0\lor \left(\text{crx}^2 (10 r-3)+2 r (2 r (8 r-3)+1)>0\land 6 r<1\right)\right)\right)\lor \text{crx}^2+4 r^2=0\lor r^2 (2 r-1) \left(16 r^4-\text{crx}^4\right)\leq 0\lor \text{crx}^2 r^2+4 r^4>0\lor \text{crx}^8 r^2+64 r^8\leq 4 \text{crx}^4 \left(4 r^6+r^4\right)\right)\right)\right)\land \left(r=0\lor \text{crx}^2\geq 4 r^2\lor \left(\left(\text{crx}=0\lor \text{crx}^2 \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right)\geq 0\lor \left(\left(\text{crx}^2 r (6 r-1)\leq 0\lor \left(r>0\land \text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\leq 0\right)\right)\land \left(\text{crx}^2 r (6 r-1)<0\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\neq 0\right)\right)\lor \text{crx}^2 r (2 r-1)\geq 0\lor \left(\text{crx}^8 r^2+16 \text{crx}^2 r^6\leq 4 \text{crx}^4 \left(\text{crx}^2+1\right) r^4\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\land \left(\text{crx}^2 \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right) \left(\text{crx}^2+4 r^2\right)^2\geq 0\lor \left(\text{crx}\neq 2 r\land \text{crx}^2\neq 2 r\land \text{crx}+2 r\neq 0\land \text{crx}^2+2 r\neq 0\right)\right)\right)\right)\land \left(6 r=1\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+4 r^2\geq 0\lor r \left(\text{crx}^2 (10 r-3)+2 r (2 r (8 r-3)+1)\right)>0\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+8 r^2 (6 r (3 r-1)+1)<0\right)\land \left(\text{crx}=2 r\lor \left(0<r<\frac{1}{6}\land \text{crx}^4+4 \text{crx}^2 r (8 r-3)+8 r^2 (6 r (3 r-1)+1)>0\right)\lor \text{crx}+2 r=0\lor \text{crx}^2+4 r^2=0\lor \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^4-16 r^4\right)\geq 0\right)\land \left(\text{crx}=2 r\lor 0<r<\frac{1}{6}\lor \text{crx}^4+4 \text{crx}^2 r (8 r-3)+8 r^2 (6 r (3 r-1)+1)<0\lor \text{crx}+2 r=0\lor \text{crx}^2+4 r^2=0\lor \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^4-16 r^4\right)\geq 0\right)\land \left(\text{crx}=2 r\lor r \left(\text{crx}^2 (10 r-3)+2 r (2 r (8 r-3)+1)\right)<0\lor \text{crx}+2 r=0\lor \text{crx}^2+4 r^2=0\lor 0<r<\frac{1}{6}\lor r^2 (2 r-1) \left(16 r^4-\text{crx}^4\right)<0\lor 4 \text{crx}^4 \left(4 r^6+r^4\right)\leq \text{crx}^8 r^2+64 r^8\right)\land \left(\text{crx}=2 r\lor r (2 r-1) \left(\text{crx}^2-2 r\right)<0\lor \text{crx}+2 r=0\lor \text{crx}^2+4 r^2=0\lor \left(\text{crx}^2 (3-10 r)<2 r (2 r (8 r-3)+1)\land 6 r<1\right)\lor r^2 (2 r-1) \left(16 r^4-\text{crx}^4\right)<0\lor 4 \text{crx}^4 \left(4 r^6+r^4\right)\leq \text{crx}^8 r^2+64 r^8\right)\land \left(6 r<1\lor 2 r^2\geq r\lor \left(\text{crx}^2-1\right) r^2 \left(\text{crx}^2-4 r^2\right)\geq 0\lor \left(\left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^4+4 \text{crx}^2 r (3 r-2)+8 r^2 \left(8 r^2-4 r+1\right)\neq 0\right)\land \left((\text{crx}\neq 2 r\land \text{crx}+2 r\neq 0)\lor \text{crx}^2 r^2+4 r^4\leq 0\right)\right)\lor \left(r^2 (2 r-1) \left(4 r^2-\text{crx}^2\right)\leq 0\land \left(\text{crx}^2+4 r^2\neq 0\lor \text{crx}^2 (5 r-2)+2 r \left(6 r^2-4 r+1\right)\neq 0\right)\right)\right)\right)\right)
\end{equation*}
$$

若在命题 (21) 成立的条件下, 进一步让 (18) 加入大前提, 则期待 (16) 能保证 (11). 用 MMA 验证,

```wolfram{.line-numbers}
Resolve[ForAll[{cx, sx, cy}, cx^2 + sx^2 == 1 && -1 <= cy <= 1, (crx*sx)^2 + (2*cx*r + 2*cy*r - 2*r - 2*r + 1)^2 <= 1], Reals]
FullSimplify[%, 0 < Abs[crx] <= 2 r && Abs[-2*r - 2*r + 1] + 2 r + 2 r <= 1]
```

得

$$
\tag{23}
\text{Cr}_x^2 + 36 r^2 \leq 12 r.
$$

即 MMA 认为, 在 (16)(18)(21) 成立的大前提下, 命题 (11) 等价于命题 (23), 而后者并非恒成立; 这就否定了 (16) 对 (11) 的充分性, 表明 <b style="color: blue;">MMA 不支持本节的主要结果 (16). 那么, 问题出在哪? 请求读者帮助!</b> 文末附有[联系方式](#contact-information).

### Problem 1.3

&ensp; &emsp; 本节讨论格式 (3) 的精度. 格式 (3a) 的 modified equation

$$
\begin{equation*}
    \tag{24}
    \begin{aligned}
        & \qquad B_t + \boldsymbol{u} \cdot \nabla B - \nabla \cdot (D \nabla B) + K_{\text{r}} B \\
        & = -\frac{\tau}{2} B_{t^2} - \frac{h_x^2}{6} u B_{x^3} - \frac{h_y^2}{6} v B_{y^3} + \frac{h_x^2}{12} D_x B_{x^4} + \frac{h_y^2}{12} D_y B_{y^4} + \mathcal{O}(\tau^2 + h_x^2 + h_y^2).
    \end{aligned}
\end{equation*}
$$

试图把 (24) 右端的时间误差项主部系数 $B_{t^2}$ 用空间导数表示, 所以让 (24) 分别被 $\partial_t$ 和 $\nabla$ 算子作用, 得

$$
\begin{equation*}
    \tag{25}
    B_{t^2} + \boldsymbol{u} \cdot \nabla B_t - \nabla \cdot (D \nabla B_t) + K_{\mathrm{r}} B_t = \mathcal{O}(\tau + h_x^2 + h_y^2),
\end{equation*}
$$

$$
\begin{equation*}
    \tag{26}
    \nabla B_{t} + (\boldsymbol{u} \cdot \nabla) \nabla B - \nabla \left[ \nabla \cdot (D \nabla B) \right] + K_{\mathrm{r}} \nabla B = \mathcal{O}(\tau + h_x^2 + h_y^2).
\end{equation*}
$$

由 (25)(26) 和 (24) 得

$$
\begin{equation*}
    \tag{27}
    \begin{aligned}
        B_{t^2} & = (\boldsymbol{u} \cdot \nabla) \left[ \left( (\boldsymbol{u} \cdot \nabla) - K_{\text{r}} \right) B - \nabla \cdot (D \nabla B) \right] \\
        & - K_{\mathrm{r}} \left[ \boldsymbol{u} \cdot \nabla B - \nabla \cdot (D \nabla B) + K_{\text{r}} B \right] \\
        & - \nabla \cdot \left[ D (\boldsymbol{u} \cdot \nabla) \nabla B - D \nabla \left[ \nabla \cdot (D \nabla B) \right] + K_{\mathrm{r}} D \nabla B \right] \\
        & + \mathcal{O}(\tau + h_x^2 + h_y^2).
    \end{aligned}
\end{equation*}
$$

Eq. (27) 有些复杂, 感兴趣的读者可尝试将其代入 (24), 整理出只含空间导数的截断误差项主部, 看有无可能通过适当选取网格参数 $\tau, h_x, h_y$ 来消去截断误差中的 $\mathcal{O}(\tau)$ 项. 我对此不抱多少期望, 所以我说格式 (3a) 是 $\mathcal{O}(\tau + h_x^2 + h_y^2)$ 精度的, so does (3b), and their combination (3).

### Problem 1.4

&ensp; &emsp; 将格式 (3) 写成迭代形式,

$$
\begin{equation*}
    \tag{28}
    \begin{aligned}
        \begin{bmatrix*}
            B \\
            O
        \end{bmatrix*}^{n+1}_{i,j} & = (1 - 2 r_x - 2 r_y) \begin{bmatrix*}
            B\\
            O
        \end{bmatrix*}^n_{i,j} - \tau \begin{bmatrix*}
            K_{\mathrm{r}} & \\
            K_{\mathrm{r}} & K_{\mathrm{a}}
        \end{bmatrix*} \begin{bmatrix*}
            B \\
            O
        \end{bmatrix*}^{n}_{i,j} + \tau \begin{bmatrix}
            q_{i,j} \\
            K_{\mathrm{a}} O_{\mathrm{sat}}
        \end{bmatrix} \\
        & + \sum_{k \in \{ -1, 1 \}} \left[
            \left( r_x - \frac{\mathrm{Cr}_x}{2} \mathop{\mathrm{sgn}}{k} \right) \begin{bmatrix}
                B \\
                O
            \end{bmatrix}^n_{i+k, j} + \left( r_y - \frac{\mathrm{Cr}_y}{2} \mathop{\mathrm{sgn}}{k} \right) \begin{bmatrix}
                B \\
                O
            \end{bmatrix}^n_{i, j+k}
        \right],
    \end{aligned}
\end{equation*}
$$

其中 BOD 的源项

$$
\begin{equation*}
    \tag{29}
    q_{i,j} = \sum_{(x_0, y_0)} \frac{Q_{\mathrm{B}}(x_0, y_0)}{h_x h_y} \mathrm{I}_{i,j} (x_0, y_0),
\end{equation*}
$$

where $(x_0, y_0)$ 是 BOD 的一个源点, $Q_{\mathrm{B}}(x_0, y_0)$ 是这源的强度 (加载率), $\mathrm{I}_{i, j}(x_0, y_0)$ 是依 (4) 定义的示性函数. 在网格边界, 按 (7)(8') 设置 halo points. 用 [Matlab](https://www.mathworks.com/products/matlab.html) 实现数值求解二维 Streeter-Phelps 模型 (1)(2) 的 FTCS 方案 (28).

&ensp; &emsp; 在 (28) 中, 令

$$
\begin{equation*}
    \tag{30}
    \begin{bmatrix*}
        B \\
        O
    \end{bmatrix*}^{n+1} = \begin{bmatrix*}
        B \\
        O
    \end{bmatrix*}^{n} =: \begin{bmatrix*}
        B^{*} \\
        O^{*}
    \end{bmatrix*},
\end{equation*}
$$

得数值稳态 (steady-state) 解 $B^*, O^*$ 满足的方程

$$
\begin{equation*}
    \tag{31}
    \begin{aligned}
        \tau \begin{bmatrix}
            q_{i,j} \\
            K_{\mathrm{a}} O_{\mathrm{sat}}
        \end{bmatrix} & = (2 r_x + 2 r_y) \begin{bmatrix*}
            B^{*} \\
            O^{*}
        \end{bmatrix*}_{i,j} + \tau \begin{bmatrix*}
            K_{\mathrm{r}} & \\
            K_{\mathrm{r}} & K_{\mathrm{a}}
        \end{bmatrix*} \begin{bmatrix*}
            B^{*} \\
            O^{*}
        \end{bmatrix*}_{i,j} \\
        & - \sum_{k \in \{ -1, 1 \}} \left[
            \left( r_x - \frac{\mathrm{Cr}_x}{2} \mathop{\mathrm{sgn}}{k} \right) \begin{bmatrix}
                B^{*} \\
                O^{*}
            \end{bmatrix}_{i+k, j} + \left( r_y - \frac{\mathrm{Cr}_y}{2} \mathop{\mathrm{sgn}}{k} \right) \begin{bmatrix}
                B^{*} \\
                O^{*}
            \end{bmatrix}_{i, j+k}
        \right].
    \end{aligned}
\end{equation*}
$$

要求解 (31), 可先求解 (31) 中关于 $\{ B^*_{i,j} \}$ 的常系数线性方程. 引入 halo points (7)(8') 后, 这方程的系数矩阵是确定的方阵, 且含大量零元素 (是 sparse matrix (稀疏矩阵)), 适合利用 Matlab 的 [`sparse()`](https://www.mathworks.com/help/matlab/ref/sparse.html) 函数来存储. 求出 $B^*$ 后, 从 (31) 和 (7)(8') 可建立关于 $\{ O^*_{i,j} \}$ 的常系数线性方程, 从而确定 FTCS 格式 (28) 的数值稳态解 (30). 我们期待数值解有渐近行为

$$
\begin{equation*}
    \tag{32}
    \begin{bmatrix*}
        B \\
        O
    \end{bmatrix*}^n \to \begin{bmatrix*}
        B^{*} \\
        O^{*}
    \end{bmatrix*}, \quad (n \to \infty).
\end{equation*}
$$

&ensp; &emsp; 为确保 FTCS 格式 (28) 稳定且收敛, 应适当选取网格参数 $h_x, h_y$ 和 $\tau$. 根据第 [1.3](#problem-13) 节的分析, 当 $v = 0, \, D_x = D_y = D$, 格式稳定的一个充分不必要条件 (19) 等价于

$$
\begin{equation*}
    \tag{33}
    h_x = h_y =: h \le {2 D / |u|}, \quad \tau \le \left(4 D / h^2 + \max{\{ K_{\mathrm{a}}, K_{\mathrm{r}} \}} \right)^{-1}.
\end{equation*}
$$

由 (33), 我们选择 $h = 5 \times 10^{-2} \, \text{m}, \, \tau = 1.24 \times 10^{-3} \, \text{s}$.

&ensp; &emsp; 我们测试了以下四种排放方案:

1. **单点排放** (single-point discharge). 所有 BOD 都排放在 $(x, y) = (10, 15) \, \text{m}$ 处.
2. **两点排放** (two-points discharge). BOD 被平均地排放在 $(x, y) = (10, 15), (70, 15)$ 两个点上.
3. **线分布排放** (line-distributed discharge). BOD 排放源均匀分布在线段 $x \in [10, 70], \, y = 15$ 上.
4. **面分布排放** (face-distributed discharge). BOD 排放源均匀分布在区域 $(x, y) \in [10, 70] \times [10, 20]$ 上.

测试结果示于下面 8 个图. 每个图都标注了相应时刻的 BOD 或 DO 的最值. 发现以下特点:

1. 对于点排放 (单点, 两点) 情形, BOD 最高值逐渐升高至稳态值, 且位于最下游的排放源点. DO 最低值逐渐下降至稳态值, 且位置随时间推移有向下游移动的趋势, 最终稳定在排放源点的下游.
2. 对于均匀分布排放 (线, 面) 情形, BOD 最高值逐渐升高至稳态值, 且位置随时间推移有从源区中心向下游移动的趋势, 最终稳定在源区下游边界的内侧. DO 最低值逐渐下降至稳态值. 随时间推移, DO 最低值点有从源区中心向下游移动的趋势, 且相对于 BOD 最高值点, 有从滞后 (即位于相对上游) 到超前的移动趋势.
3. 采取分布排放, 可以大幅降低 BOD 峰值, 因为相比点排放可以大幅降低排放源强度. 但未必能使 DO 最低值达到最高. 物理解释: 点排放引起 DO 在小范围内有快速下降的倾向, 分布排放引起 DO 在大范围内有缓慢下降的倾向. 若平流速度充分大或点排放源强度相对足够小, 则对于点排放情形, 水团受点排放源的影响时间较小, 引起的 DO 损失足够被水-气氧气通量来补偿; 而对于分布排放, 水团在较长时间内持续地被 BOD 反应掉, 足以产生比点排放情形更低的 DO 最低值.

![fig_1_BOD_single-point](fig/fig_1_BOD_single-point.jpg "single-point discharge: BOD")

![fig_1_DO_single-point](fig/fig_1_DO_single-point.jpg "single-point discharge: DO")

![fig_1_BOD_two-points](fig/fig_1_BOD_two-points.jpg "two-points discharge: BOD")

![fig_1_DO_two-points](fig/fig_1_DO_two-points.jpg "two-points discharge: DO")

![fig_1_BOD_line-distributed](fig/fig_1_BOD_line-distributed.jpg "line-distributed discharge: BOD")

![fig_1_DO_line-distributed](fig/fig_1_DO_line-distributed.jpg "line-distributed discharge: DO")

![fig_1_BOD_face-distributed](fig/fig_1_BOD_face-distributed.jpg "face-distributed discharge: BOD")

![fig_1_DO_face-distributed](fig/fig_1_DO_face-distributed.jpg "face-distributed discharge: DO")

### Problem 1.5

&ensp; &emsp; 第 [4](#problem-14) 节 8 个图的子图 (d) 展示了四种排放情形下 BOD 和 DO 的稳态分布. 可见, 最值点都出现在河道中心线 ($y = 15$) 上. 下图展示了河道中心线断面上 BOD 和 DO 的稳态分布, 结果支持第 [4](#problem-14) 节的分析.

![fig_2_steady_state_at_central_line](fig/fig_2.svg "steady state of BOD and DO at river central line.")

### Problem 1.6

&ensp; &emsp; 第 [4](#problem-14) 和 [5](#problem-15) 节的结果表明, 单点排放 BOD 方案导致的 DO 稳态值不满足 Environmental Protection Agency (EPA) 的标准 ($\ge 6 \, \mathrm{mg} \, \mathrm{L}^{-1}$). 若改用两点排放方案, 或分布 (线, 面) 排放方案, 则可达到 EPA 标准.

### Problem 1.7

&ensp; &emsp; 本文得到的计算结果, 与真实情形相比, 存在模型误差、观测误差、截断误差 (方法误差) 和舍入误差这四种误差. 以下方法或有助于改善之.

1. 完善数学模型. 考虑在经典版本的 S-P 模型 (1)(2) 中, 设法引入对动力过程、物理过程、生物地球化学过程的更精细描述. 这样可以减少模型误差.
2. 改善初始场、边界场对观测的描述质量, 选用更接近真实的模型参数. 这样可以减少观测误差.
3. 改用有更高阶精度离散格式. 本文对模型采取 FTCS 离散格式 (28), 这是一个时间一阶、空间二阶的显式方案 (第 [3](#problem-13) 节). 可以尝试时间上的梯形法, Runge-Kutta 法等有更高阶精度的方案. 这样可以减少截断误差.
4. 加密网格, 缩小时间步长. 这样可以减少截断误差.
5. 用支持更高精度数据表示的计算机 (但现在基本都是 64 位计算机, 似乎没听说要推广有更宽数据总线的 CPU ?), 和编程语言. 这样可以减少舍入误差?

### Problem 1.8

&ensp; &emsp; 若要提高数值计算速度, 可考虑以下方法:

1. 加长时间步, 同时加粗网格以保证格式稳定.
2. 换用稳定性更好的离散格式, 从而用更大的时间步.
3. 改善算法的实现. 例如求稳态解时, 考虑线性方程的稀疏性, 优化存储和计算 (已采用). 又如, 换用计算速度更快的编程语言, 如 Fortran. 又如, 在编程中通过向量化, 以及其他巧妙的程序设计技术 (包括但不限于, 以空间换时间, 减少冗余操作), 来优化程序性能. 可惜作者的算法和程序设计水平有限.
4. 用更高性能的计算平台, 更优的并行计算方案.
5. 不知是否存在比显式格式 (28) 时间复杂度更低的格式? 若有, 采用之.

## Acknowledgement

&ensp; &emsp; I am grateful to my fellow PhD student, Dan Cao (曹丹), whose insightful discussions ... (关于命题 (11) 的讨论促使首次正式使用 MMA, 关于 solid wall 边界如何引入的讨论使增进物理理解)

## Contact Information

- **Author:** Guorui Wei (危国锐)
- **E-mail:** [313017602@qq.com](mailto:313017602@qq.com)
- **Website:** <https://github.com/grwei>

## Appendix

### Solid wall B.C. 对 Fick 扩散的约束 (ChatGPT mini 4o)

> In physics, how does a solid wall boundary constrain Fick's diffusion? Please explain it both physically and mathematically. (13 NOV 2024)

Fick's law of diffusion describes the process by which particles spread from regions of high concentration to low concentration. In the context of a solid wall boundary, the diffusion is constrained because the particles cannot pass through the wall, effectively limiting their movement in that direction.

Let's break it down both physically and mathematically:

#### Physical Interpretation

Imagine a particle diffusing in a medium, such as gas or liquid. According to Fick’s law, the particle will naturally move from regions of high concentration to low concentration, attempting to reach equilibrium.

When the particle encounters a boundary (like a solid wall), the wall acts as an impermeable surface (assuming no permeability or flux through the wall). The particles cannot continue their motion beyond the wall. This constraint leads to what is called the "no-flux" condition at the boundary, meaning that the net flux of particles perpendicular to the wall must be zero. Physically, this reflects the fact that no particles can escape across the wall, but they may scatter or reflect back, depending on the wall's nature (for example, if it’s a reflective wall).

#### Mathematical Interpretation

Mathematically, the constraint can be expressed in terms of the **boundary condition** for the diffusion equation.

Fick’s first law states that the diffusion flux \( J \) is related to the concentration gradient \( \nabla C \) by:

\[
J = -D \nabla C
\]

where:

- \( J \) is the diffusion flux (amount of substance diffusing per unit area per unit time),
- \( D \) is the diffusion coefficient (which depends on the medium and particle properties),
- \( \nabla C \) is the gradient of concentration \( C \).

For a one-dimensional problem, where diffusion occurs along the \( x \)-axis, the flux becomes:

\[
J(x) = -D \frac{dC}{dx}
\]

Now, at a solid wall boundary (say, located at \( x = 0 \)), the boundary condition reflects the fact that no particles can cross the wall. This is mathematically expressed as the **no-flux condition**:

\[
\left. \frac{dC}{dx} \right|_{x=0} = 0
\]

This condition means that the derivative of the concentration with respect to \( x \) (the spatial direction normal to the wall) is zero at the wall. In other words, the concentration gradient at the wall is zero, implying that no flux occurs through the boundary. This reflects the physical idea that the particles cannot cross the wall, and thus there is no change in concentration at the boundary in the direction perpendicular to the wall.

#### Summary

- **Physically**, the solid wall creates a boundary that particles cannot cross, causing them to either accumulate at the wall or be reflected back into the medium.
- **Mathematically**, this is expressed by the no-flux boundary condition \( \frac{dC}{dx} = 0 \) at the location of the wall, which ensures that there is no net particle flow through the wall.

This no-flux condition is critical in solving diffusion problems that involve solid boundaries, as it shapes the concentration profile near the wall.

### Matlab code for problem 1

```Text
.
├── @DoSag2D
│   ├── DoSag2D.m
│   ├── bSrc.m
│   ├── getSolution.m
│   ├── prepareParams.m
│   ├── reset.m
│   ├── solve2D.m
│   └── solve2DSteady.m
├── fig
├── hw5_1.m
├── inc
├── lib
└── log

6 directories, 8 files
```

#### hw5_1.m

```matlab {.line-numbers}
%% hw5_1.m
% Description: Solving the Streeter-Phelps model in a 2-D river
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Dec. 1, 2024
% Last modified: Dec. 4, 2024
%

clc; clear; close all

%%% 程序参数

FLAG_EXPORT_GRAPH = true;
RESULTS_PATH = ".\lib\hw5_1_20241203_080420.mat";
LOG_FILE_PATH = ".\log\log_" + string(datetime("now"), "yyyyMMdd_HHmmss", "en_US") + ".log";

%%% 环境准备

if ~isfolder(".\log\")
    mkdir(".\log\")
end
LOG_FILE_ID = fopen(LOG_FILE_PATH,'a');
fprintf(LOG_FILE_ID, "%s\n\tStart.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"));

%% 1. 求解

if ~isfile(RESULTS_PATH)
    t_start = tic;
    warning("%s\n\tresults file:\n\t%s\n\tdoes not exist, 正在重新求解.\n\t求解完成后, 请相应更改 script 中的 RESULTS_PATH.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH)
    fprintf(LOG_FILE_ID, "%s\n\tresults file:\n\t%s\n\tdoes not exist, 正在重新求解.\n\t求解完成后, 请相应更改 script 中的 RESULTS_PATH.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH);
    RESULTS_PATH = hw5_1_4();
    warning("%s\n\t已重新求解, 耗时 %.1f min. 请将 script 中的 RESULTS_PATH 更改为:\n\t%s\n", toc(t_start) / 60, string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH)
    fprintf(LOG_FILE_ID, "%s\n\t已重新求解, 耗时 %.1f min. 请将 script 中的 RESULTS_PATH 更改为:\n\t%s\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), toc(t_start) / 60, RESULTS_PATH);
end

%% 2. 绘图

load(RESULTS_PATH, "solvers");

time_inds = [1,2,5,6];
test_names = [
    "single-point";
    "two-points";
    "line-distributed";
    "face-distributed"
];
var_names = [
    "BOD";
    "DO";
];

%%% 2.2 各试验, 河道中心线上的 BOD, DO 稳态分布

fig_name = sprintf("fig_2");
varSteadyDistYline(solvers, 15, fig_name, FLAG_EXPORT_GRAPH);

%%% 2.1 各试验的 BOD, DO 时空演变图

for var_ind = 1:length(var_names)
    for test_ind = 1:length(solvers)
        title_str = sprintf("%s (%s discharge)", var_names(var_ind), test_names(test_ind));
        fig_name = sprintf("fig_1_%s_%s", var_names(var_ind), test_names(test_ind));
        varSpatialTemporalDist(solvers{test_ind}, var_ind, time_inds, title_str, fig_name, FLAG_EXPORT_GRAPH);
    end
end

%% 3. 程序结束

fprintf(LOG_FILE_ID, "%s\n\tFinished.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"));
fclose(LOG_FILE_ID);

% clc; clear; close all;

%% 

function [] = varSteadyDistYline(solvers, y_query, fig_name, flag_export_graph, fig_dir)
    arguments
        solvers 
        y_query 
        fig_name 
        flag_export_graph = false;
        fig_dir = ".\fig\"
    end
    test_names = [
        "single-point discharge";
        "two-points discharge";
        "line-distributed discharge";
        "face-distributed discharge"
    ];

    t_fig = figure(Name=fig_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 18, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, length(solvers), 1, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
    ylabel(t_TCL, "biochemical Oxygen demand (BOD) (mg/L)", Interpreter="latex", FontSize=10, Color="#0072BD");

    Y_IND = NaN;
    for ind_ = 1:length(solvers)
        if isnan(Y_IND)
            [~, Y_IND] = min(solvers{ind_}.y_grid - y_query, [], "ComparisonMethod", "abs");
        end

        t_axes = nexttile(t_TCL, ind_);
        plot(t_axes, solvers{ind_}.x_grid, solvers{ind_}.f_list{end}{1}(:, Y_IND), LineWidth=.8);
        hold(t_axes, "on")
        set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", ...
            Tag=sprintf("(%c) %s", char('a' - 1 + ind_), test_names(ind_)));
        xline(t_axes, [10, 70], "--", Color="#A9A9A9")
        yyaxis right
        plot(t_axes, solvers{ind_}.x_grid, solvers{ind_}.f_list{end}{2}(:, Y_IND), LineWidth=.8);
        set(t_axes, YColor = "#D95319", YLim = [5.5, 8], XLimitMethod="tight");
        hold(t_axes, "off")
        yyaxis left
        t_axes.YColor = "#0072BD";
        yyaxis right
        
        if ind_ < length(solvers)
            t_axes.XTickLabel = {};
        end
    end
    title(t_TCL, sprintf("steady state ($y$ = %.3g m)", solvers{ind_}.y_grid(Y_IND)), Interpreter="latex", FontSize=10);

    linkaxes(findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}), 'xy');

    % add right ylabel
    % cb = colorbar(t_axes, Visible="off", TickLabels={});
    % cb.Position(3) = .5*cb.Position(3);
    % cb.Layout.Tile = "east";
    annotation(t_fig, "textbox", String="dissolved Oxygen (DO) (mg/L)", Color="#D95319", FitBoxToText=false, ...
        Position=[sum(t_TCL.InnerPosition([1,3])), sum(t_TCL.InnerPosition([2,4])), .89*(sum(t_TCL.InnerPosition([2,4])) - t_TCL.InnerPosition(2)), 1.2*(sum(t_TCL.OuterPosition([1,3])) - sum(t_TCL.InnerPosition([1,3])))], ...
        FontSize=10, Rotation=-90, Interpreter="latex", LineStyle="none", HorizontalAlignment="center", VerticalAlignment="top");

    % 子图编号
    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        ind_ = t_axes.Layout.Tile;
        var_data{2} = solvers{ind_}.f_list{end}{2}(:, Y_IND); % DO
        var_data{1} = solvers{ind_}.f_list{end}{1}(:, Y_IND); % BOD
        [BOD_max, x_ind_BOD] = max(var_data{1});
        [DO_min, x_ind_DO] = min(var_data{2});
        str_{3} = sprintf("min(DO) = %.3g ($x$ = %.3g)", DO_min, solvers{ind_}.x_grid(x_ind_DO));
        str_{2} = sprintf("max(BOD) = %.3g ($x$ = %.3g)", BOD_max, solvers{ind_}.x_grid(x_ind_BOD));
        str_{1} = "\bf" + t_axes.Tag;
        t_txt_box = annotation(t_fig, "textbox", String=str_{1}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*20, 10*2], 10*20, 10*2];
        t_txt_box.Units = UNIT_ORIGINAL;

        t_txt_box = annotation(t_fig, "textbox", String=str_(2:end), Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]).*[1, 0], .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*20, 0], 10*20, 10*3];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_export_graph
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        print(t_fig, fig_dir + t_fig.Name + ".svg", "-vector", "-dsvg")
        % !很不合理的现象: 下面这行语句, 导出的图片, 与从显示器看到的不同.
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".jpg", Resolution=800, BackgroundColor="none");
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".pdf", ContentType="vector", BackgroundColor="none");
        % close(t_fig)
    end
end

%%

function [] = varSpatialTemporalDist(solver, var_ind, time_inds, title_str, fig_name, flag_export_graph, fig_dir)
    arguments
        solver      DoSag2D
        var_ind     {isscalar} = 1;
        time_inds   {isvector} = [1,2,5,6];
        title_str   string = "BOD (single-point discharge)";
        fig_name    string = "fig_1_BOD_single_point";
        flag_export_graph = false;
        fig_dir     string = ".\fig\"
    end

    % var_name = "var_name";
    if var_ind == 1
        var_name = "BOD";
    elseif var_ind == 2
        var_name = "DO";
    else
        error("Unknown variable name.")
    end

    % 下载并应用 colormap

    if ~isfolder(".\inc\")
        mkdir(".\inc\")
    end
    if ~isfile(".\inc\cmocean_oxy.rgb")
        system('curl "https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/cmocean_oxy.rgb" >> ".\inc\cmocean_oxy.rgb"');
    end
    if ~isfile(".\inc\cmocean_matter.rgb")
        system('curl "https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/cmocean_matter.rgb" >> ".\inc\cmocean_matter.rgb"');
    end
    if var_name == "DO"
        col_map = uint8(readmatrix(".\inc\cmocean_oxy.rgb", FileType="text", NumHeaderLines=2, ExpectedNumVariables=3));
        cb_str = "dissolved Oxygen (DO) (mg/L)";
        ct_levels = [6, 6.5, 7, 7.5];
        clim_ = [5.5, 8];
    elseif var_name == "BOD"
        col_map = uint8(readmatrix(".\inc\cmocean_matter.rgb", FileType="text", NumHeaderLines=2, ExpectedNumVariables=3));
        cb_str = "biochemical Oxygen demand (BOD) (mg/L)";
        ct_levels = 5;
        clim_ = [];
    else
        error("未定义 col_map");
    end

    t_fig = figure(Name=fig_name);
    colormap(t_fig, col_map)

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 18, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % 决定 clim_
    if isempty(clim_)
        clim_ = [0, 1.1*max(solver.f_list{time_inds(end)}{var_ind}, [], "all")];
    end

    % plot
    t_TCL = tiledlayout(t_fig, length(time_inds), 1, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
    ylabel(t_TCL, "$y$ (m)", Interpreter="latex", FontSize=10);
    title(t_TCL, title_str, Interpreter="latex", FontSize=10);

    for ind_ = 1:length(time_inds)
        var_data = solver.f_list{time_inds(ind_)}{var_ind}.';

        t_axes = nexttile(t_TCL, ind_);
        pcolor(t_axes, solver.x_grid, solver.y_grid, var_data, FaceColor="interp", EdgeColor="none");
        hold(t_axes, "on")
        if var_name == "BOD"
            ct_levels = quantile(var_data, [.25, .50, .75, .90], "all");
        end
        [C, h] = contour(t_axes, solver.x_grid, solver.y_grid, var_data, ct_levels, LineWidth=.6, ShowText = "on", LabelFormat="%.3g", LabelColor="#808080", EdgeColor="#808080");
        clabel(C, h, FontSize=7.5, Interpreter="latex");
        hold(t_axes, "off")
        clim(t_axes, clim_);
        set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", YLimitMethod="tight", ...
            Tag=sprintf("(%c) $t$ = %.2f s", char('a' - 1 + ind_), solver.t_list(time_inds(ind_))));
        axis(t_axes, "equal")
        if ind_ < length(time_inds)
            t_axes.XTickLabel = {};
        end
    end

    cb = colorbar(t_axes, TickLabelInterpreter = "latex", FontSize = 10);
    cb.Layout.Tile = "east";
    set(cb.Label, String=cb_str, Interpreter = "latex", FontSize = 10);

    linkaxes(findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}), 'xy');
    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        ind_ = t_axes.Layout.Tile;
        var_data = solver.f_list{time_inds(ind_)}{var_ind};
        [min_val, min_ind_lin] = min(var_data, [], "all", "linear");
        [min_ind_x, min_ind_y] = ind2sub(size(var_data), min_ind_lin);
        [max_val, max_ind_lin] = max(var_data, [], "all", "linear");
        [max_ind_x, max_ind_y] = ind2sub(size(var_data), max_ind_lin);
        str_{1} = t_axes.Tag;
        str_{2} = sprintf("min. = %.2f at (%.2f, %.2f),\nmax. = %.2f at (%.2f, %.2f)", min_val, solver.x_grid(min_ind_x), solver.y_grid(min_ind_y), max_val, solver.x_grid(max_ind_x), solver.y_grid(max_ind_y));

        t_txt_box = annotation(t_fig, "textbox", String=str_{1}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*30, 10*1.5], 10*30, 10*1.5];
        t_txt_box.Units = UNIT_ORIGINAL;

        t_txt_box = annotation(t_fig, "textbox", String=str_{2}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]).*[1, 0], .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*30, 0], 10*30, 10*3];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_export_graph
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        % print(t_fig, fig_dir + t_fig.Name + ".svg", "-vector", "-dsvg") % 可能导致程序崩溃
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".jpg", Resolution=800, BackgroundColor="none");
        % exportgraphics(t_fig, fig_dir + t_fig.Name + ".pdf", ContentType="vector", BackgroundColor="none"); % 可能导致程序崩溃
        % close(t_fig)
    end
end

%%

function results_path = hw5_1_4()
    %%% 0. 参数设置

    Q_B_SUM = 70;
    LIB_DIR = ".\lib\";

    if ~isfolder(LIB_DIR)
        mkdir(LIB_DIR);
    end

    %%% 1. 求解

    solvers = cell([1, 4]);
    solvers_path = solvers;

    % 1.1 单点排放

    desc_str = "单点排放.\n";
    Q_B_loc = [10, 15; 70, 15]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = [Q_B_SUM; 0];    % BOD 点源强度

    [solvers{1}, solvers_path{1}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.2 两点排放

    desc_str = "两点排放.\n";
    Q_B_loc = [10, 15; 70, 15]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{2}, solvers_path{2}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.3 分布排放

    desc_str = "线分布排放.\n";

    x_src = solvers{1}.x_grid(solvers{1}.x_grid >= 10 & solvers{1}.x_grid <= 70);
    Q_B_loc = [x_src, mean(solvers{1}.range_y)*ones([length(x_src), 1])]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{3}, solvers_path{3}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.4 面分布排放

    desc_str = "面分布排放.\n";
    x_src = solvers{1}.x_grid(solvers{1}.x_grid >= 10 & solvers{1}.x_grid <= 70);
    y_src = solvers{1}.y_grid(solvers{1}.y_grid >= 10 & solvers{1}.y_grid <= 20);
    [X, Y] = ndgrid(x_src, y_src);
    Q_B_loc = [X(:), Y(:)]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{4}, solvers_path{4}] = createSolvers(desc_str, Q_B_loc, Q_B);

    %%% 2. 保存结果

    results_path = LIB_DIR + sprintf("hw5_1_%s.mat", string(datetime("now"), "yyyyMMdd_HHmmss", "en_US"));
    save(results_path, "solvers");
    delete(solvers_path{:});
end

%% 求解

function [solver, solver_path] = createSolvers(desc_str, Q_B_loc, Q_B, t_query, delta_x, delta_y, delta_t, t_start, range_x, range_y, velocity, diffu_coeff, K_r, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1)
    arguments (Input)
        desc_str string = "description of this solver\n"
        Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
        Q_B = [70; 0];    % BOD 点源强度
        t_query = [0, 15.64, 41.576, 84.26, 112.036, +Inf]; % [0, 15.64, 41.576, 84.26, 112.036, +Inf];
        delta_x = .05; % .05 | .1 | .2 | 1
        delta_y = delta_x;
        delta_t = .00124; % .001249 | .00499 | .0199 | .4950
        t_start = 0;
        range_x = [0, 150];
        range_y = [0, 30];
        velocity = [.4, 0];
        diffu_coeff = [.5, .5];
        K_r = .01;
        K_a = .02;
        O_s = 8;
        init_func = {@(x, y) zeros(size(x)); @(x, y) O_s*ones(size(x))};
        bndry_func_type_x0 = "Dirichlet";
        bndry_func_x0 = {@(y, t) zeros(size(y)); @(y, t) O_s*ones(size(y))};
        bndry_func_type_x1 = "Neumann";
        bndry_func_x1 = {@(y, t) zeros(size(y)); @(y, t) zeros(size(y))};
        bndry_func_type_y0 = "Neumann";
        bndry_func_y0 = {@(x, t) zeros(size(x)); @(x, t) zeros(size(x))};
        bndry_func_type_y1 = bndry_func_type_y0;
        bndry_func_y1 = bndry_func_y0;
    end
    arguments (Output)
        solver DoSag2D
        solver_path {isfile}
    end

    %%% 1. 输入合法性检查

    % 1.1 FTCS 格式稳定性

    % Cr = delta_t * velocity ./ [delta_x, delta_y];
    r = diffu_coeff .* delta_t ./ [delta_x, delta_y].^2;
    if any(abs(1 - 2*sum(r) - [K_r, K_a]*delta_t) + 2*sum(r) > 1)
        warning("格式的稳定性存疑. 建议 delta_t < %g", 1 / (4*diffu_coeff(1)/delta_x^2 + max([K_a, K_r])));
    end

    %%% 2. 求解

    solver = DoSag2D(desc_str, t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query);
    solver.solve("time_forward_Euler_space_central");
    % [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = solver.getSolution();
    solver_path = solver.saveSnapShot();
end
```

#### DoSag2D.m

```matlab {.line-numbers}
%% DoSag2D.m
% Description: Solving the Streeter-Phelps model in a 2-D river
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Nov. 17, 2024
% Last modified: Dec. 4, 2024
%

%%

classdef DoSag2D < handle
    properties (Constant)
        PARAMS_REQUIRED = ["t_start", "range_x", "range_y", "velocity", "diffu_coeff", "K_r", "Q_B_loc", "Q_B", "K_a", "O_s", "init_func", "bndry_func_x0", "bndry_func_x1", "bndry_func_y0", "bndry_func_y1", "bndry_func_type_x0", "bndry_func_type_x1", "bndry_func_type_y0", "bndry_func_type_y1", "delta_t", "delta_x", "delta_y", "t_query"];
        SCHEME_NAME = ["time_forward_Euler_space_central"];
    end

    properties (Access = public)
        desc_str
    end

    properties (GetAccess = public, SetAccess = immutable)
        objID
        LOG_FILE_PATH
    end

    methods (Static, Access=public)
        params_struct = prepareParams(t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query);
    end

    methods (Access=public)
        function obj = DoSag2D(desc_str, varargin)
            arguments
                desc_str string = "description of this solver\n"
            end
            arguments (Repeating)
                varargin
            end

            obj.objID = "DoSag_" + string(datetime("now"), "yyyyMMdd_HHmmss_SS", "en_US");
            obj.desc_str = desc_str;

            if ~isfolder(".\log\")
                mkdir(".\log\")
            end
            obj.LOG_FILE_PATH = sprintf(".\\log\\%s.log", obj.objID);
            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tSolver created.\n\tDescription: %s\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), desc_str);
            fclose(LOG_FILE_ID);
            
            if nargin < 2
                obj.reset(DoSag2D.prepareParams());
                return
            end
            obj.reset(DoSag2D.prepareParams(varargin{:}));
        end

        obj = reset(obj, params_struct);
    
        [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = getSolution(obj);

        function file_path = saveSnapShot(obj, file_path)
            arguments
                obj
                file_path = ".\lib\" + obj.objID + ".mat";
            end
            [filepath, ~, ~] = fileparts(file_path);
            if ~isfolder(filepath)
                mkdir(filepath);
            end
            save(file_path, "obj");
        end

        function obj = solve(obj, scheme_name)
            arguments (Input)
                obj 
                scheme_name string = "time_forward_Euler_space_central";
            end

            assert(length(obj.velocity) == 2, "仅支持 2-D 问题");
            [obj.f_list, obj.x_grid, obj.y_grid, obj.t_list] = obj.solve2D(scheme_name);
            obj.saveSnapShot();
        end
    end

    methods (Access=private)
        f_steady = solve2DSteady(obj, scheme_name);

        [f_list, x_grid, y_grid, t_list] = solve2D(obj, scheme_name);

        B_src = bSrc(obj);
    end % end of `methods (Access=private)`

    properties (GetAccess = public, SetAccess=private)
        % solver params
        t_start
        range_x
        range_y
        velocity % [u, v]
        diffu_coeff % [D_x, D_y]
        K_r
        Q_B_loc % [x1,y1; x2, y2;] BOD 点源位置
        Q_B     % BOD 点源强度. [Q1; Q2;]
        K_a
        O_s
        init_func % {B_func; O_func}
        bndry_func_x0 % x = range_x(1) 边界 {B_func; O_func}
        bndry_func_x1 % x = range_x(2) 边界 right boundary
        bndry_func_y0 % y = range_y(1) 边界
        bndry_func_y1 % y = range_y(2) 边界
        bndry_func_type_x0 % x = range_x(1) 边界类型. "Dirichlet" | "Neumann"
        bndry_func_type_x1 % x = range_x(2) 边界类型.
        bndry_func_type_y0 % y = range_y(1) 边界类型.
        bndry_func_type_y1 % y = range_y(2) 边界类型.
        delta_t
        delta_x
        delta_y
        t_query

        % 中间结果
        B_src   % (i, j) 元素表示 x = x_grid{i}, y = y_grid{j} 格点上的 BOD 源项 (加载率): \sum_{(x_0, y_0)} Q(x_0, y_0) I_{i,j}(x_0, y_0) / (h_x h_y), 即该格点所代表的那个有限体积元中的源
        f_steady  % f_steady{1}, f_steady{2} 分别表示 BOD, DO 的稳态解

        % results
        x_grid
        y_grid
        t_list
        f_list  % f_list{k}{1}(i, j)、f_list{k}{2}(i, j) 分别是 t_list(k) 时刻 x = x_grid(i), y = y_grid(j) 格点的 BOD、DO 值.
        scheme_name
    end
end % end of `classdef DoSag2D < handle`
```

#### bSrc.m

```matlab {.line-numbers}
function B_src = bSrc(obj)
    % Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
    % Q_B = [70; 0];    % BOD 点源强度

    x_bin_edges = [obj.x_grid - obj.delta_x / 2; obj.x_grid(end) + obj.delta_x / 2];
    x_bin_ind = discretize(obj.Q_B_loc(:, 1), x_bin_edges);

    y_bin_edges = [obj.y_grid - obj.delta_y / 2; obj.y_grid(end) + obj.delta_y / 2];
    y_bin_ind = discretize(obj.Q_B_loc(:, 2), y_bin_edges);

    B_src = zeros([length(obj.x_grid), length(obj.y_grid)]);
    for src_ind = 1:length(obj.Q_B)
        B_src(x_bin_ind(src_ind), y_bin_ind(src_ind)) = B_src(x_bin_ind(src_ind), y_bin_ind(src_ind)) + obj.Q_B(src_ind) / obj.delta_x / obj.delta_y;
    end
end
```

#### getSolution.m

```matlab {.line-numbers}
function [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = getSolution(obj)
    if isempty(obj.f_list)
        error("Problem not solved yet. Please call obj.solve(scheme_name) before calling get_solution().")
    end

    B_list = arrayfun(@(t_ind) obj.f_list{t_ind}{1}, 1:length(obj.f_list), UniformOutput=false);
    O_list = arrayfun(@(t_ind) obj.f_list{t_ind}{2}, 1:length(obj.f_list), UniformOutput=false);
    x_grid = obj.x_grid;
    y_grid = obj.y_grid;
    t_list = obj.t_list;
    scheme_name = obj.scheme_name;
end
```

#### prepareParams.m

```matlab {.line-numbers}
function params_struct = prepareParams(t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query)
    arguments
        t_start = 0;
        range_x (1, 2) = [0, 150];
        range_y (1, 2) = [0, 30];
        velocity (1, 2) = [.4, 0];
        diffu_coeff (1, 2) = [.5, .5];
        K_r {isscalar} = .01;
        Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
        Q_B = [70; 0];    % BOD 点源强度
        K_a {isscalar} = .02;
        O_s = 8;
        init_func = {@(x, y) zeros(size(x)); @(x, y) O_s*ones(size(x))};
        bndry_func_type_x0 = "Dirichlet";
        bndry_func_x0 = {@(y, t) zeros(size(y)); @(y, t) O_s*ones(size(y))};
        bndry_func_type_x1 = "Neumann";
        bndry_func_x1 = {@(y, t) zeros(size(y)); @(y, t) zeros(size(y))};
        bndry_func_type_y0 = "Neumann";
        bndry_func_y0 = {@(x, t) zeros(size(x)); @(x, t) zeros(size(x))};
        bndry_func_type_y1 = bndry_func_type_y0;
        bndry_func_y1 = bndry_func_y0;
        delta_t {isscalar} = .01;
        delta_x {isscalar} = .2;
        delta_y {isscalar} = .2;
        t_query {isvector} = [0, 1, 2, +Inf];
    end

    % 求解区域
    params_struct.t_start = t_start;
    params_struct.range_x = sort(range_x, "ascend");
    params_struct.range_y = sort(range_y, "ascend"); 

    params_struct.velocity = velocity;
    params_struct.diffu_coeff = diffu_coeff;
    params_struct.K_r = K_r;
    params_struct.Q_B_loc = Q_B_loc;
    params_struct.Q_B = Q_B;
    params_struct.K_a = K_a;
    params_struct.O_s = O_s;

    params_struct.init_func = init_func;
    params_struct.bndry_func_type_x0 = bndry_func_type_x0;
    params_struct.bndry_func_x0 = bndry_func_x0;
    params_struct.bndry_func_type_x1 = bndry_func_type_x1;
    params_struct.bndry_func_x1 = bndry_func_x1;
    params_struct.bndry_func_type_y0 = bndry_func_type_y0;
    params_struct.bndry_func_y0 = bndry_func_y0;
    params_struct.bndry_func_type_y1 = bndry_func_type_y1;
    params_struct.bndry_func_y1 = bndry_func_y1;

    params_struct.delta_t = delta_t;
    params_struct.delta_x = delta_x;
    params_struct.delta_y = delta_y;
    params_struct.t_query = unique(t_query(t_query >= t_start), "sorted");
end
```

#### reset.m

```matlab {.line-numbers}
function obj = reset(obj, params_struct)
    field_missing_flag = ~isfield(params_struct, DoSag2D.PARAMS_REQUIRED);
    if any(field_missing_flag)
        warning("Invalid params_struct. Params not changed. Missing field: %s", join(DoSag2D.PARAMS_REQUIRED(field_missing_flag), ", "));
        return
    end

    params_struct = DoSag2D.prepareParams(params_struct.t_start, params_struct.range_x, params_struct.range_y, params_struct.velocity, params_struct.diffu_coeff, params_struct.K_r, params_struct.Q_B_loc, params_struct.Q_B, params_struct.K_a, params_struct.O_s, params_struct.init_func, params_struct.bndry_func_type_x0, params_struct.bndry_func_x0, params_struct.bndry_func_type_x1, params_struct.bndry_func_x1, params_struct.bndry_func_type_y0, params_struct.bndry_func_y0, params_struct.bndry_func_type_y1, params_struct.bndry_func_y1, params_struct.delta_t, params_struct.delta_x, params_struct.delta_y, params_struct.t_query);
    for field_name = DoSag2D.PARAMS_REQUIRED
        obj.(field_name) = params_struct.(field_name);
    end

    obj.x_grid = (obj.range_x(1): obj.delta_x: obj.range_x(2)).';
    obj.y_grid = (obj.range_y(1): obj.delta_y: obj.range_y(2)).';
    obj.B_src = obj.bSrc();

    obj.t_list = [];
    obj.f_steady = [];
    obj.f_list = [];
end
```

#### solve2D.m

```matlab {.line-numbers}
function [f_list, x_grid, y_grid, t_list] = solve2D(obj, scheme_name)
    arguments (Input)
        obj
        scheme_name = "time_forward_Euler_space_central";
    end

    x_grid = obj.x_grid;
    y_grid = obj.y_grid;
    [X, Y] = ndgrid(obj.x_grid, obj.y_grid);
    f_list = cell(size(obj.t_query));
    t_list = nan(size(obj.t_query));
    
    if ismember(scheme_name, "time_forward_Euler_space_central")
        obj.scheme_name = scheme_name;

        Cr = obj.velocity * obj.delta_t ./ [obj.delta_x, obj.delta_y];
        r = obj.diffu_coeff * obj.delta_t ./ [obj.delta_x, obj.delta_y].^2;
        
        f_last = cellfun(@(i_func) i_func(X, Y), obj.init_func, UniformOutput=false);
        % f_last{1} = obj.init_func{1}(X, Y); % BOD
        % f_last{2} = obj.init_func{2}(X, Y); % DO
        f_now = f_last;
        t_last = obj.t_start;
        t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

        % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
        if isempty(t_query_next_idx) || t_query_next_idx > 1
            t_list(1) = t_last;
            f_list{1} = f_last;

            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tInitial value added.\n\tt_start = %g\n", datetime('now'), obj.desc_str, scheme_name, t_last);
            fclose(LOG_FILE_ID);
            fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tInitial value added.\n\tt_start = %g\n", datetime('now'), obj.desc_str, scheme_name, t_last);
        end

        %%% 若请求 +Inf 时刻, 则创建稳态解, obj.f_steady
        if any(obj.t_query == +Inf)
            f_steady = obj.solve2DSteady(scheme_name);
            t_list(obj.t_query == +Inf) = +Inf;
            f_list{obj.t_query == +Inf} = f_steady;

            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tSteady state solved.\n\tt_steady = +Inf, DO_min (sat) = %.3g (%.3g),\n", datetime('now'), obj.desc_str, scheme_name, min(f_steady{2}, [], "all"), obj.O_s);
            fclose(LOG_FILE_ID);
            fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tSteady state solved.\n\tt_steady = +Inf, DO_min (sat) = %.3g (%.3g),\n", datetime('now'), obj.desc_str, scheme_name, min(f_steady{2}, [], "all"), obj.O_s);
        end

        %%% 下面求解有限时刻的解.

        T_END = obj.t_query(obj.t_query < Inf);
        T_END = T_END(end);
        
        % 用于记录进度
        REC_TIME_INT = 5; % seconds
        time_start = tic;
        t_now_last_rec = t_last;

        % 测试
        if ~isfolder(".\log\")
            mkdir(".\log\")
        end
        thres = [.5, (obj.O_s - 6) / (obj.O_s - min(f_steady{2}, [], "all")), .99, .999];
        flag_test = true(size(thres));

        while t_last < T_END
            t_now = t_last + obj.delta_t;

            % 输出进度
            if toc(time_start) > REC_TIME_INT
                test_prog = (obj.O_s - min(f_last{2}, [], "all")) / (obj.O_s - min(f_steady{2}, [], "all"));
                time_remained = REC_TIME_INT * (T_END - t_now) / (t_now - t_now_last_rec);
                fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), Cr = (%s), r = (%s), test_prog = %.2f %%,\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), obj.desc_str, scheme_name, time_remained, join(string(size(X)), ", "), join(string(Cr), ", "), join(string(r), ", "), test_prog*100, t_now, obj.t_query(t_query_next_idx), T_END);
                time_start = tic;
                t_now_last_rec = t_now;
            end
            
            % add halo points, advance one time step
            if scheme_name == "time_forward_Euler_space_central"
                switch obj.bndry_func_type_x0
                    case "Dirichlet"
                        halo_upper = cellfun(@(bndry_func, f) 2*bndry_func(obj.y_grid.', t_last) - f(2, :), obj.bndry_func_x0, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_upper = cellfun(@(bndry_func, f) f(2, :) - 2*obj.delta_x*bndry_func(obj.y_grid.', t_last), obj.bndry_func_x0, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_upper = cellfun(@(f) f(end-1, :), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_x1
                    case "Dirichlet"
                        halo_lower = cellfun(@(bndry_func, f) 2*bndry_func(obj.y_grid.', t_last) - f(end-1, :), obj.bndry_func_x1, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_lower = cellfun(@(bndry_func, f) f(end-1, :) + 2*obj.delta_x*bndry_func(obj.y_grid.', t_last), obj.bndry_func_x1, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_lower = cellfun(@(f) f(2, :), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_y0
                    case "Dirichlet"
                        halo_left = cellfun(@(bndry_func, f) 2*bndry_func(obj.x_grid, t_last) - f(:, 2), obj.bndry_func_y0, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_left = cellfun(@(bndry_func, f) f(:, 2) - 2*obj.delta_y*bndry_func(obj.x_grid, t_last), obj.bndry_func_y0, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_left = cellfun(@(f) f(:, end-1), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_y1
                    case "Dirichlet"
                        halo_right = cellfun(@(bndry_func, f) 2*bndry_func(obj.x_grid, t_last) - f(:, end-1), obj.bndry_func_y1, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_right = cellfun(@(bndry_func, f) f(:, end-1) + 2*obj.delta_y*bndry_func(obj.x_grid, t_last), obj.bndry_func_y1, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_right = cellfun(@(f) f(:, 2), f_last, UniformOutput=false) ;
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end
                
                f_last_expanded = cellfun(@(f_l, h_left, h_right, h_upper, h_lower) [[NaN; h_left; NaN], [h_upper; f_l; h_lower], [NaN; h_right; NaN]], f_last, halo_left, halo_right, halo_upper, halo_lower, UniformOutput=false);
                f_now{1} = (1 - 2*r(1) - 2*r(2) - obj.delta_t*obj.K_r) * f_last{1} ...
                         + (r(1) - Cr(1)/2) * f_last_expanded{1}(3:end, 2:end-1) ...
                         + (r(1) + Cr(1)/2) * f_last_expanded{1}(1:end-2, 2:end-1) ...
                         + (r(2) - Cr(2)/2) * f_last_expanded{1}(2:end-1, 3:end) ...
                         + (r(2) + Cr(2)/2) * f_last_expanded{1}(2:end-1, 1:end-2) ...
                         + obj.delta_t * obj.B_src;
                f_now{2} = (1 - 2*r(1) - 2*r(2) - obj.delta_t*obj.K_a)*f_last{2} ...
                + (r(1) - Cr(1)/2) * f_last_expanded{2}(3:end, 2:end-1) ...
                + (r(1) + Cr(1)/2) * f_last_expanded{2}(1:end-2, 2:end-1) ...
                + (r(2) - Cr(2)/2) * f_last_expanded{2}(2:end-1, 3:end) ...
                + (r(2) + Cr(2)/2) * f_last_expanded{2}(2:end-1, 1:end-2) ...
                + obj.delta_t * (obj.K_a * obj.O_s - obj.K_r * f_last{1});
            else
                error("格式 %s 正在开发", scheme_name)
            end

            % store solution at queried time.
            if ~((obj.t_query(t_query_next_idx) - t_last)*(obj.t_query(t_query_next_idx) - t_now) > 0)
                if abs(obj.t_query(t_query_next_idx) - t_now) < abs(obj.t_query(t_query_next_idx) - t_last)
                    t_list(t_query_next_idx) = t_now;
                    f_list{t_query_next_idx} = f_now;
                else
                    t_list(t_query_next_idx) = t_last;
                    f_list{t_query_next_idx} = f_last;
                end
                t_query_next_idx = t_query_next_idx + 1;
            end

            % prepare for next time step
            t_last = t_now;
            f_last = f_now;

            % 测试
            test_prog = (obj.O_s - min(f_last{2}, [], "all")) / (obj.O_s - min(f_steady{2}, [], "all"));
            for test_ind = 1:length(thres)
                if flag_test(test_ind) && test_prog > thres(test_ind)
                    flag_test(test_ind) = false;
                    LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
                    fprintf(LOG_FILE_ID, "%s\n\tt = %g, m = %.3g (m_s = %.3g), prog = %.3f %%\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), t_last, min(f_last{2}, [], "all"), min(f_steady{2}, [], "all"), 100*test_prog);
                    fprintf("%s\n\tt = %g, m = %.3g (m_s = %.3g), prog = %.3f %%\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), t_last, min(f_last{2}, [], "all"), min(f_steady{2}, [], "all"), 100*test_prog);
                    fclose(LOG_FILE_ID);
                end
            end
        end % end of `while t_now < obj.t_query(end)`

        LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
        fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), obj.desc_str, scheme_name, join(string(t_list), ", "));
        fclose(LOG_FILE_ID);
        fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), obj.desc_str, scheme_name, join(string(t_list), ", "));
        return
    end % end of `if ismember(scheme_name, "time_forward_Euler_space_central")`

    error("Invalid scheme_name. Please choose one of the following options: %s", join(DoSag2D.SCHEME_NAME, ", "));
end % end of `function [f_list, x_grid, t_list] = solve_2D(obj, scheme_name)`
```

#### solve2DSteady.m

```matlab {.line-numbers}
function f_steady = solve2DSteady(obj, scheme_name)
    arguments (Input)
        obj
        scheme_name = "time_forward_Euler_space_central";
    end

    if ~isempty(obj.f_steady)
        f_steady = obj.f_steady;
        return
    end
    f_steady = cell([2, 1]);

    if scheme_name == "time_forward_Euler_space_central"
        Cr = obj.velocity * obj.delta_t ./ [obj.delta_x, obj.delta_y];
        r = obj.diffu_coeff * obj.delta_t ./ [obj.delta_x, obj.delta_y].^2;

        %%% 1. 先求 B_steady
        
        % 1.1 创建方程组的系数矩阵 (sparse) 和常向量 (full)

        sz = [length(obj.x_grid), length(obj.y_grid)];
        c_vec{2} = nan(prod(sz), 1);
        c_vec{1} = c_vec{2};
        i_list{2} = nan(prod(sz - 2)*5 + sum(sz - 2)*6 + 8, 1);
        i_list{1} = i_list{2};
        j_list = i_list;
        v_list = i_list;
        cnt = 0;

        % 1.1.1 建立 halo_points
        switch obj.bndry_func_type_x0
            case "Dirichlet"
                halo_upper = {2, [2, -1]}; % halo point 是 边界函数 以及 法向格点值 的 线性组合, 记作 {边界函数 系数, [格点 1 索引, 系数 1; 格点 2 索引, 系数 2]}
            case "Neumann"
                halo_upper = {-2*obj.delta_x, [2, 1]};
            case "Periodic"
                % halo_upper = {0, [length(obj.x_grid)-1, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_x1
            case "Dirichlet"
                halo_lower = {2, [length(obj.x_grid)-1, -1]};
            case "Neumann"
                halo_lower = {2*obj.delta_x, [length(obj.x_grid)-1, 1]};
            case "Periodic"
                % halo_lower = {0, [2, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_y0
            case "Dirichlet"
                halo_left = {2, [2, -1]};
            case "Neumann"
                halo_left = {-2*obj.delta_y, [2, 1]};
            case "Periodic"
                % halo_left = {0, [length(obj.y_grid)-1, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_y1
            case "Dirichlet"
                halo_right = {2, [length(obj.y_grid)-1, -1]};
            case "Neumann"
                halo_right = {2*obj.delta_y, [length(obj.y_grid)-1, 1]};
            case "Periodic"
                % halo_right = {0, [2, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        % 1.1.2 建立 BOD 方程 (以后, 计划把 halo points 也写进方程的未知数, 便于统一各类边界条件的实现)
        for ind_x = 1:sz(1)
            for ind_y = 1:sz(2)
                ind_x_guess = ind_x + [0; 1; -1; 0; 0]; % 以后: 应加入用于表示四个 halo 的格点
                ind_y_guess = ind_y + [0; 0; 0; 1; -1];
                flag_invalid = ind_x_guess < 1 | ind_x_guess > sz(1) | ind_y_guess < 1 | ind_y_guess > sz(2);
                
                ind_x_valid = ind_x_guess(~flag_invalid);
                ind_y_valid = ind_y_guess(~flag_invalid);
                ind_lin = sub2ind(sz, ind_x_valid, ind_y_valid);
                c_vec{1}(ind_lin(1)) = obj.delta_t * obj.B_src(ind_x, ind_y) ...
                    + (ind_x == 1)*(r(1) + Cr(1)/2)*halo_upper{1}*obj.bndry_func_x0{1}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_x == sz(1))*(r(1) - Cr(1)/2)*halo_lower{1}*obj.bndry_func_x1{1}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_y == 1)*(r(2) + Cr(2)/2)*halo_left{1}*obj.bndry_func_y0{1}(obj.x_grid(ind_x), +Inf) ...
                    + (ind_y == sz(2))*(r(2) - Cr(2)/2)*halo_right{1}*obj.bndry_func_y1{1}(obj.x_grid(ind_x), +Inf);
                
                nnz_ = length(ind_lin);
                i_list{1}(cnt + 1 : cnt + nnz_) = ind_lin(1);
                j_list{1}(cnt + 1 : cnt + nnz_) = ind_lin;
                v_list_ = [ind_x_valid == ind_x & ind_y_valid == ind_y, ind_x_valid == ind_x+1, ind_x_valid == ind_x-1, ind_y_valid == ind_y+1, ind_y_valid == ind_y-1] * [ ...
                    2*sum(r) + obj.delta_t * obj.K_r; % B_{i,j}
                    -r(1) + Cr(1)/2; % B_{i+1, j}
                    -r(1) - Cr(1)/2; % B_{i-1, j}
                    -r(2) + Cr(2)/2; % B_{i, j+1}
                    -r(2) - Cr(2)/2  % B_{i, j-1}
                ];
                if ind_x == 1
                    % 以后: 用循环, 考虑用于表示 halo 的所有格点
                    v_list_(ind_x_valid == halo_upper{2}(1,1)) = v_list_(ind_x_valid == ind_x + 1) - (r(1) + Cr(1)/2)*halo_upper{2}(1,2);
                end
                if ind_x == sz(1)
                    v_list_(ind_x_valid == halo_lower{2}(1,1)) = v_list_(ind_x_valid == ind_x - 1) - (r(1) - Cr(1)/2)*halo_lower{2}(1,2);
                end
                if ind_y == 1
                    v_list_(ind_y_valid == halo_left{2}(1,1)) = v_list_(ind_y_valid == ind_y + 1) - (r(2) + Cr(2)/2)*halo_left{2}(1,2);
                end
                if ind_y == sz(2)
                    v_list_(ind_y_valid == halo_right{2}(1,1)) = v_list_(ind_y_valid == ind_y - 1) - (r(2) - Cr(2)/2)*halo_right{2}(1,2);
                end
                v_list{1}(cnt + 1 : cnt + nnz_) = v_list_;

                cnt = cnt + nnz_;
            end
        end

        if cnt ~= length(v_list{1})
            error("这是不可能的, 除非程序有误")
        end

        % 1.1.3 求解 BOD 稳态方程

        B_steady = sparse(i_list{1}, j_list{1}, v_list{1}, prod(sz), prod(sz)) \ c_vec{1};
        f_steady{1} = reshape(B_steady, sz);

        %%% 2. 再求 O_steady

        i_list{2} = i_list{1};
        j_list{2} = j_list{1};

        cnt = 0;
        for ind_x = 1:sz(1)
            for ind_y = 1:sz(2)
                ind_x_guess = ind_x + [0; 1; -1; 0; 0]; % 以后: 应加入用于表示四个 halo 的格点
                ind_y_guess = ind_y + [0; 0; 0; 1; -1];
                flag_invalid = ind_x_guess < 1 | ind_x_guess > sz(1) | ind_y_guess < 1 | ind_y_guess > sz(2);
                
                ind_x_valid = ind_x_guess(~flag_invalid);
                ind_y_valid = ind_y_guess(~flag_invalid);
                ind_lin = sub2ind(sz, ind_x_valid, ind_y_valid);
                c_vec{2}(ind_lin(1)) = obj.delta_t * (obj.K_a * obj.O_s - obj.K_r * f_steady{1}(ind_x, ind_y)) ...
                    + (ind_x == 1)*(r(1) + Cr(1)/2)*halo_upper{1}*obj.bndry_func_x0{2}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_x == sz(1))*(r(1) - Cr(1)/2)*halo_lower{1}*obj.bndry_func_x1{2}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_y == 1)*(r(2) + Cr(2)/2)*halo_left{1}*obj.bndry_func_y0{2}(obj.x_grid(ind_x), +Inf) ...
                    + (ind_y == sz(2))*(r(2) - Cr(2)/2)*halo_right{1}*obj.bndry_func_y1{2}(obj.x_grid(ind_x), +Inf);
                
                nnz_ = length(ind_lin);
                v_list_ = [ind_x_valid == ind_x & ind_y_valid == ind_y, ind_x_valid == ind_x+1, ind_x_valid == ind_x-1, ind_y_valid == ind_y+1, ind_y_valid == ind_y-1] * [ ...
                    2*sum(r) + obj.delta_t * obj.K_a; % O_{i,j}
                    -r(1) + Cr(1)/2; % O_{i+1, j}
                    -r(1) - Cr(1)/2; % O_{i-1, j}
                    -r(2) + Cr(2)/2; % O_{i, j+1}
                    -r(2) - Cr(2)/2  % O_{i, j-1}
                ];
                if ind_x == 1
                    % 以后: 用循环, 考虑用于表示 halo 的所有格点
                    v_list_(ind_x_valid == halo_upper{2}(1,1)) = v_list_(ind_x_valid == ind_x + 1) - (r(1) + Cr(1)/2)*halo_upper{2}(1,2);
                end
                if ind_x == sz(1)
                    v_list_(ind_x_valid == halo_lower{2}(1,1)) = v_list_(ind_x_valid == ind_x - 1) - (r(1) - Cr(1)/2)*halo_lower{2}(1,2);
                end
                if ind_y == 1
                    v_list_(ind_y_valid == halo_left{2}(1,1)) = v_list_(ind_y_valid == ind_y + 1) - (r(2) + Cr(2)/2)*halo_left{2}(1,2);
                end
                if ind_y == sz(2)
                    v_list_(ind_y_valid == halo_right{2}(1,1)) = v_list_(ind_y_valid == ind_y - 1) - (r(2) - Cr(2)/2)*halo_right{2}(1,2);
                end
                v_list{2}(cnt + 1 : cnt + nnz_) = v_list_;

                cnt = cnt + nnz_;
            end
        end

        if cnt ~= length(v_list{2})
            error("这是不可能的, 除非程序有误")
        end

        O_steady = sparse(i_list{2}, j_list{2}, v_list{2}, prod(sz), prod(sz)) \ c_vec{2};
        f_steady{2} = reshape(O_steady, sz);

        %%% 3.

        obj.f_steady = f_steady;
    else
        error("scheme_name: %s 的稳态解算法未实现", scheme_name)
    end
end
```
