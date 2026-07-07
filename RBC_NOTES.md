# 混合频率半参数回归中的稳健偏差校正（RBC）技术说明

> 本文档说明 `mfxtsemipar_cv3` 中 bias correction 的理论动机、两种 RBC 形式（`bc1`、`bc2`）以及为什么必须采用 stacked regression 来估计标准误。重点在原理，不在代码实现。

---

## 1. 背景：混合频率数据与半参数估计

在混合频率（mixed-frequency）设定中，某个解释变量 $u$ 在高频数据上观测，而被解释变量 $y$ 只在低频层面有定义。典型的模型可以写成

$$
y_{it} = x_{it}'\beta + g(u_{it}) + \alpha_i + \delta_t + \varepsilon_{it},
$$

其中 $g(\cdot)$ 是未知的光滑函数，$\alpha_i$、$\delta_t$ 分别是个体/时间固定效应。为了把高频的 $u$ 用于低频的 $y$，通常先把 $u$ 的某种“基函数变换”在高频层面加总到低频单元，再对低频样本做回归。

用样条（splines）或更一般的 partitioning-based series basis $B_K(u)$ 来逼近 $g(u)$，就得到级数/样条估计量

$$
\hat g_K(u) = B_K(u)' \hat\beta,
\qquad
\hat\beta = \arg\min_b \sum_{\ell}\bigl(y_\ell - B_K(u_\ell)'b - \text{controls}_\ell\bigr)^2.
$$

$K$ 表示基函数的维数（节点数 + 多项式阶数等）。$K$ 越大，逼近偏差越小，但估计方差越大。

---

## 2. 为什么需要偏差校正？

### 2.1 偏差来源

对 $g(u)$ 的最小二乘级数估计可以分解为

$$
\hat g_K(u) - g(u)
= \underbrace{B_K(u)'\hat\beta - B_K(u)'\beta}_{\text{抽样误差}}
+ \underbrace{B_K(u)'\beta - g(u)}_{\text{逼近偏差}}.
$$

第一项在 $K$ 固定时有标准的渐近正态分布；第二项是模型误设/逼近误差，其阶数取决于 $g$ 的光滑性、基的阶数和节点数。对于“IMSE 最优”的 $K$ 选择，偏差与方差通常同阶，导致基于 $\hat g_K(u)$ 的常规置信区间覆盖不足（undercover）。

### 2.2 传统做法的缺陷

两种常见做法都有问题：

1. **undersmoothing**：选比 IMSE 最优更小的 $K$，让偏差渐近可忽略。这会损失效率，且如何选择“足够小”的 $K$ 缺乏明确标准。
2. **解析偏差估计**：需要知道 $g$ 的未知高阶导数，并且依赖渐近展开，有限样本表现不稳定。

稳健偏差校正（Robust Bias Correction, RBC）的思路是：**保持 IMSE 最优的点估计，但用另一个高维/高阶基估计偏差项，从而构造有效的推断。**

---

## 3. 稳健偏差校正（RBC）原理

本节主要参考 Cattaneo, Farrell and Feng (2020a) 与 `lspartition` 软件包（Cattaneo, Farrell and Feng, 2020b）。

### 3.1 记法

设主回归使用的基为 $B(u)$（维数 $K$），偏差校正使用的更富裕基为 $\tilde B(u)$（维数 $\tilde K > K$，通常是更高阶或更多节点）。用同样的记号表示低频聚合后的设计矩阵 $B$ 与 $\tilde B$。可能还有权重矩阵 $W$（如加权最小二乘或高频加总权重）。

### 3.2 bc1：直接用高阶基重新估计

`bc1` 直接用 $\tilde B$ 代替 $B$ 做一次新的级数回归：

$$
\hat g_{\text{bc1}}(u) = \tilde B(u)' \hat\gamma,
\qquad
\hat\gamma = (\tilde B' W \tilde B)^{-1} \tilde B' W y.
$$

直观上，$\tilde B$ 比 $B$ 更灵活，因而逼近偏差更小。Cattaneo, Farrell and Feng (2020a) 证明，在适当的 rate 条件下，$\hat g_{\text{bc1}}(u)$ 的偏差阶数低于原始估计量，因而可以构造有效的置信区间。

`bc1` 的优点是概念简单；缺点是它完全抛弃原始基 $B$ 的估计结果，只用高阶基，方差通常比 `bc2` 大。

### 3.3 bc2：最小二乘偏差校正

`bc2` 的想法是：保留原始估计量 $\hat g_K(u)=B(u)'\hat\beta$，但用一个基于 $\tilde B$ 的修正项来估计并扣除偏差。具体地，先把高阶基投影到原始基空间：

$$
\hat P = (B' W B)^{-1} B' W \tilde B.
$$

然后定义校正后的预测矩阵为

$$
W_{\text{bc2}}(u) = \bigl[B(u),\; \tilde B(u) - B(u)\hat P\bigr].
$$

`bc2` 估计量为

$$
\hat g_{\text{bc2}}(u)
= B(u)'\hat\beta + \bigl[\tilde B(u) - B(u)\hat P\bigr]' \hat\gamma
= W_{\text{bc2}}(u)'
\begin{pmatrix}\hat\beta\\ \hat\gamma\end{pmatrix}.
$$

这里的 $\hat\gamma$ 来自高阶基 $\tilde B$ 的回归系数。$\tilde B(u)-B(u)\hat P$ 是把高阶基中“已被原始基解释掉”的部分去掉后剩余的修正方向；乘以 $\hat\gamma$ 就得到对原始估计偏差的估计。

Cattaneo, Farrell and Feng (2020a) 把这类估计量称为基于最小二乘投影的稳健偏差校正。在适当的正则条件下，`bc2` 的偏差阶数与 `bc1` 相同，但方差通常更小或至少不更大。

### 3.4 bc1 与 bc2 的区别

- `bc1` 完全依赖高阶基 $\tilde B$ 重新拟合；
- `bc2` 在原估计 $B\hat\beta$ 基础上加上一个“去重后的高阶修正”。

二者在渐近上等价于某些特定的正交化构造，但在有限样本中并不相等。在我们的实现中它们被同时报告，研究者可以根据覆盖率和效率权衡选择。

---

## 4. 为什么必须用 stacked regression 估计标准误？

### 4.1 独立回归的问题

`bc2` 的公式同时用到了 $\hat\beta$ 与 $\hat\gamma$。如果分别跑两个独立回归——一个用 $B$、一个用 $\tilde B$——确实可以得到相同的点估计，但**无法得到正确的标准误**：

- `bc2` 是 $\hat\beta$ 与 $\hat\gamma$ 的线性组合；
- 其方差需要完整的协方差矩阵
  $$
  \text{Var}\begin{pmatrix}\hat\beta\\ \hat\gamma\end{pmatrix}
  =
  \begin{pmatrix}
  \Sigma_{\beta\beta} & \Sigma_{\beta\gamma}\\
  \Sigma_{\gamma\beta} & \Sigma_{\gamma\gamma}
  \end{pmatrix}.
  $$

两个独立回归只能给出 $\Sigma_{\beta\beta}$ 和 $\Sigma_{\gamma\gamma}$，而缺少 $\Sigma_{\beta\gamma}$。忽略协方差会导致 `bc2` 的标准误有偏。

### 4.2 Stacked regression 的等价性

为了同时估计 $\hat\beta$、$\hat\gamma$ 以及它们的协方差，可以把两个方程纵向堆叠（stacked regression）。

构造一个方程指示变量 $\text{eq}$：

- 当 $\text{eq}=0$（主方程）时，被解释变量为 $y$，解释变量为 $B$ 和控制变量；
- 当 $\text{eq}=1$（偏差校正方程）时，被解释变量同样为 $y$，解释变量为 $\tilde B$ 和同样的控制变量。

把所有方程 0 的变量在方程 1 中置零，反之亦然，得到如下联合回归：

$$
y = B_{\text{eq}=0}\beta + \tilde B_{\text{eq}=1}\gamma + X_{\text{eq}=0}\delta_0 + X_{\text{eq}=1}\delta_1 + \text{FE}_{\text{eq}=0} + \text{FE}_{\text{eq}=1} + \text{error}.
$$

关键点是：

1. **主方程与 BC 方程的固定效应必须不同**：把固定效应与 $\text{eq}$ 交乘，例如 `id^eq`、`time^eq`，这样两个方程各估各的 FE，等价于独立回归。
2. **控制变量也要按方程分开**：两个方程中对控制变量的系数允许不同。
3. **聚类变量保持原变量**：每个低频单元在堆叠数据中出现两次，按原聚类变量聚类可以让 cluster-robust 方差估计自动包含跨方程的协方差 $\Sigma_{\beta\gamma}$。

### 4.3 标准误的 delta method

得到堆叠回归的协方差矩阵 $\hat\Sigma$ 后，任意评估点 $u$ 的 `bc2` 估计量可以写成

$$
\hat g_{\text{bc2}}(u) = W_{\text{bc2}}(u)' \hat\theta,
\qquad
\hat\theta = (\hat\beta', \hat\gamma')'.
$$

其标准误为

$$
\widehat{\text{SE}}\bigl(\hat g_{\text{bc2}}(u)\bigr)
= \sqrt{W_{\text{bc2}}(u)' \hat\Sigma W_{\text{bc2}}(u)}.
$$

因为 $\hat\Sigma$ 来自堆叠回归，它同时反映 $B$ 和 $\tilde B$ 系数的不确定性及其协方差，从而给出正确的推断。

---

## 5. 与 `lspartition` 的对应关系

`lspartition`（Cattaneo, Farrell and Feng, 2020b）实现了 partitioning-based least squares 的非参数回归，并内置三种 RBC 策略。我们的实现借鉴其思想，但针对混合频率半参数回归做了调整：

- 主回归是在低频单元上进行的，但 $B$ 与 $\tilde B$ 都在高频层面构造后再加总；
- 偏差校正基 $\tilde B$ 的阶数通常取 `degree + 1`，与 `lspartition` 中通过 `m.bc` 提升阶数的做法一致；
- `bc1` / `bc2` 的定义与 `lspartition` 中的 `bc` 选项对应；
- 标准误通过堆叠回归一次性获得，而不是对两个独立回归结果做后期拼接。

---

## 6. 主要参考文献

1. **Cattaneo, M. D., M. H. Farrell, and Y. Feng (2020a).** “Large Sample Properties of Partitioning-Based Series Estimators.” *Annals of Statistics*, 48(3), 1718–1741.  
   —— 系列/样条 RBC 的理论基础，给出 `bc1` 与 `bc2` 的渐近性质。

2. **Cattaneo, M. D., M. H. Farrell, and Y. Feng (2020b).** “lspartition: Partitioning-Based Least Squares Regression.” *R Journal*, 12(1), 172–187.  
   —— `lspartition` 软件包介绍，包含实现细节与经验示例。

3. **Cattaneo, M. D., and M. H. Farrell (2013).** “Optimal Convergence Rates, Bahadur Representation, and Asymptotic Normality of Partitioning Estimators.” *Journal of Econometrics*, 174(2), 127–143.  
   —— partitioning-based 估计量的收敛速度与 Bahadur 表示。

4. **Calonico, S., M. D. Cattaneo, and M. H. Farrell (2018).** “On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference.” *Journal of the American Statistical Association*, 113(522), 767–779.  
   —— RBC 在核回归/RD 中的覆盖率分析。

5. **Calonico, S., M. D. Cattaneo, and M. H. Farrell (2022).** “Coverage Error Optimal Confidence Intervals for Local Polynomial Regression.” *Bernoulli*, 28(4), 2998–3022.  
   —— 局部多项式框架下覆盖误差最优的 RBC 置信区间。

6. **Andrews, D. W. K. (1991).** “Asymptotic Normality of Series Estimators for Nonparametric and Semiparametric Regression Models.” *Econometrica*, 59(2), 307–345.  
   —— 级数估计量渐近正态性的经典文献。

7. **Newey, W. K. (1997).** “Convergence Rates and Asymptotic Normality for Series Estimators.” *Journal of Econometrics*, 79(1), 147–168.  
   —— 级数估计的收敛速度与标准误理论。

---

## 7. 一句话总结

> 在样条级数估计中，IMSE 最优节点选择会导致不可忽略的逼近偏差。RBC 用一个更高阶/更富裕的基来估计并校正这一偏差；`bc1` 直接改用高阶基，`bc2` 在原估计上加入最小二乘投影修正。为了正确计算 `bc2` 的标准误，必须通过 stacked regression 同时估计两个方程，并把固定效应、控制变量按方程分开，聚类变量保持原变量，以保留系数间的协方差。

---

## 补充：Stacked Regression 与 Joint Regression 的关系

### CFF 的原始写法是 joint regression

Cattaneo, Farrell and Feng (2020a) 的理论构造本质上是一个**增广最小二乘（augmented/joint regression）**：把低阶基 $B$ 与高阶/偏差校正基 $\tilde B$ 放在同一个回归方程里，得到联合系数 $(\hat\beta, \hat\gamma)$，再用同一个方差–协方差矩阵做 delta-method 推断。`lspartition::lsprobust()` 也是按这一思路实现的（`bc1`、`bc2`、`bc3` 都基于同一次拟合）。因此，**如果目标是严格复现 CFF/`lspartition`，joint regression 更直接、更一致。**

### Stacked regression 不是 CFF 的原始写法，但可以作为等价实现

Stacked regression 把两个方程纵向堆叠：

- `eq=0`：主方程，使用 $B$ 和控制变量；
- `eq=1`：偏差校正方程，使用 $\tilde B$ 和相同的控制变量；
- 固定效应、控制变量都与 `eq` 交乘，聚类变量保持原变量。

这种写法在以下意义下是合理的：

1. **对 `bc1` 完全等价**：`bc1` 的定义就是单独用 $\tilde B$ 做一次级数回归，stacked regression 的 `eq=1` 正好就是这个回归。
2. **对 `bc2` 渐近等价**：`bc2` 的校正项 $[\tilde B(u)-B(u)\hat P]'\hat\gamma$ 估计的是同一 leading bias。只要 $\tilde B$ 足够灵活、$K$ 与 $\tilde K$ 满足 CFF 的 rate 条件，stacked 版 `bc2` 与 joint 版 `bc2` 具有相同的渐近分布。
3. **标准误可以通过聚类得到**：按原聚类变量聚类后，cluster-robust VCV 会保留 $\hat\beta$ 与 $\hat\gamma$ 之间的协方差，从而 delta-method 标准误是正确的。

### 但两者在有限样本中不完全相同

需要注意，stacked regression 与 joint regression 的**系数参数化不同**，finite-sample 表现可能不一样：

- **共线性处理方式不同**：对 B-spline 等基，$\tilde B$ 的列空间通常包含 $B$ 的列空间。joint regression 如果把 $B$ 与 $\tilde B$ 同时放入同一方程，会出现完全共线性；`lspartition` 通过 `proj=TRUE` 把高阶基的 leading approximation error 投影到低阶空间后再进入回归。Stacked regression 则把 $B$ 与 $\tilde B$ 放在不同方程，天然避免了同方程内的共线性，但这也意味着 $\hat\gamma$ 来自独立的 `y ~ \tilde B` 回归，而不是 joint regression 中经过正交化后的系数。
- **标准误的估计方式不同**：joint regression 使用同一次回归的 EHW/HC 方差；stacked regression 使用跨方程聚类的 robust 方差。两者是同一渐近方差的不同有限样本估计，数值上可能有差异。

### 实际建议

- **如果研究目标是和 CFF/`lspartition` 严格保持一致**（例如审稿人要求完全对应 `lspartition` 的 `bc2`/`bc3`），应采用 joint regression，并像 `lspartition` 那样处理投影/共线性问题（对 splines 尤其要注意 `proj=TRUE` 或改用 `bc3` plug-in）。
- **如果研究目标是保留原始主回归的系数解释，同时做一个渐近等价的 RBC**，stacked regression 是 valid 的。在我们的实现中，stacked regression 的主方程系数与原 `mfxtsemipar` 独立回归完全相同，这是它的主要优势。
- **折中方案**：也可以先独立跑主回归得到 $\hat\beta$，然后把 $\tilde B$ 关于 $B$、控制变量和固定效应取残差，再与 $B$ 一起进入一次 joint regression 估计 $\hat\gamma$。这样 $B$ 的系数仍等于原主回归系数，而 `bc2` 的构造更接近 CFF 的 joint regression。不过这一步会重新引入同方程共线性的处理，需要额外小心。

### 一句话总结

> **CFF 的方法论基于 joint regression；stacked regression 是后者的一个有效但非原始的实现，对 `bc1` 完全等价，对 `bc2` 渐近等价，finite-sample 细节不同。是否“更一致”取决于你的目标：严格复现 CFF 选 joint regression；保持原独立回归系数解释选 stacked regression。**

---

## 补充二：实现选项 `bc_est`

在 `mfxtsemipar_cv3` 中，可以通过 `bc_est` 选择两种实现方式：

- `bc_est = "stacked"`（默认）：纵向堆叠两个方程，主方程用 $B$，BC 方程用 $\tilde B$，FE/控制变量与方程指示交乘。`bc1` 与 `bc2` 通常不同；系数与分别独立跑主回归、BC 回归完全相同。

- `bc_est = "joint"`：单次正交化联合回归，方程为 $y \sim B + (\tilde B - B P) + X + \text{FE}$。由于 $\tilde B - B P$ 与 $B$ 正交，$B$ 的系数仍等于原主回归系数；`bc1` 与 `bc2` 在此参数化下数值相同（都等于联合拟合值/高阶预测）。当 $\tilde B$ 与 $B$ 存在共线性时，`fixest` 会自动 drop 部分 residualized BC 列，这是正常的。

两种方法在适当的正则条件下都是渐近有效的；选择哪一种取决于你更看重“与独立回归系数完全一致”（`stacked`）还是“与 CFF/`lspartition` 的单次回归实现风格一致”（`joint`）。

---

## 补充三：Bootstrap 标准误

`mfxtsemipar_cv3` 还可以通过 `brep` 参数计算 wild-bootstrap 标准误：

```r
mfxtsemipar_cv3(..., brep = 200)
```

- `brep = 0`（默认）：使用解析的 cluster-robust 标准误。
- `brep > 0`：额外计算 wild-bootstrap 标准误，输出到 `gen_bc1_boot_se`（`bc1`）和 `gen_bc_boot_se`（`bc2`），同时返回 `vcov_bc1_boot` 和 `vcov_bc_boot`。

### 两种 `bc_est` 下的 bootstrap 策略

**`bc_est = "stacked"`**：
- `bc1` 的 bootstrap SE 来自对单独 BC 回归 `y ~ \tilde B + X + FE` 的 wild bootstrap。
- `bc2` 的 bootstrap SE 来自对堆叠回归的 wild bootstrap，以保留 $\hat\beta$ 与 $\hat\gamma$ 之间的协方差。

**`bc_est = "joint"`**：
- `bc1` 与 `bc2` 使用同一次正交化联合回归的 bootstrap VCV，因为二者在该参数化下是同一拟合值。

扰动方式与 `mfxtsemipar_cv.R` 一致：有聚类变量时按聚类抽 Rademacher 符号，无聚类时按观测抽。

---

## 补充四：调节 BC 模型的灵活性

除了默认的 `bc_degree = degree + 1`、同节点数之外，`mfxtsemipar_cv3` 还提供两个参数来进一步增加 BC 模型的灵活性：

- `bc_degree`：BC 样条的多项式阶数。默认比主模型高一阶；如果真实 $g(u)$ 非常光滑，可以设得更高（例如 `bc_degree = degree + 2`）。
- `bc_nknots`：BC 样条的内部节点数。默认 `NULL` 表示和主模型相同；设为更大的整数可以让 BC 模型在更多节点上更灵活。

示例：

```r
mfxtsemipar_cv3(
  ...,
  degree = 2,
  bc_degree = 4,   # 更高阶
  bc_nknots = 6    # 更多节点
)
```

这对应 `lspartition` 里的 `m.bc` 和 `bnknot`（`same = FALSE`）选项。增加 `bc_degree` 或 `bc_nknots` 都会让 BC 空间更大、偏差更小，但也会增加方差；标准误（解析或 bootstrap）会相应反映这一点。

---

## 补充五：CV 选择的是 `g(u)` 的 IMSE 最优吗？

### 精确表述：CV 最小化的是 `M g(u)` 的 IMSE

设 $M$ 为把控制变量 $x$、个体固定效应 $\alpha_i$ 和时间固定效应 $\delta_t$ 全部 partial out 的投影矩阵（即 Frisch–Waugh–Lovell 残差化矩阵）。对低频样本应用 $M$ 后，模型变为

$$
\tilde y = M y \approx M g(u) + M\varepsilon.
$$

用样条基 $B_K(u)$ 对 $\tilde y$ 做回归并做交叉验证时，第 $k$ 折的验证残差平方和为

$$
\bigl\| M y^{(k)} - M\hat g_K^{(-k)}(u^{(k)}) \bigr\|^2,
$$

其中 $\hat g_K^{(-k)}$ 在训练折上用 partial-out 后的数据估计。对完整 CV 目标取期望并忽略误差交叉项，最小化的对象是

$$
\bigl\| M g(u) - M\hat g_K(u) \bigr\|^2.
$$

也就是说，**CV 选出来的是在给定 partial-out 矩阵 $M$ 下，让 $M\hat g(u)$ 最接近 $M g(u)$ 的节点数**，而不是无条件 $g(u)$ 的 IMSE 最优。你的判断是对的，但需要更精确：它是在 $M$ 已经把 $x$ 和 FE 去掉之后，比较 partial-out 后的真实 $g$ 与 partial-out 后的预测 $g$。

### 什么时候 $M g(u) \approx g(u)$？

如果控制变量 $x$ 与 $u$ 不相关，且固定效应在 $u$ 的支集上大致均匀，则 $M$ 对 $g(u)$ 的影响很小，此时 CV 近似最小化 $g(u)$ 的 IMSE。但如果 $x$ 与 $u$ 相关（例如 $u$ 的某些区域被 $x$ 系统性解释），则 $M g(u)$ 会与 $g(u)$ 有差异：CV 会更重视那些“在 partial out $x$ 之后仍然难以解释”的 $u$ 区域。

这与 Robinson (1988) 的部分线性模型精神一致：非参数部分是在已经把参数部分去掉之后的残差中估计的，因此光滑参数选择也自然针对**条件非参数分量** $M g(u)$ 而非原始 $g(u)$。

### 与 CFF (2020a) 的 plug-in IMSE 选择的区别

Cattaneo, Farrell and Feng (2020a) 的 IMSE-optimal knot 选择是基于纯非参数回归

$$
y_i = g(u_i) + \varepsilon_i
$$

推导的 plug-in 公式，其中需要估计 $g$ 的高阶导数和误差方差（参见 `lspartition` 的 `k.select = 'imse-rot'` 或 `imse-dpi`）。这个公式在我们的面板设定里会失效或不直接适用，原因包括：

1. **误差结构更复杂**：原始误差 $\varepsilon_{it}$ 之外还有个体/时间固定效应；partial out 后残差 $\tilde\varepsilon_{it}$ 的方差结构与被吸收的 FE 有关，不再是 i.i.d. 噪声。
2. **混合频率加总**：低频 $y$ 对应多个高频 $u$ 的聚合；IMSE 的定义需要在加总后的低频单元上进行，CFF 的标量 $u_i$ 公式不直接对应。
3. **协变量影响**：如果 $x$ 与 $u$ 相关，partial out 步骤会改变有效样本的协方差结构，从而影响 IMSE 的偏-方差权衡。
4. **聚类结构**：面板数据通常存在聚类相关，plug-in IMSE 需要针对聚类误差重新推导，目前没有现成公式。

因此，**直接套用 `lspartition` 的 plug-in IMSE 节点选择到面板混合频率设定是站不住脚的。**

### CV 的有限样本解释

CV 虽然目标函数正确，但它仍然是样本层面的经验准则，不是总体 IMSE：

- 它对折（fold）的划分敏感；
- 当 $K$ 候选集合很小时，CV 曲线可能很平坦，导致节点选择不稳定；
- 它最小化的是预测 MSE，而不是理论上更精致的 coverage-error optimal 或 IMSE-optimal 准则。

但在没有可靠 plug-in 的复杂面板设定下，**CV 是一种可行的、以数据驱动方式逼近 $\bigl\|M g(u) - M\hat g_K(u)\bigr\|^2$ 最小的方案**。这也与 Robinson (1988) 的 partialing-out 精神一致：先把参数部分去掉，再对剩余的非参数部分用标准方法选择光滑参数。

### 参考 Robinson (1988)、FWL 与 CV 文献

- **Robinson, P. M. (1988).** “Root-$N$-Consistent Semiparametric Regression.” *Econometrica*, 56(4), 931–954.  
  —— 部分线性模型中 partial out 参数部分、再估计非参数部分的经典框架。

- **Frisch, R., and F. V. Waugh (1933).** “Partial Time Regressions as Compared with Individual Trends.” *Econometrica*, 1(4), 387–401；以及 **Lovell, M. C. (1963).** “Seasonal Adjustment of Economic Time Series.” *Journal of the American Statistical Association*, 58(304), 993–1010.  
  —— Frisch–Waugh–Lovell 定理：partial out 参数部分后，非参数系数的估计等价于在残差化后的变量上回归。

- **Härdle, W., P. Hall, and J. S. Marron (1988).** “How Far Are Automatically Chosen Regression Smoothing Parameters from Their Optimum?” *Journal of the American Statistical Association*, 83(404), 86–95.  
  —— 交叉验证在非参数回归中的最优性与收敛速度。

- **Li, Q., and J. S. Racine (2007).** *Nonparametric Econometrics: Theory and Practice*. Princeton University Press.  
  —— 第 3–4 章讨论 CV、IMSE 与样条节点选择。

- **Cattaneo, Farrell and Feng (2020a)**：前面已引用，给出 plug-in IMSE 的理论公式，但适用于纯截面非参数设定。

### 一句话总结

> `mfxtsemipar_cv3` 的 CV 目标函数是在 partial-out 矩阵 $M$ 作用之后，最小化 $\bigl\|M g(u) - M\hat g_K(u)\bigr\|^2$。当 $x$ 与 $u$ 大致不相关时，这近似等于 $g(u)$ 的 IMSE 最优；否则它针对的是“在 partial out 其他因素后，由样条预测的 $g$ 与真实 $g$ 的差距”。CFF (2020a) 的 plug-in IMSE 公式不能直接搬到面板混合频率设定，CV 是在这种复杂数据结构下一个合理且可操作的替代。
