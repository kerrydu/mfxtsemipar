# `mfxtsemipar2_cv` 模型与 R 代码说明

## 1. 模型背景

`mfxtsemipar2_cv` 用于估计**混合频率（mixed-frequency）半参数面板模型**，并同时得到：

- **短期反应曲线**（short-run response curve）
- **长期反应曲线**（long-run response curve）

典型数据结构中：

- `hf`（high-frequency）：高频数据，例如日内、日度观测，包含面板标识 `id`、时间层级 `tl`（如年、月）以及半参数变量 `uvar`。
- `lf`（low-frequency）：低频数据，每个 `id × tl` 只有一行，包含被解释变量 `y` 和低频控制变量 `x`。

高频样条基在 `hf` 上构造，然后按 `id × tl` 汇总到 `lf`，再与低频数据合并后进行回归。

---

## 2. 模型设定

对于低频观测 $(i, t)$，模型可写成：

$$
y_{it} = X_{it}'\beta + B_{it}'\theta + \bar{B}_i'\delta + \bar{X}_i'\eta + \alpha'\text{FE} + \varepsilon_{it}
$$

其中：

- $B_{it} = \sum_{j \in (i,t)} B(u_{ij})$：高频样条基在 `tl` 内的求和。
- $\bar{B}_i = \frac{1}{T_i}\sum_t B_{it}$：个体 $i$ 的样条基时间均值。
- $\bar{X}_i$：控制变量的个体均值。
- $\text{FE}$：通过 `absorb` 吸收的其他固定效应（例如时间 FE）。

### 2.1 短期与长期反应曲线

- **短期反应曲线**：$g_{SR}(u) = B(u)'\theta$
  - 反映同一单位内，$u$ 的当期变化对 $y$ 的影响。
- **长期反应曲线**：$g_{LR}(u) = B(u)'(\theta + \delta)$
  - 反映跨单位的、持久的 $u$ 差异对 $y$ 的影响。

通过引入个体均值 $\bar{B}_i$，模型把“组内（within）”和“组间（between）”效应分离开来。这与 Mundlak / correlated random effects 的思想一致。

---

## 3. 个体固定效应的处理

与 `mfxtsemipar_cv`（把 `id` 直接放入 `absorb`）不同，`mfxtsemipar2_cv` **不把 `id` 作为固定效应吸收**，而是：

1. 对每个回归变量（样条基、主回归控制变量、partial-out 变量）生成其按 `id` 的个体均值；
2. 把这些个体均值作为额外的回归变量放进模型。

这样做的好处：

- 可以同时识别短期（within）和长期（between）反应曲线。
- 保留了对新个体的预测能力（不像传统 FE 那样依赖训练集里的 `id`）。
- 与 Stata `mfxtsemipar2_cv.ado` 的行为一致。

> **注意**：通常不要把 `id` 再放入 `absorb`，否则可能与 `cvars` 产生完全共线性。

---

## 4. 交叉验证与结数选择

### 4.1 CV 折的生成

- 若未提供 `cvgroup`，则按 `id` 聚类生成 `nfold` 个折。
- 每个折包含若干完整个体，保证训练集和验证集在个体层面不重叠。

### 4.2 结数选择

对从 `minnk` 到 `maxnk` 的每个候选结数：

1. 在 `hf` 上生成样条基；
2. 按 `id × tl` 求和并合并到 `lf`；
3. 生成 `cvars`（个体均值）；
4. 用 CV 计算 RMSE，选择使 CV RMSE 最小的结数（或 `sopt = TRUE` 时选第一个局部最小）。

### 4.3 `iabsorb` 的作用

如果模型中通过 `absorb` 吸收了时间 FE 等粗层级固定效应，CV 时训练折估计出的 FE 需要加到验证折的预测中。`iabsorb = TRUE` 会：

- 在训练折上用 `fixest::feols` 估计模型；
- 用 `predict(fit, newdata = val, fixef = TRUE)` 提取验证折上的 FE；
- 从验证折残差中扣除这部分 FE。

若验证折出现训练折中没有的 FE 层级，则回退为 0（与 Stata 行为一致）。

---

## 5. 稳健纠偏（Robust Bias Correction, RBC）

R 实现继承了 `mfxtsemipar_bc.R` 的纠偏框架，默认 `bias_correct = TRUE`。

### 5.1 思想

用一组更高阶的样条基 $\tilde{B}(u)$（通常 `bc_degree = degree + 1`）去估计并修正主样条基的偏差。

### 5.2 两种估计方式

- **`bc_est = "stacked"`（默认）**：把主方程和纠偏方程堆叠成系统方程，分别估计 $\theta$ 和纠偏系数 $\gamma$。
- **`bc_est = "joint"`**：用正交化后的纠偏基 $\tilde{B} - BP$ 与主基一起做一个联合回归，其中 $P = (B'WB)^{-1}B'W\tilde{B}$。

### 5.3 两种纠偏预测

- **bc1**：直接用高阶基预测，$\hat{g}_{BC1}(u) = \tilde{B}(u)'\gamma$。
- **bc2**（默认）：主基预测 + 正交化纠偏，
  $$
  \hat{g}_{BC2}(u) = B(u)'\theta + \left[\tilde{B}(u) - B(u)P\right]'\gamma
  $$

在 `bc_est = "joint"` 时，bc1 与 bc2 数值相同。

### 5.4 短期与长期的纠偏曲线

- 短期：使用 within 系数 $\theta$ 和纠偏系数 $\gamma$。
- 长期：在短期基础上再加上 between 系数 $\delta$：
  $$
  \hat{g}_{LR}(u) = B(u)'(\theta + \delta) + \left[\tilde{B}(u) - B(u)P\right]'\gamma
  $$

---

## 6. Uniform Confidence Bands (UCB)

当 `ucb = TRUE` 且 `bias_correct = TRUE` 时，函数会对纠偏后的短期和长期反应曲线计算 **uniform confidence band**。

### 6.1 计算原理

UCB 基于 wild bootstrap 协方差矩阵与多元正态模拟的 sup-t 临界值：

1. 在默认的 200 个等距网格点（或用户提供的 `ucb_grid`）上计算曲线的点态估计与标准误；
2. 从估计系数的（wild bootstrap 或解析）协方差矩阵中抽取多元正态随机向量；
3. 对每次模拟计算网格上绝对 t 统计量的最大值，得到 sup-t 分布的样本；
4. 取该样本的 `ucb_level` 分位数作为临界值 `crit_sr` 与 `crit_lr`；
5. 将临界值乘以每个高频观测的纠偏方差标准误，得到
   $$\text{UB}(u) = \hat{g}(u) + \text{crit} \times \widehat{SE}_{BC}(u),$$
   $$\text{LB}(u) = \hat{g}(u) - \text{crit} \times \widehat{SE}_{BC}(u).$$

> 注意：严格来说，sup-t 临界值对网格点上的曲线同时有效；把它应用到所有高频观测是一种实用近似。

### 6.2 参数

| 参数 | 说明 |
|------|------|
| `ucb` | 是否计算 UCB，默认 `FALSE` |
| `ucb_level` | 覆盖水平，默认 `0.95` |
| `ucb_sim_reps` | sup-t 模拟次数，默认 `2000` |
| `ucb_grid` | 自定义评估网格，默认 `NULL`（在 `uvar` 范围内生成 200 个点） |

### 6.3 输出

- `res$fitted` 中增加 `<gen>_sr_lb`、`<gen>_sr_ub`、`<gen>_lr_lb`、`<gen>_lr_ub`。
- `res$ucb` 包含评估网格、临界值 `crit_sr`/`crit_lr`、覆盖水平与模拟次数。
- `predict(res, newdata = grid, uvar = "x3", ucb = TRUE)` 会返回 `lb_sr`、`ub_sr`、`lb_lr`、`ub_lr`。

---

## 7. R 函数用法

### 7.1 基本调用

```r
library(data.table)
source("R/mfxtsemipar2_cv.R")

res <- mfxtsemipar2_cv(
  hf = hf,          # 高频数据
  lf = lf,          # 低频数据
  y = "y",          # 被解释变量
  x = "x2",         # 低频控制变量
  uvar = "x3",      # 半参数变量（高频）
  id = "id",        # 面板标识
  tl = "t",         # 时间层级变量
  gen = "fitted",   # 输出拟合值前缀
  type = "poly",
  degree = 1,
  center = 0,
  absorb = ~ t,     # 吸收时间 FE，不把 id 放进去
  partialout = "all",
  maxnk = 5,
  minnk = 2,
  nfold = 5,
  seed = 123,
  iabsorb = TRUE
)
```

### 7.2 主要参数

| 参数 | 说明 |
|------|------|
| `hf` / `lf` | 高频 / 低频数据框或 data.table |
| `y` | 被解释变量名（在 `lf` 中） |
| `x` | 低频控制变量名向量 |
| `uvar` | 半参数变量名（在 `hf` 中） |
| `id` / `tl` | 面板标识和时间层级变量 |
| `gen` | 输出拟合值变量名前缀 |
| `hfcov` | 高频控制变量名向量（自动平均到 LF） |
| `type` | 样条类型：`"poly"`、`"bs"`、`"ms"`、`"is"`、`"ibs"` |
| `maxnk` / `minnk` | 最大 / 最小结数 |
| `absorb` | 固定效应公式或字符向量；**通常不放 `id`** |
| `partialout` | `NULL`、`"all"` 或变量名向量 |
| `iabsorb` | 是否在 CV 中传播训练折 FE 到验证折 |
| `bias_correct` | 是否做稳健纠偏 |
| `bc_type` | `"bc2"`（默认）或 `"bc1"` |
| `bc_est` | `"stacked"`（默认）或 `"joint"` |
| `brep` | wild bootstrap 次数，默认 0 |
| `ucb` | 是否计算 uniform confidence band |
| `ucb_level` | UCB 覆盖水平，默认 0.95 |
| `ucb_sim_reps` | sup-t 模拟次数，默认 2000 |
| `ucb_grid` | UCB 评估网格，默认 `NULL` |

### 7.3 返回值

`res` 是一个 `mfxtsemipar2_cv` 对象，主要字段：

- `res$nknots`：选定的结数。
- `res$knots`：结的位置。
- `res$cv_mse`：各结数下的 CV RMSE。
- `res$rmse`：最终拟合的样本内 RMSE。
- `res$coef` / `res$vcov`：最终模型系数和方差矩阵。
- `res$coef_bc` / `res$vcov_bc`：纠偏后的系数和方差矩阵（若 `bias_correct = TRUE`）。
- `res$fitted`：高频级别的拟合值 data.table，包含：
  - `<gen>_sr`、`<gen>_lr`：纠偏后的短期、长期拟合值
  - `<gen>_sr_se`、`<gen>_lr_se`：原始标准误
  - `<gen>_raw_sr`、`<gen>_raw_lr`：未纠偏的原始拟合值
  - 若 `bias_correct = TRUE`：`<gen>_bc1_sr`、`<gen>_bc1_lr`、`<gen>_bc_sr_se`、`<gen>_bc_lr_se` 等
  - 若 `ucb = TRUE`：`<gen>_sr_lb`、`<gen>_sr_ub`、`<gen>_lr_lb`、`<gen>_lr_ub`
- `res$ucb`：UCB 网格、临界值与模拟信息（若 `ucb = TRUE`）。
- `res$estimation`：完整的 `fixest` 估计对象。

### 7.4 预测新网格

```r
grid <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 100))
pred <- predict(res, newdata = grid, uvar = "x3")

# pred 包含：u, g_sr, g_lr, se_sr, se_lr
head(pred)

# 同时输出 UCB
pred_ucb <- predict(res, newdata = grid, uvar = "x3", ucb = TRUE)
# pred_ucb 额外包含：lb_sr, ub_sr, lb_lr, ub_lr
```

---

## 8. 输出解读示例

```r
print(res)
```

典型输出：

```
Mixed-frequency semiparametric regression with short-/long-run curves and cross-validation
  Selected knots: 3
  Boundary knots: -4.0594, 3.3891
  Knot locations: -0.6834, -0.0196, 0.6483
  Minimum CV RMSE: 4.578577
  Final fit RMSE: 4.152158
  Simple optimal knots: 3
  iabsorb: TRUE
  Bias correction: bc2 (degree 2)
  BC estimation: stacked
```

- `Selected knots`：CV 选出的结数。
- `Simple optimal knots`：CV 曲线的第一个局部最小值（`sopt = TRUE` 时使用）。
- `Final fit RMSE`：最终模型在全部样本上的 RMSE。
- `Bias correction`：当前使用的纠偏类型和阶数。

---

## 9. 与 Stata 原命令的对应

| Stata 选项 | R 参数 |
|------------|--------|
| `uvar()` | `uvar` |
| `id()` | `id` |
| `tl()` | `tl` |
| `gen()` | `gen` |
| `absorb()` | `absorb` |
| `partialout` / `partialout1()` | `partialout` |
| `iabsorb` | `iabsorb` |
| `type()` | `type` |
| `maxnk()` / `minnk()` | `maxnk` / `minnk` |
| `sopt` | `sopt` |
| `brep()` | `brep` |
| `cluster()` | `cluster` |
| `weights` | `weights` |

主要增强：

- R 版默认带稳健纠偏（Stata 原命令没有）。
- R 版显式输出短期和长期两条反应曲线。
- `partialout` 统一了 Stata 的 `partialout` 和 `partialout1(...)`。

---

## 10. 注意事项

1. **不要把 `id` 放进 `absorb`**：个体异质性已由 `cvars` 处理，再放 `id` 会导致共线性。
2. **`iabsorb` 的适用场景**：当 `absorb` 包含时间 FE 等粗层级 FE 时打开；若 FE 就是 `id` 本身，传播 FE 到验证折没有意义。
3. **低频数据 `lf` 必须唯一**：每个 `id × tl` 只能有一行。
4. **样条类型依赖**：`"bs"`、`"ms"`、`"is"`、`"ibs"` 需要 `splines2` 包。
5. **纠偏方程需要更多数据**：当样本较小时，高阶纠偏基可能导致共线性，可适当降低 `bc_degree` 或 `bc_nknots`。
6. **UCB 计算量**：`ucb = TRUE` 会进行大量多元正态模拟，全样本高频观测多时建议先用较小的 `ucb_sim_reps` 测试。
