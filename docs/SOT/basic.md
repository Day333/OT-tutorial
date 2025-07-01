# 1 OT是什么？

OT是Optimal Transport的缩写，是一种数学工具，用于解决优化问题。 Optimal Transport（OT）是研究如何以最小成本将一个概率分布“搬运”成另一个概率分布的问题。

## 形式化
给定两个概率分布 $\mu$ 和 $\nu$，以及一个代价函数 $c(x, y)$，OT 问题是寻找一个搬运方案（联结 $\gamma$），使得总搬运成本最小：

$$
\inf_{\gamma \in \Gamma_{\mu, \nu}} \int c(x, y)\, d\gamma(x, y)
$$

其中：

- $\Gamma_{\mu, \nu}$ 是所有联结 $\gamma$ 的集合（即边缘分布为 $\mu, \nu$ 的联合分布）；

- $c(x, y)$ 是从 $x$ 搬到 $y$ 的代价，例如欧氏距离的平方。

# 2 Monge 问题是啥？
Monge 是指法国数学家 Gaspard Monge，也是“最优传输问题”的最早提出者。他在 1781 年发表的一篇论文中首次提出了经典的“搬运沙子”的问题，因此这个问题后来就被称为 Monge 问题（Monge's Optimal Transport Problem）。

**Monge 问题（Monge OT）** 是最优传输问题的最初形式：

给定两个分布 $\mu$ 和 $\nu$，找到一个确定性搬运函数 $T$，把 $\mu$ 搬到 $\nu$，使总搬运成本最小。

🧱 数学表达：
给定：
源分布 $\mu$（比如沙子在哪里）；
目标分布 $\nu$（坑在哪里）；
成本函数 $c(x, T(x))$（搬运 $x$ 到 $T(x)$ 的代价）；

Monge 的问题是：

$$
\inf_{T : T_\# \mu = \nu} \int c(x, T(x))\, d\mu(x)
$$

其中：

- $T_{\#} \mu = \nu$ 表示 $T$ 是将 $\mu$ 推送成 $\nu$ 的映射（pushforward measure）；

- $c(x, T(x))$ 是搬运成本；

- $T$ 是传输映射（transport map）。

Monge 问题的限制：
Monge 要求：
- 每个 $x$ 只能搬到一个 $T(x)$，也就是说搬运是确定性的。

这导致有些时候 Monge 的传输映射 $T$ 不存在，比如：

- 如果 $\mu = \delta_0$（全部集中在 $x=0$），而 $\nu = \frac{1}{2} \delta_{-1} + \frac{1}{2} \delta_1$（目标是一分为二），就找不到 $T$，因为一个点不能一分为二搬到两个地方。

# 3 Kantorovich 问题是什么？

Kantorovich（列昂尼德·康托洛维奇） 是苏联著名数学家、经济学家，1940 年代提出了 Kantorovich 最优传输问题，即最早的 Monge 问题的“放松”版本。 Kantorovich 最优传输问题（KOT） 是对 Monge 问题的放宽（relaxation）：

不再要求一个确定的搬运映射 $T$，而是允许一个“搬运计划” $\gamma(x, y)$，表示从 $x$ 搬运到 $y$ 的比例。

设：

- 源分布：$\mu$；
- 目标分布：$\nu$；
- 成本函数：$c(x, y)$（从 $x$ 搬到 $y$ 的代价）；

Kantorovich 问题是：

$$
\inf_{\gamma \in \Gamma_{\mu, \nu}} \int_{\mathbb{R}^d \times \mathbb{R}^d} c(x, y)\, d\gamma(x, y)
$$

其中：

- $\gamma$ 是一个联合分布（也叫 **联结 coupling**），满足边缘为 $\mu$ 和 $\nu$：
  
  $$
  \gamma(A \times \mathbb{R}^d) = \mu(A), \quad \gamma(\mathbb{R}^d \times B) = \nu(B)
  $$

- $\Gamma_{\mu, \nu}$ 表示所有这样的 $\gamma$ 的集合。


# 4 Wasserstein 距离（Wasserstein Distance） 的 连续形式 与 离散形式

Wasserstein 距离来源于 Kantorovich 最优传输问题（KOT），是该问题在特定成本函数下的最小值。

## 1. Wasserstein 距离概述

给定两个概率分布 $\mu$ 和 $\nu$，以及一个搬运成本函数 $c(x, y) = \|x - y\|^p$，$p$ 阶 Wasserstein 距离定义为：

$$
W_p(\mu, \nu) := \left( \inf_{\gamma \in \Gamma_{\mu, \nu}} \int_{\mathbb{R}^d \times \mathbb{R}^d} \|x - y\|^p \, d\gamma(x, y) \right)^{1/p}
$$

其中：
- $\Gamma_{\mu, \nu}$ 表示所有边缘为 $\mu$ 和 $\nu$ 的联合分布（联结 coupling）；
- $p \geq 1$；
- 最常用的是 $p=1$ 和 $p=2$，分别称为 $W_1$ 和 $W_2$ 距离。

---

## 🌊 2. 连续形式（Continuous Form）

### 💡 定义：

若 $\mu, \nu$ 是 $\mathbb{R}^d$ 上的连续概率测度，满足 $\int \|x\|^p \, d\mu(x) < \infty$（“总搬运代价的期望值”），则 $p$ 阶 Wasserstein 距离定义为：

$$
W_p^p(\mu, \nu) = \inf_{\gamma \in \Gamma_{\mu, \nu}} \int_{\mathbb{R}^d \times \mathbb{R}^d} \|x - y\|^p \, d\gamma(x, y)
$$

然后取 $p$ 次方根得到 $W_p(\mu, \nu)$。

---

### 📌 常见情形（欧几里得空间）：

- **$W_1$ 距离**：
  $$
  W_1(\mu, \nu) = \inf_{\gamma \in \Gamma_{\mu, \nu}} \int \|x - y\| \, d\gamma(x, y)
  $$

- **$W_2$ 距离**：
  $$
  W_2^2(\mu, \nu) = \inf_{\gamma \in \Gamma_{\mu, \nu}} \int \|x - y\|^2 \, d\gamma(x, y)
  $$

---

## 📦 3. 离散形式（Discrete Form）

当 $\mu$ 和 $\nu$ 是有限支持的离散分布时，Wasserstein 距离转化为线性规划问题。

### 🧮 离散输入：

设：
- $\mu = \sum_{i=1}^m p_i \delta_{x_i}$；
- $\nu = \sum_{j=1}^n q_j \delta_{y_j}$；

其中：
- $x_i, y_j \in \mathbb{R}^d$ 是支持点；
- $p = (p_1, ..., p_m)$，$q = (q_1, ..., q_n)$ 是对应的概率权重。

定义代价矩阵 $C \in \mathbb{R}^{m \times n}$，其中：
$$
C_{ij} = \|x_i - y_j\|^p
$$

---

### 🧩 Wasserstein 距离变为线性规划：

引入耦合矩阵 $P \in \mathbb{R}_+^{m \times n}$（$\gamma$ 的离散版本），满足边缘约束：

$$
P \mathbf{1}_n = p, \quad P^\top \mathbf{1}_m = q
$$

则 Wasserstein 距离为：

$$
W_p^p(\mu, \nu) = \min_{P \in \mathbb{R}_+^{m \times n}} \langle C, P \rangle 
\quad \text{s.t.} \quad P \mathbf{1}_n = p, \quad P^\top \mathbf{1}_m = q
$$

即最小化：

$$
\sum_{i=1}^m \sum_{j=1}^n C_{ij} P_{ij}
$$

这就是经典的 **运输问题（transportation problem）线性规划模型**。


注释：

| From ↓ / To →  | \$y\_1\$ | \$y\_2\$ | \$y\_3\$ | 行和（= \$p\_i\$） |
| -------------- | -------- | -------- | -------- | -------------- |
| \$x\_1\$       | 0.2      | 0.1      | 0.0      | 0.3            |
| \$x\_2\$       | 0.3      | 0.0      | 0.1      | 0.4            |
| \$x\_3\$       | 0.0      | 0.1      | 0.2      | 0.3            |
| 列和（= \$q\_j\$） | 0.5      | 0.2      | 0.3      |                |

$P \mathbf{1}_n = p$ 和 $P^\top \mathbf{1}_m = q$ 表示：搬出去的不能多于源分布，搬进来的不能少于目标分布，是离散最优传输中的“质量守恒”约束。