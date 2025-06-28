# 1 Optimal Transport

## 1.1 最优传输问题

# 1.1 最优传输问题（The Optimal Transport Problem）

## 🚜 Monge 问题的起源

在 1781 年的备忘录中，Monge 提出了一个经典问题：

> 如何将一堆沙子运送到一个同样体积的坑中，且使运输成本最小？

我们可以将这个问题建模为概率测度之间的转换问题。令：
- $\mu$：表示沙堆的位置分布；
- $\nu$：表示坑的位置分布。

由于体积守恒，不妨设 $\mu$, $\nu$ 都是定义在 $\mathbb{R}^d$ 上的概率测度（总质量为 1）。我们也可以引入随机变量：
- $X \sim \mu$
- $Y \sim \nu$

### 1.1.1 The Monge and Kantorovich problems

回到我们的沙子类比，运输这堆沙子意味着找到一个（可测的）函数，称为传输映射（transport map）$T : \mathbb{R}^d \rightarrow \mathbb{R}^d$，它指明位于 $x \in \mathbb{R}^d$ 的沙子应被移动到 $T(x) \in \mathbb{R}^d$。

为了使这个传输映射真正完成任务（即填满那个坑），我们需要确保 $X \sim \mu$ 时有 $T(X) \sim \nu$。

我们说 $T$ 将 $\mu$ 推送为 $\nu$，或者说 $\nu$ 是 $\mu$ 通过 $T$ 得到的推送测度，并记作 $T_\# \mu = \nu$。这是我们的约束条件。

现在来考虑我们的目标函数，回忆 Monge 的问题是要最小化运输沙子的成本。

衡量成本的方法有很多（如体力、燃料消耗等），为了简化表述，我们用沙子移动的欧几里得距离来衡量。

位于位置 $x$ 的沙子移动的距离为 $\|T(x) - x\|$。

因此，平均搬运代价为：
$$
\int \|T(x) - x\| \mu(dx)。
$$

因此，Monge 形式的最优传输问题是要在满足 $T_\# \mu = \nu$ 的条件下，最小化上述目标：
$$
\inf_{T : T_\# \mu = \nu} \int \|T(x) - x\| \mu(dx)。
$$

注意，有许多种方式可以选择传输成本。

通常，我们考虑一般形式的成本函数 $c(X, T(X))$，其中 $c(x, y)$ 衡量将 $x \in \mathbb{R}^d$ 传输到 $y \in \mathbb{R}^d$ 的成本。

在这个更一般的框架中，我们甚至可以允许 $X$ 和 $Y$ 定义在两个不同的空间上，不一定是 $\mathbb{R}^d$。

在这些讲义中，我们主要关注 $c(x, y) = \|x - T(x)\|$ 或 $c(x, y) = \|x - T(x)\|^2$ 的情形，它们导出 Wasserstein 距离。

空间 $\mathbb{R}^d$ 也可以替换为更复杂的空间，如黎曼流形，但这一般超出了本讲义的讨论范围（除了第 5.6 节以外）。

虽然 Monge 问题很容易表述，但我们需要提出几个问题：

- 总是存在这样的传输映射吗？

- 如果存在一个最小化器，它是否唯一？如何刻画它？注意，我们的约束不是凸的，这使得回答这个问题相当困难。

我们可以通过一个简单的例子来看，Monge 的传输映射 $T$ 并不总是存在。

假设我们只考虑一维的情形，也就是说 $d = 1$，一切都在一条数轴上。

我们把原始的“沙堆”集中在数轴上的一个点——$x = 0$，也就是说所有的质量都集中在这一个位置，用数学表达就是 $\mu = \delta_0$。

而目标“沙坑”并不是在某一个点上，而是被分成了两半：一半在 $x = -1$，另一半在 $x = 1$。这就是目标分布 $\nu = \frac{1}{2} \delta_{-1} + \frac{1}{2} \delta_1$ 的意思。它表示：“沙子最后要平均分布在 $x = -1$ 和 $x = 1$ 两点”。

Monge 要求：一个点只能搬到一个地方
在 Monge 的最初设定中，每一点上的沙子只能被搬到一个特定的位置，也就是说：

如果你在 $x = 0$，你必须决定：把沙子全都搬到 $T(0) = -1$，或者 $T(0) = 1$，不能一半搬左边、一半搬右边。

数学上说，这种“只能搬到一个点”的规则就是所谓的确定性映射 $T$。

因此，在这种情形下，Monge 的方式是无法完成目标的，因为你无法把一个点的沙子“一分为二”。

也就是说，没有任何符合要求的函数 $T$ 可以把原始的 $\mu$ 映射成目标的 $\nu$。

直觉上，我们希望：
$$
T(0) =
\begin{cases}
-1 & \text{概率 } \frac{1}{2} \\
1 & \text{概率 } \frac{1}{2}
\end{cases}
$$

并且对所有 $x \ne 0$ 有 $T(x) = x$。

这样的 $T$ 不是一个函数，而是一个马尔可夫核（Markov kernel）：它给 $\mathbb{R}$ 中的每个点 $x$ 分配一个概率分布。

第二个问题在近两个世纪内一直没有令人满意的答案，直到苏联数学家 Leonid Kantorovich 在其一篇突破性的两页论文中提出了这个问题的松弛形式，该形式正好允许马尔可夫核，如上面的例子所示。

Kantorovich 放宽了限制：不再要求你一定得找到一个确定性的搬运路径（即函数 $T$），而是允许你从 $x$ 出发，以一定比例搬到多个 $y$ 上。

这种更宽松的“搬沙计划”叫做一个 联结（coupling），记作 $\gamma(x, y)$，它满足两个要求：

- 从整体看，有多少沙子从某个 $x$ 出发，加起来就是 $\mu(x)$；

- 搬到每个 $y$ 上的总沙子量，加起来就是 $\nu(y)$。

设 $\mu, \nu$ 是 $\mathbb{R}^d$ 上的两个概率测度，$\gamma$ 是这两个分布之间的联结，也就是说，是 $\mathbb{R}^d \times \mathbb{R}^d$ 上的联合分布，其第一边缘是 $\mu$，第二边缘是 $\nu$：对任意 Borel 集 $A \subset \mathbb{R}^d$ 有
$$
\gamma(A \times \mathbb{R}^d) = \mu(A) \quad \text{and} \quad \gamma(\mathbb{R}^d \times A) = \nu(A)。
$$

$\gamma$ 是一个“搬运方案”，它不能凭空创造或丢失沙子；对所有从 $x \in A$ 出发、搬到任意 $y$ 的质量加总起来，必须等于 $\mu$ 在 $A$ 上的总质量，搬到 $y$ 的总量等于 $\nu$。

“联结”这个术语的含义是，虽然 $X \sim \mu$ 与 $Y \sim \nu$ 本来是彼此无关的随机变量，但联结让它们生活在同一个概率空间中，并描述它们之间的概率依赖关系。

在这些讲义中，我们用 $\Gamma_{\mu, \nu}$ 表示 $\mu$ 和 $\nu$ 的所有联结的集合。

令 $c : \mathbb{R}^d \times \mathbb{R}^d \to [0, \infty)$ 是一个可测的成本函数。

**Kantorovich 形式的最优传输问题**如下：
$$
\inf_{\gamma \in \Gamma_{\mu, \nu}} \int c(x, y) \, \gamma(dx, dy) \tag{KOT}
$$

## 1.1.2 Couplings（耦合）

为了更好地理解 Kantorovich 问题，研究耦合集合 $\Gamma_{\mu,\nu}$ 会很有帮助。在数学上，耦合是一个联合分布 $\gamma$，它的边缘分布分别是 $\mu$ 和 $\nu$。

或许最简单的一种耦合是独立耦合（independent coupling）$\gamma = \mu \otimes \nu$，即假设 $X \sim \mu$ 与 $Y \sim \nu$ 是独立的：对于任意 Borel 集 $A, B \subseteq \mathbb{R}^d$，
$$
\gamma(A \times B) = \mu(A) \cdot \nu(B).
$$

下述命题总结了关于 $\Gamma_{\mu,\nu}$ 的一些初步事实。

**命题 1.1**：设 $\mu, \nu$ 是 $\mathbb{R}^d$ 上的两个概率测度。其耦合集合 $\Gamma_{\mu,\nu}$（即其耦合的集合）在弱收敛拓扑下是非空、凸且紧的。

注意：非空：至少存在一个耦合（比如独立耦合）。“凸”问题容易求最小值。“紧” = “不会跑到无穷远 + 总是有极小值点”。

**证明**：

- **非空**性由独立耦合 $\mu \otimes \nu$ 的存在性保证，因此 $\Gamma_{\mu,\nu} \neq \emptyset$。
  
- **凸性**：设 $\gamma_0, \gamma_1 \in \Gamma_{\mu,\nu}$，对任意 $\lambda \in (0,1)$，定义 $\gamma_\lambda = (1 - \lambda)\gamma_0 + \lambda\gamma_1$。对任意 Borel 集 $A \subseteq \mathbb{R}^d$，
  $$
  \gamma_\lambda(A \times \mathbb{R}^d) = (1 - \lambda)\mu(A) + \lambda \mu(A) = \mu(A),
  $$
  同理 $\gamma_\lambda(\mathbb{R}^d \times A) = \nu(A)$，因此 $\gamma_\lambda \in \Gamma_{\mu,\nu}$，说明 $\Gamma_{\mu,\nu}$ 是凸的。

- **紧性**：根据 Prokhorov 定理，如果一个测度集是紧的，只需证明其**紧致性**与**闭性**。由于 $\mu$ 和 $\nu$ 是概率测度，对于任意 $\varepsilon > 0$，存在紧集 $K \subset \mathbb{R}^d$，使得 $\mu(K^c) + \nu(K^c) < \varepsilon$，则 $K \times K$ 也是紧的，且对任意 $\gamma \in \Gamma_{\mu,\nu}$，
  $$
  \gamma((K \times K)^c) \leq \gamma(\mathbb{R}^d \times K^c) + \gamma(K^c \times \mathbb{R}^d) = \mu(K^c) + \nu(K^c) < \varepsilon,
  $$
  说明 $\Gamma_{\mu,\nu}$ 是紧的。

  又因为对所有有界连续函数 $f$ 有：
  $$
  \int f(x)\, \gamma(dx,dy) = \int f\, d\mu, \quad \int f(y)\, \gamma(dx,dy) = \int f\, d\nu,
  $$
  根据弱收敛定义 $\Gamma_{\mu,\nu}$ 是闭的。

因此，根据 Prokhorov 定理，$\Gamma_{\mu,\nu}$ 是紧的。$\square$

---

**图 1.1**：（左）两个高斯混合分布的独立耦合；（右）$X \sim \mathcal{N}(0,1)$ 与 $Y \sim \chi^2_1$ 的确定性耦合（$Y = X^2$）。

---

耦合 $\gamma \in \Gamma_{\mu,\nu}$ 捕捉了 $X \sim \mu$ 和 $Y \sim \nu$ 两个随机变量之间的依赖关系。如上所述，一种选择是认为 $X$ 与 $Y$ 独立，这对应独立耦合。图 1.1（左）展示了两个高斯混合分布之间的独立耦合。

在独立耦合的另一端，考虑 $X \sim \mathcal{N}(0,1)$，$Y \sim \chi^2_1$，注意到 $Y$ 与 $X^2$ 同分布。于是我们可以构造确定性耦合 $Y = X^2$，写作：
$$
\gamma(dx, dy) = \mu(dx) \delta_{x^2}(dy).
$$
这个耦合在图 1.1（右）中展示，是一个退化耦合（即完全由 $X$ 决定 $Y$）。

为了进一步探索耦合，设 $X \sim \mathcal{N}(0,1)$ 和 $Y \sim \mathcal{N}(0,1)$，那么我们可以构造如下的联合分布：
$$
\begin{pmatrix} X \\ Y \end{pmatrix} \sim \mathcal{N} \left( \begin{pmatrix} 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix} \right), \quad \rho \in [-1,1],
$$
即相关系数为 $\rho$ 的二维高斯耦合。见图 1.2 展示。


## 1.1.3 离散最优传输（Discrete Optimal Transport）

当 $\mu$ 和 $\nu$ 是两个离散分布时，这种情形具有特别的实用意义。例如，$\mu,\nu$ 可以是定义在点云上的经验测度。

我们考虑如下情形：

$$
\mu = \sum_{i=1}^m p_i \delta_{x_i}, \quad \nu = \sum_{j=1}^n q_j \delta_{y_j},
$$

其中 $\delta_{x_i}$ 和 $\delta_{y_j}$ 表示在点 $x_i$ 和 $y_j$ 上的狄拉克测度（即把概率质量集中在某个点上）。

注释：$\delta_{x_i}$ 是狄拉克测度，表示“所有质量都集中在 $x_i$ 上”；所以这个 $\mu$ 是把概率 $p_i$ 分别分配到点 $x_i$ 上，构成离散分布。

在这种情况下，$X \sim \mu$ 和 $Y \sim \nu$ 的耦合 $\gamma$ 可以用一个非负矩阵 $P \in \mathbb{R}_+^{m \times n}$ 来刻画，其中 $P_{ij} = \gamma(X = x_i, Y = y_j)$。

耦合的边缘约束 $\gamma \in \Gamma_{\mu,\nu}$ 可以表示为：

- 对所有 $i \in [m]$，有：
  $$
  \sum_{j \in [n]} P_{ij} = p_i
  $$
- 对所有 $j \in [n]$，有：
  $$
  \sum_{i \in [m]} P_{ij} = q_j
  $$

引入 $\mathbf{1}_m, \mathbf{1}_n$ 分别表示长度为 $m$ 和 $n$ 的全 1 向量，则这些约束可以简洁地写为：
- $P \mathbf{1}_n = p$
- $P^\top \mathbf{1}_m = q$

其中：
- $p = (p_1, ..., p_m)^\top$
- $q = (q_1, ..., q_n)^\top$

表格中形象展示了上面公式所表达的意思：

|             | $Y = y_1$ | $Y = y_2$ | $Y = y_3$ | 行和 $\sum_j P_{ij}$ |
|-------------|-----------|-----------|-----------|-----------------------|
| $X = x_1$   | $P_{11}$  | $P_{12}$  | $P_{13}$  | $p_1$                |
| $X = x_2$   | $P_{21}$  | $P_{22}$  | $P_{23}$  | $p_2$                |
| 列和 $\sum_i P_{ij}$ | $q_1$     | $q_2$     | $q_3$     |                       |


与耦合类似，代价也可以由一个 $m \times n$ 的矩阵 $C$ 表示，其中：

$$
C_{ij} = c(x_i, y_j)
$$

于是，**Kantorovich 最优传输问题（KOT）** 等价于以下线性规划问题：

$$
\min_{P \in \mathbb{R}_+^{m \times n}} \sum_{i,j \in [n]} C_{ij} P_{ij} \quad \text{s.t.} \quad P \mathbf{1}_n = p,\quad P^\top \mathbf{1}_m = q
$$

也可以更紧凑地写为：

$$
\min_{P \in \mathbb{R}_+^{m \times n}} \langle C, P \rangle \quad \text{s.t.} \quad P \mathbf{1}_n = p,\quad P^\top \mathbf{1}_m = q
$$

其中：
$$
\langle C, P \rangle = \text{tr}(C^\top P)
$$
是矩阵上的 Frobenius 内积（即元素对应相乘后求和）。

特别地，当 $m = n$ 且所有权重 $p_i, q_j$ 都等于 $1/n$ 时，所有满足条件的耦合矩阵 $P$（乘上常数n）就是所谓的**双随机矩阵（doubly stochastic matrix）**集合，也称为 **Birkhoff 多面体**：

注释：双随机矩阵的每行每列都加起来等于 1，表示均匀质量的匹配；它们构成一个特殊的凸多面体。

$$
\text{Birk} := \left\{ \gamma \in \mathbb{R}_+^{n \times n} : \gamma \mathbf{1}_n = \mathbf{1}_n,\ \mathbf{1}_n^\top \gamma = \mathbf{1}_n^\top \right\} \tag{1.2}
$$

此时，Kantorovich 最优传输问题简化为：

$$
\min_{P \in \text{Birk}} \langle C, P \rangle \tag{1.3}
$$

Birkhoff 多面体的极点是**置换矩阵**（permutation matrices）：

- 它们是二值矩阵 $\pi \in \{0,1\}^{n \times n}$，每一行和每一列都恰好有一个元素为 1。

特别地，线性规划几何理论告诉我们：

> 在线性目标函数 $\langle C, P \rangle$ 下，其最优解通常出现在多面体的极点上。

因此，(1.3) 的解可以取为一个 $n \times n$ 的置换矩阵 $\pi$，它表示某种确定性的映射关系（transport plan）：

注释：置换矩阵表示“完美一一对应”的匹配方案。比如：第1行第3列为1，表示第1个点被送去第3个点。

- 它告诉我们每一个 $x_i$ 被送往哪个 $y_j$。

这也被称为**Monge 问题**的解（完全确定性的匹配）。

正如我们将在后续中看到的那样，某些耦合（特别是退化的）在最优传输的几何中扮演着关键角色。
