# 📘 OT-tutorial: Optimal Transport 入门与实践教程

本项目旨在以简明通俗的方式系统介绍 Optimal Transport（最优传输）理论及其在机器学习中的应用。适合对数学基础有一定了解、希望深入掌握 OT 理论与实战的学习者。

---

项目导航  
- 项目主页：https://github.com/Day333/OT-tutorial.git
- 教程作者：Batch128  
- 参考资料：  
  - coderlemon 博客：https://coderlemon17.github.io/posts/2022/07-16-ot/  
  - YouTube 教学：https://www.youtube.com/watch?v=kjOBJP7gglw&list=PLJ6garKOlK2qKVhRm6UwvcQ46wK-ciHbl  
  - Notes on Optimal Transport: https://michielstock.github.io/posts/2017/2017-11-5-OptimalTransport/#notes_on_optimal_transport
  - 文献1: Peyré, Gabriel, and Marco Cuturi. "Computational optimal transport: With applications to data science." Foundations and Trends® in Machine Learning 11.5-6 (2019): 355-607.
  - 文献2: Séjourné, Thibault, Gabriel Peyré, and François-Xavier Vialard. "Unbalanced optimal transport, from theory to numerics." Handbook of Numerical Analysis 24 (2023): 407-471.

---

## 📚 教程大纲

0. 引言：大海的聚会 🌊

1. 第一章：OT 起源与基础直觉  
  1.1 OT 是什么？为什么重要？  
  1.2 从 Monge 问题到 Kantorovich 放松  
  1.3 Wasserstein 距离的几何意义与直观理解  
  📓 Notebook：用 `scipy.optimize` 实现基础 OT 问题

2. 第二章：数学基础与理论推导  
  2.1 Monge/Kantorovich 两种形式  
  2.2 Wasserstein 距离的性质与距离空间  
  2.3 对偶形式（Dual formulation）与强对偶  
  📓 Notebook：构造 cost matrix 计算 OT 距离

3. 第三章：Sinkhorn 算法与快速求解  
  3.1 熵正则化的动机与效果  
  3.2 Sinkhorn-Knopp 算法推导  
  3.3 数值稳定性优化（log-domain 实现）  
  📓 Notebook：使用 `POT` / `geomloss` 实现 Sinkhorn OT

4. 第四章：OT 在机器学习中的应用  
  4.1 生成模型中的 OT（WGAN, WAE 等）  
  4.2 OT 用于 domain adaptation（CORAL, JDOT）  
  4.3 Word Mover's Distance 在 NLP 中的使用  
  📓 Notebook：图像迁移 / 句子匹配实战演示

5. 第五章：高维分布与连续映射上的 OT  
  5.1 Brenier map 与单射映射理论  
  5.2 OT PDE 形式：Monge-Ampère 方程简介  
  5.3 高维近似方法：Sliced OT、Projection Robust OT  
  📓 Notebook：图像风格迁移与 Sliced OT 实践

6. 第六章：图上的 OT 与离散结构  
  6.1 图结构下的最优匹配问题  
  6.2 Gromov-Wasserstein 距离原理与直觉  
  📓 Notebook：演示图结构迁移与匹配任务

7. 第七章：实践进阶与扩展  
  7.1 OT Python 库快速指南（POT, GeomLoss, OTT-JAX）  
  7.2 实战案例：图像匹配 / 域适配 / 分布对齐  
  7.3 拓展阅读与研究前沿  
  📓 Notebook：OTT-JAX 实战对比多个算法性能

---

🗂 项目结构
```
OT-tutorial/
├── notebooks/ # 每章 Jupyter Notebook
├── docs/ # 教学文档与公式推导
├── resources/ # 推荐论文、图书与视频
├── slides/ # 教学演示幻灯片（可选）
├── requirements.txt # Python 依赖
├── README.md # 项目首页
└── LICENSE # 开源许可证
```
---

## 🔧 快速开始

```bash
git clone https://github.com/Day333/OT-tutorial.git
cd OT-tutorial
pip install -r requirements.txt
jupyter notebook