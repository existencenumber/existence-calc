
---

#### 2. JOSS 论文 (paper.md)

这是 JOSS 的核心提交文件，用于介绍你的软件。我已经帮你填好了相关的项目地址和元数据，直接复制使用即可。

```markdown
---
title: 'Collapse Monster: A Universal Divergent Series Regularization Engine Based on Existence Number Theory'
authors:
  - name: Wanpeng Xu
    orcid: 0009-0008-6677-4582
    affiliation: "Independent Researcher"
tags:
  - Python
  - divergent series
  - regularization
  - symbolic computation
  - mathematical physics
date: 2026-05-27
bibliography: paper.bib
---

# Summary

`Collapse Monster` is a Python library and web application that automatically regularizes divergent infinite series by mapping them through the nine-domain dual graph of Existence Number Theory [@xu2026]. The program identifies the asymptotic growth type of a given series, searches for the shortest path (at most 3 steps) to a convergent domain, and applies the corresponding analytic continuation (zeta regularization, Abel summation, Euler summation, Borel summation, or Jones polynomial regularization for super-exponential series) to produce a finite ``collapsed'' value.

# Statement of Need

Divergent series appear ubiquitously in theoretical physics (perturbative quantum field theory, string theory, statistical mechanics), but traditional mathematical software treats them as invalid operations. Existing regularization methods (Ramanujan, Abel, Borel, etc.) are implemented as isolated procedures in computer algebra systems, with no unifying principle. `Collapse Monster` fills this gap: it operationalizes the Dual Diameter Theorem of Existence Number Theory, which guarantees that any divergent dynamic number can be mapped to a convergent representation within at most three domain transformations. The tool provides a single, automated interface for all classical summation methods and extends regularization to super-exponential divergences via braided domain (multi-layer Borel) techniques.

# Features

- Accepts standard summation notation (e.g., `Sum(n**2,(n,1,oo))`)
- Automatically classifies divergence type and identifies the original computational domain
- Implements BFS path-finding on the nine-domain dual graph
- Supports zeta, Dirichlet eta, Abel, Euler, Borel, and multi-layer Borel regularizations
- Handles super-exponential series (e.g., $n^n$) up to depth 4
- Provides a Flask web interface for interactive use
- Returns the mapping path, number of steps, and the collapsed value

# Example Usage

Input: Sum(n,(n,1,oo))       → Collapsed: -1/12 (zeta(-1))
Input: Sum(2**n,(n,0,oo))    → Collapsed: -1 (Abel)
Input: Sum(factorial(n),(n,0,oo)) → Collapsed: 0.596347 (Borel)
Input: Sum((-1)**(n+1)/n,(n,1,oo)) → Collapsed: ln(2) (Dirichlet eta)
Input: Sum(n**n,(n,1,oo))    → Collapsed: finite (braided domain)

# References

@book{xu2026,
  title={Existence Number Theory: A Mathematical Language Centered on Process},
  author={Xu, Wanpeng},
  year={2026},
  publisher={Preprint}
}
