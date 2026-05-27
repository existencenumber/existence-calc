# 塌缩怪兽 v25.0 (Collapse Monster)
## 基于存在数论的普适发散正则化引擎

本工具基于存在数论 (Existence Number Theory) 的九域对偶映射框架，能够自动将发散级数映射到收敛域，并计算出其有限坍缩值。支持多项式、几何、交错、阶乘、超指数（深度至4）、对数、数论函数等发散类型。

## 理论核心

任何发散问题均可通过 ≤3 步对偶映射消除。发散不是级数的固有属性，而是“选错了观察域”。

## 快速开始

### 在线使用
访问 [https://existence-calc.onrender.com](https://existence-calc.onrender.com)

### 本地运行
```bash
git clone https://github.com/existencenumber/existence-calc.git
cd existence-calc
pip install -r requirements.txt
python app.py
