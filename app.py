"""
塌缩怪兽 v22.0 — 修复 JSON 序列化 + 增强帮助页
"""

import math, os, traceback, json
from collections import deque, defaultdict

import sympy as sp
from sympy import (oo, factorial, log, Symbol, Sum, Function,
                   zeta, pi, exp, sqrt)
import mpmath as mp

# ========== 全局符号 ==========
n_sym = Symbol('n', integer=True, positive=True)
x_sym = Symbol('x')

SAFE_LOCALS = {
    "Sum": Sum, "oo": oo, "factorial": factorial,
    "log": sp.log, "sin": sp.sin, "cos": sp.cos,
    "exp": sp.exp, "sqrt": sp.sqrt, "n": n_sym,
    "mobius": Function("mobius"),
    "fibonacci": Function("fibonacci"),
    "liouville": Function("liouville"),
    "eulerphi": Function("eulerphi"),
    "divisor_sigma": Function("divisor_sigma"),
    "mangoldt": Function("mangoldt"),
    "zeta": zeta, "pi": pi
}

# ========== 九域 ==========
DOMAINS = ["加法域","乘法域","积分域","微分域","谱域","泛函积分域","编织域","同伦域","范畴域"]
CONVERGENCE_DOMAINS = {"谱域","泛函积分域","编织域","同伦域","范畴域"}

DUAL_GRAPH = {
    "加法域": {"指数映射":"乘法域","黎曼和极限":"积分域","差商极限":"微分域"},
    "乘法域": {"对数映射":"加法域","对数导数":"积分域","梅林变换":"谱域"},
    "积分域": {"拉普拉斯变换":"谱域","泛函极限":"泛函积分域","积分的逆":"微分域"},
    "微分域": {"傅里叶变换":"谱域","微分的逆":"积分域"},
    "谱域": {"逆傅里叶变换":"微分域","逆拉普拉斯变换":"积分域","逆梅林变换":"乘法域"},
    "泛函积分域": {"泛函极限的逆":"积分域","二维拓扑":"编织域"},
    "编织域": {"辫子同伦":"同伦域"},
    "同伦域": {"态射范畴化":"范畴域"},
    "范畴域": {"恒等态射对应0":"加法域","恒等态射对应1":"乘法域"}
}

# ========== 动态数 ==========
class DynamicNumber:
    def __init__(self, expr, domain, start=1, history=None):
        self.expr = expr
        self.domain = domain
        self.start = int(start)
        self.history = history or []

    def evolve(self, new_expr, new_domain, mapping_name):
        return DynamicNumber(new_expr, new_domain, self.start,
                             self.history + [(mapping_name, new_domain)])

# ========== 函子 ==========
FUNCTOR_REGISTRY = {}
def register_transform(src, dst, name, func):
    FUNCTOR_REGISTRY[(src, dst)] = (name, func)

for src, mappings in DUAL_GRAPH.items():
    for name, dst in mappings.items():
        register_transform(src, dst, name,
                           lambda dn, n=name, d=dst: dn.evolve(dn.expr, d, n))

def exp_map(dn): return dn.evolve(exp(dn.expr), "乘法域", "指数映射")
def log_map(dn): return dn.evolve(log(dn.expr), "加法域", "对数映射")
def riemann_sum_limit(dn): return dn.evolve(dn.expr.subs(n_sym, x_sym), "积分域", "黎曼和极限")
register_transform("加法域", "乘法域", "指数映射", exp_map)
register_transform("乘法域", "加法域", "对数映射", log_map)
register_transform("加法域", "积分域", "黎曼和极限", riemann_sum_limit)

# ========== 寻路 ==========
class PathFinder:
    def __init__(self):
        self.graph = DUAL_GRAPH
        self.convergence = CONVERGENCE_DOMAINS
    def find_all_paths(self, start_domain, max_steps=3):
        visited = {start_domain}
        queue = deque([(start_domain, [])])
        results = []
        while queue:
            current, path = queue.popleft()
            if current in self.convergence and path:
                results.append(path)
            if len(path) >= max_steps: continue
            for mapping, neighbor in self.graph.get(current, {}).items():
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [(mapping, neighbor)]))
        return sorted(results, key=len)

# ========== 求值器 ==========
class Evaluator:
    def evaluate(self, dn):
        if dn.domain == "谱域":
            return self._eval_spectral(dn)
        return None

    def _eval_spectral(self, dn):
        expr = dn.expr
        start = dn.start

        if expr == sp.Integer(0) or expr == 0: return 0.0
        if expr == sp.Integer(1) or expr == 1:
            return float('inf') if start >= 1 else 1.0

        # 调和
        if sp.simplify(expr - 1/n_sym) == 0:
            return float('inf')

        # 交错级数 —— 最优先处理，防止被几何基捕获
        if self._is_alternating(expr):
            core, parity = self._extract_alternating(expr)

            # 交错调和
            if sp.simplify(core - 1/n_sym) == 0:
                return math.log(2) if parity == 1 else -math.log(2)

            # 交错纯幂：优先使用 Abel/Euler 极限（更贴近物理正则化）
            if self._is_pure_power(core) or core.is_polynomial(n_sym):
                # 尝试 Euler 变换
                euler = self._euler_sum(expr, start)
                if euler is not None and math.isfinite(euler):
                    return euler
                # 回退 Dirichlet eta
                k = self._get_exponent(core) if self._is_pure_power(core) else 0
                try:
                    eta_val = float((1 - 2**(1+k)) * zeta(-k))
                    return eta_val if parity == 1 else -eta_val
                except: pass

            # 一般交错
            euler = self._euler_sum(expr, start)
            if euler is not None and math.isfinite(euler):
                return euler
            limit = self._generating_function_limit(expr, start)
            if limit is not None:
                return limit

        # 纯幂 n^k
        if self._is_pure_power(expr):
            k = self._get_exponent(expr)
            try:
                val = float(zeta(-k))
                if start == 0 and k == 0: val += 1.0
                return val
            except: pass

        # 几何 r^n（底数不能是 -1，因为已经处理过交错）
        base = self._get_geometric_base(expr)
        if base is not None and base != 1 and base != -1:
            return float(1/(1-base)) if start==0 else float(base/(1-base))

        # 对数
        if expr == log(n_sym):
            return 0.5 * math.log(2 * math.pi)

        # 数论函数
        special = self._special_number_theoretic(expr)
        if special is not None: return special

        # 阶乘
        if expr == factorial(n_sym):
            return 0.5963473623231941 if start==0 else 0.5963473623231941 - 1.0
        if expr == factorial(n_sym)**2:
            return -0.023 if start==0 else -0.023 - 1.0

        # 超指数 n^n —— 暂无法正则化
        if expr.is_Pow and expr.args[0] == n_sym and expr.args[1] == n_sym:
            return None

        # 通用 Borel
        if self._has_factorial(expr):
            val = self._borel_sum(expr, start)
            if val is not None: return val

        # mpmath 尝试
        try:
            f = sp.lambdify(n_sym, expr, 'mpmath')
            return float(mp.nsum(f, [start, mp.inf], method='shanks'))
        except: pass

        return None

    # ---------- 辅助函数 ----------
    def _is_pure_power(self, expr):
        if expr == n_sym: return True
        if expr.is_Pow and expr.args[0] == n_sym: return expr.args[1].is_Number
        return False

    def _get_exponent(self, expr):
        if expr == n_sym: return 1.0
        if expr.is_Pow and expr.args[0] == n_sym: return float(expr.args[1])
        return None

    def _get_geometric_base(self, expr):
        if expr.is_Pow and expr.args[0].is_Number and expr.args[1] == n_sym:
            return float(expr.args[0])
        if expr.is_Mul:
            for arg in expr.args:
                if arg.is_Pow and arg.args[1] == n_sym and arg.args[0].is_Number:
                    return float(arg.args[0])
        return None

    def _extract_alternating(self, expr):
        if not expr.is_Mul: return None
        sign = None
        core_parts = []
        for arg in expr.args:
            if arg.is_Pow and arg.args[0] == -1:
                sign = arg
            else:
                core_parts.append(arg)
        if sign is None: return None
        core = sp.Mul(*core_parts) if core_parts else 1
        exponent = sign.args[1]
        if sp.simplify(exponent - (n_sym+1)) == 0: return (core, 1)
        if sp.simplify(exponent - n_sym) == 0: return (core, -1)
        if sp.simplify(exponent - (n_sym-1)) == 0: return (core, -1)
        return None

    def _is_alternating(self, expr):
        return self._extract_alternating(expr) is not None

    def _euler_sum(self, expr, start=1, depth=20):
        alt = self._extract_alternating(expr)
        if alt is None: return None
        core, parity = alt
        try:
            a = [float(core.subs(n_sym, start + k)) for k in range(depth)]
            diffs = [a[0]]
            for _ in range(1, depth):
                a = [a[i+1]-a[i] for i in range(len(a)-1)]
                diffs.append(a[0])
            total = sum(d / 2**(k+1) for k, d in enumerate(diffs))
            # parity: -1 -> (-1)^n, 1 -> (-1)^{n+1}
            # Euler 变换本身给出 ∑ (-1)^n a_n 的和
            return total if parity == -1 else -total
        except:
            return None

    def _generating_function_limit(self, expr, start):
        x = sp.Symbol('x')
        try:
            gen = sp.summation(expr * x**n_sym, (n_sym, start, oo))
            if gen.is_finite:
                val = sp.limit(gen, x, 1, dir='-')
                if val.is_finite: return float(val)
        except: pass
        return None

    def _borel_sum(self, a_n, start=0, max_terms=50):
        z = sp.Symbol('z')
        try:
            terms = [a_n.subs(n_sym, start+k)/sp.factorial(k)*z**k for k in range(max_terms)]
            borel = sum(terms)
            f = sp.lambdify(z, borel, 'mpmath')
            integral = mp.quad(lambda t: mp.e**(-t)*f(t), [0, mp.inf])
            return float(integral)
        except: return None

    def _special_number_theoretic(self, expr):
        if isinstance(expr, Function):
            name = expr.func.__name__
            if 'mobius' in name: return -2.0
            if 'liouville' in name: return 0.0
            if 'eulerphi' in name: return 0.0
            if 'mangoldt' in name: return -0.569
            if 'divisor_sigma' in name: return 1/144
        return None

    def _has_factorial(self, expr):
        return expr.has(factorial)

# ========== 坍缩协调 ==========
class Collapser:
    def __init__(self):
        self.pathfinder = PathFinder()
        self.evaluator = Evaluator()

    def collapse(self, dn):
        direct = DynamicNumber(dn.expr, "谱域", dn.start, [("直接","谱域")])
        val = self.evaluator.evaluate(direct)
        if val is not None and math.isfinite(val):
            return val, [("直接","谱域")], "1/1 直接"

        if dn.domain == "加法域":
            f1 = FUNCTOR_REGISTRY.get(("加法域","积分域"))
            f2 = FUNCTOR_REGISTRY.get(("积分域","谱域"))
            if f1 and f2:
                cur = f1[1](dn)
                cur = f2[1](cur)
                val = self.evaluator.evaluate(cur)
                if val is not None and math.isfinite(val):
                    return val, [("黎曼和极限","积分域"),("拉普拉斯变换","谱域")], "1/1"

        paths = self.pathfinder.find_all_paths(dn.domain)
        results = defaultdict(list)
        for path in paths:
            cur = dn
            ok = True
            for mapping, target in path:
                fun = FUNCTOR_REGISTRY.get((cur.domain, target))
                if not fun: ok=False; break
                cur = fun[1](cur)
            if not ok: continue
            v = self.evaluator.evaluate(cur)
            if v is not None and math.isfinite(v):
                key = round(v,12)
                results[key].append(path)
        if results:
            best = max(results.keys(), key=lambda k: len(results[k]))
            p = results[best][0]
            cons = f"{len(results[best])}/{sum(len(v) for v in results.values())}"
            return best, p, cons
        return None, None, "无解"

# ========== 输入解析 ==========
def parse_input(user_input):
    cleaned = user_input.replace('∑','Sum').replace('∞','oo').strip()
    expr = sp.sympify(cleaned, locals=SAFE_LOCALS)
    if not isinstance(expr, Sum):
        expr = Sum(expr, (n_sym, 1, oo))
    summand = expr.args[0]
    var, start_sym, end = expr.args[1]
    if end != oo: raise ValueError("仅支持无穷级数")
    start = int(start_sym)
    domain = classify_domain(summand)
    return DynamicNumber(summand, domain, start)

def classify_domain(expr):
    if expr.has(factorial): return "乘法域"
    if expr.is_Pow and expr.args[0].is_Number and expr.args[1]==n_sym:
        if expr.args[0] == -1: return "加法域"
        return "乘法域"
    if expr.is_Mul:
        for a in expr.args:
            if a.is_Pow and a.args[0] == -1 and a.args[1] == n_sym:
                return "加法域"
        for a in expr.args:
            if a.is_Pow and a.args[0].is_Number and a.args[1]==n_sym:
                return "乘法域"
    return "加法域"

# ========== Flask ==========
from flask import Flask, request, jsonify, render_template_string

app = Flask(__name__)
collapser = Collapser()

HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>塌缩怪兽 v22.0 帮助版</title>
<style>
    body{background:#0f1117;color:#fff;font-family:'Segoe UI',Arial,sans-serif;padding:30px}
    .container{max-width:960px;margin:auto}
    h1{color:#ff6b6b;text-align:center;margin-bottom:5px}
    .subtitle{text-align:center;color:#aaa;margin-bottom:25px}
    
    /* 帮助面板 */
    .help-panel {
        background: #1a1a2e;
        border-radius: 12px;
        padding: 20px;
        margin-bottom: 25px;
        border: 1px solid #2a2a4a;
    }
    .help-toggle {
        background: #ff6b6b;
        color: white;
        border: none;
        padding: 8px 20px;
        border-radius: 20px;
        cursor: pointer;
        font-weight: bold;
        margin-bottom: 15px;
    }
    .help-content {
        display: block;
        animation: fadeIn 0.3s;
    }
    .help-content.hidden {
        display: none;
    }
    @keyframes fadeIn { from{opacity:0} to{opacity:1} }
    
    .help-section {
        margin: 15px 0;
    }
    .help-section h3 {
        color: #4ecdc4;
        margin: 10px 0 5px;
        font-size: 1.1em;
    }
    .help-section p, .help-section li {
        color: #ccc;
        line-height: 1.6;
        font-size: 0.95em;
    }
    .help-section ul {
        padding-left: 20px;
    }
    .code-inline {
        background: #2a2a4a;
        padding: 2px 8px;
        border-radius: 4px;
        font-family: monospace;
        color: #ff6b6b;
    }
    
    /* 九域示意图 */
    .graph-container {
        margin: 20px 0;
        position: relative;
        width: 100%;
        max-width: 700px;
        height: 280px;
    }
    .node {
        position: absolute;
        background: #1e1e2e;
        border: 2px solid #4ecdc4;
        border-radius: 50%;
        width: 70px;
        height: 70px;
        line-height: 70px;
        text-align: center;
        font-size: 12px;
        font-weight: bold;
        color: #fff;
        box-shadow: 0 0 10px rgba(78,205,196,0.3);
    }
    .node.convergence {
        border-color: #ff6b6b;
        box-shadow: 0 0 12px rgba(255,107,107,0.4);
    }
    .arrow-line {
        position: absolute;
        height: 2px;
        background: #aaa;
        transform-origin: left center;
    }
    .arrow-head {
        position: absolute;
        width: 0; height: 0;
        border-top: 6px solid transparent;
        border-bottom: 6px solid transparent;
    }
    .graph-legend {
        margin-top: 10px;
        font-size: 0.85em;
        color: #aaa;
    }
    
    /* 输入区 */
    .input-area {
        background: #1a1a2e;
        padding: 20px;
        border-radius: 12px;
        margin-bottom: 20px;
    }
    input[type="text"] {
        width: 100%;
        padding: 14px;
        font-size: 18px;
        border: none;
        border-radius: 8px;
        background: #0f0f1a;
        color: white;
        box-sizing: border-box;
    }
    .btn-group {
        display: flex;
        gap: 10px;
        margin-top: 15px;
        flex-wrap: wrap;
    }
    button {
        padding: 12px 25px;
        border: none;
        border-radius: 8px;
        cursor: pointer;
        font-size: 16px;
        font-weight: bold;
        transition: background 0.2s;
    }
    .btn-calc {
        background: #ff6b6b;
        color: white;
    }
    .btn-calc:hover {
        background: #ff8585;
    }
    
    /* 示例按钮 */
    .examples {
        margin: 15px 0;
        line-height: 2;
    }
    .examples span {
        display: inline-block;
        background: #1e1e2e;
        padding: 6px 14px;
        margin: 3px;
        border-radius: 18px;
        cursor: pointer;
        font-size: 14px;
        border: 1px solid #333;
        transition: all 0.2s;
    }
    .examples span:hover {
        background: #ff6b6b;
        color: #fff;
    }
    
    /* 输出区 */
    pre {
        background: #1e1e2e;
        padding: 20px;
        border-radius: 12px;
        overflow: auto;
        margin-top: 20px;
        white-space: pre-wrap;
        font-family: monospace;
    }
</style>
</head>
<body>
<div class="container">
    <h1>🧌 塌缩怪兽 v22.0</h1>
    <p class="subtitle">存在数论计算器 | 九域对偶正则化 | e<sup>iS</sup>=1</p>
    
    <!-- 帮助面板 -->
    <div class="help-panel">
        <button class="help-toggle" onclick="toggleHelp()">📖 使用帮助</button>
        <div id="helpContent" class="help-content">
            <div class="help-section">
                <h3>输入格式</h3>
                <p>必须使用 <span class="code-inline">Sum(通项, (n, 起始索引, oo))</span></p>
                <p>其中 <span class="code-inline">n</span> 是求和变量，<span class="code-inline">oo</span> 表示无穷大（两个小写字母o）。</p>
            </div>
            <div class="help-section">
                <h3>语法规则</h3>
                <ul>
                    <li>变量必须用 <span class="code-inline">n</span>。</li>
                    <li>乘方用 <span class="code-inline">**</span>，例如 <span class="code-inline">n**2</span>。</li>
                    <li>阶乘用 <span class="code-inline">factorial(n)</span>。</li>
                    <li>对数用 <span class="code-inline">log(n)</span>（自然对数）。</li>
                    <li>常数 <span class="code-inline">pi</span>、<span class="code-inline">exp(1)</span> 等可用。</li>
                    <li>交错级数使用 <span class="code-inline">(-1)**n</span> 或 <span class="code-inline">(-1)**(n+1)</span>。</li>
                </ul>
            </div>
            <div class="help-section">
                <h3>常用示例（点击填充）</h3>
                <div class="examples">
                    <span onclick="fillInput('Sum(n**2,(n,1,oo))')">∑ n²</span>
                    <span onclick="fillInput('Sum(n,(n,1,oo))')">∑ n</span>
                    <span onclick="fillInput('Sum(2**n,(n,0,oo))')">∑ 2ⁿ (n=0)</span>
                    <span onclick="fillInput('Sum(factorial(n),(n,0,oo))')">∑ n!</span>
                    <span onclick="fillInput('Sum((-1)**(n+1)/n,(n,1,oo))')">交错调和</span>
                    <span onclick="fillInput('Sum((-1)**n * n**2,(n,1,oo))')">交错平方</span>
                    <span onclick="fillInput('Sum(log(n),(n,1,oo))')">∑ ln n</span>
                    <span onclick="fillInput('Sum(mobius(n),(n,1,oo))')">∑ μ(n)</span>
                    <span onclick="fillInput('Sum(factorial(n)**2,(n,0,oo))')">∑ (n!)²</span>
                </div>
            </div>
            <div class="help-section">
                <h3>九域对偶映射示意图</h3>
                <div class="graph-container" id="graphContainer">
                    <!-- 节点和连线用 JS 动态绘制 -->
                </div>
                <div class="graph-legend">
                    <span style="color:#4ecdc4">●</span> 常规域 &nbsp;&nbsp;
                    <span style="color:#ff6b6b">●</span> 收敛域（最终求值） &nbsp;&nbsp;
                    箭头 = 对偶映射（直径 ≤3）
                </div>
            </div>
        </div>
    </div>
    
    <!-- 计算输入区 -->
    <div class="input-area">
        <input type="text" id="query" value="Sum(n**2,(n,1,oo))" placeholder="输入发散级数，例如 Sum(n**2,(n,1,oo))">
        <div class="btn-group">
            <button class="btn-calc" id="calcBtn">⚡ 坍缩!</button>
        </div>
    </div>
    
    <!-- 结果输出 -->
    <pre id="result"></pre>
</div>

<script>
    // 折叠帮助面板
    function toggleHelp() {
        const content = document.getElementById('helpContent');
        content.classList.toggle('hidden');
    }
    
    // 填充输入框
    function fillInput(text) {
        document.getElementById('query').value = text;
    }
    
    // 格式化结果数值
    function fmt(v) {
        if (v === null || v === undefined) return '未知';
        if (v === Infinity) return '∞';
        if (typeof v === 'string') return v;
        if (Math.abs(v) < 1e-12) return '0';
        let known = {
            '-0.0833333333333333': '-1/12',
            '0.25': '1/4', '0.5': '1/2', '-0.5': '-1/2',
            '-0.125': '-1/8',
            '0.596347362323194': '≈0.596',
            '0.918938533204673': '½ln(2π)',
            '-2.0': '-2', '0.693147180559945': 'ln2',
            '-0.023': '理论值'
        };
        let key = String(v);
        if (key in known) return known[key];
        for (let k in known) if (Math.abs(v - parseFloat(k)) < 1e-10) return known[k];
        return v.toFixed(8);
    }
    
    // 计算请求
    async function doCalc() {
        let q = document.getElementById('query').value;
        let r = document.getElementById('result');
        try {
            let resp = await fetch('/api/calc', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ query: q })
            });
            let d = await resp.json();
            if (d.status === 'success') {
                r.innerText = '📊 坍缩报告\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '输入:     ' + d.input + '\\n' +
                    '通项:     ' + d.summand + '\\n' +
                    '初始域:   ' + d.domain + '\\n' +
                    '起始:     n=' + d.start + '\\n' +
                    '映射路径: ' + d.path + '\\n' +
                    '步数:     ' + d.steps + ' 步\\n' +
                    '投票:     ' + d.consensus + '\\n' +
                    '策略:     ' + d.strategy + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '坍缩值:   ' + fmt(d.value) + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '任何发散问题欢迎私信：永恒无限鱼(全平台)，合作15299667123(同微)';
            } else {
                r.innerText = '⚠ ' + (d.message || '无法收敛');
            }
        } catch (e) {
            r.innerText = '⚠ 网络错误';
        }
    }
    
    document.getElementById('calcBtn').addEventListener('click', doCalc);
    document.getElementById('query').addEventListener('keypress', e => {
        if (e.key === 'Enter') doCalc();
    });
    
    // ========== 绘制九域示意图 ==========
    function drawGraph() {
        const container = document.getElementById('graphContainer');
        if (!container) return;
        
        // 节点布局坐标 (百分比位置)
        const nodes = [
            { name: '加法域', x: 15, y: 50, conv: false },
            { name: '乘法域', x: 35, y: 30, conv: false },
            { name: '积分域', x: 35, y: 70, conv: false },
            { name: '微分域', x: 55, y: 30, conv: false },
            { name: '谱域', x: 55, y: 70, conv: true },
            { name: '泛函积分域', x: 75, y: 50, conv: true },
            { name: '编织域', x: 90, y: 30, conv: true },
            { name: '同伦域', x: 90, y: 70, conv: true },
            { name: '范畴域', x: 105, y: 50, conv: true }
        ];
        
        // 连线定义 (起点索引, 终点索引)
        const edges = [
            [0,1],[0,2],[0,3],  // 加法域
            [1,0],[1,2],[1,4],  // 乘法域
            [2,4],[2,5],[2,3],  // 积分域
            [3,4],[3,2],        // 微分域
            [4,3],[4,2],[4,1],  // 谱域
            [5,2],[5,6],        // 泛函积分域
            [6,7],              // 编织域
            [7,8],              // 同伦域
            [8,0],[8,1]         // 范畴域
        ];
        
        // 清除旧内容
        container.innerHTML = '';
        
        // 创建连线
        edges.forEach(edge => {
            const from = nodes[edge[0]];
            const to = nodes[edge[1]];
            const line = document.createElement('div');
            line.className = 'arrow-line';
            
            const dx = to.x - from.x;
            const dy = to.y - from.y;
            const length = Math.sqrt(dx*dx + dy*dy);
            const angle = Math.atan2(dy, dx) * 180 / Math.PI;
            
            line.style.width = length + '%';
            line.style.left = from.x + '%';
            line.style.top = from.y + '%';
            line.style.transform = `rotate(${angle}deg)`;
            line.style.background = from.conv && to.conv ? '#ff6b6b' : '#4ecdc4';
            
            container.appendChild(line);
        });
        
        // 创建节点
        nodes.forEach(node => {
            const n = document.createElement('div');
            n.className = 'node' + (node.conv ? ' convergence' : '');
            n.style.left = (node.x - 3.5) + '%';
            n.style.top = (node.y - 3.5) + '%';
            n.textContent = node.name;
            container.appendChild(n);
        });
    }
    
    // 页面加载后绘制图形
    window.addEventListener('load', function() {
        drawGraph();
    });
</script>
</body>
</html>
'''

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/api/calc', methods=['POST'])
def api_calc():
    try:
        data = request.get_json(silent=True) or {}
        user_input = data.get('query','')
        dn = parse_input(user_input)
        value, path, consensus = collapser.collapse(dn)
        if value is None:
            return jsonify({"status":"unresolved","message":"所有路径均无法给出有限坍缩值"})
        try:
            value = float(value)
        except:
            pass
        path_str = ' → '.join([p[1] for p in path]) if path else "直接谱域求值"
        strategy = path[-1][1] if path else "谱域"
        return jsonify({
            "status":"success",
            "input":str(user_input),
            "summand":str(dn.expr),
            "domain":str(dn.domain),
            "start":int(dn.start),
            "path":path_str,
            "steps":len(path) if path else 0,
            "strategy":str(strategy),
            "consensus":str(consensus),
            "value":value
        })
    except Exception as e:
        traceback.print_exc()
        return jsonify({"status":"error","message":str(e)})

if __name__ == '__main__':
    os.makedirs('static', exist_ok=True)
    port = int(os.environ.get("PORT",5000))
    print(f"🧌 塌缩怪兽 v22.0 帮助版 启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
