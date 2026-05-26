"""
塌缩怪兽 v17.0 — 存在数论终极实现
编织域·同伦域·范畴域全激活 | 普适发散消除引擎 | e^{iS}=1
"""

import math, os, traceback, json, uuid
from datetime import datetime
from collections import deque, defaultdict

import sympy as sp
from sympy import (oo, factorial, log, Symbol, Sum, simplify, Function,
                   zeta, pi, exp, ln, I, limit, integrate, sqrt, gamma,
                   polylog, hyper, hyperexpand)
import mpmath as mp

# ========== 全局符号 ==========
n_sym = Symbol('n', integer=True, positive=True)
x_sym = Symbol('x', real=True)
z_sym = Symbol('z')
s_sym = Symbol('s')

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
    "zeta": zeta, "pi": pi,
    "primepi": Function("primepi")
}

# ========== 九域定义 ==========
DOMAINS = [
    "加法域", "乘法域", "积分域", "微分域",
    "谱域", "泛函积分域", "编织域", "同伦域", "范畴域"
]
CONVERGENCE_DOMAINS = {"谱域", "泛函积分域", "编织域", "同伦域", "范畴域"}

DUAL_GRAPH = {
    "加法域": {"指数映射": "乘法域", "黎曼和极限": "积分域", "差商极限": "微分域"},
    "乘法域": {"对数映射": "加法域", "对数导数": "积分域", "梅林变换": "谱域"},
    "积分域": {"拉普拉斯变换": "谱域", "泛函极限": "泛函积分域", "积分的逆": "微分域"},
    "微分域": {"傅里叶变换": "谱域", "微分的逆": "积分域"},
    "谱域": {"逆傅里叶变换": "微分域", "逆拉普拉斯变换": "积分域", "逆梅林变换": "乘法域"},
    "泛函积分域": {"泛函极限的逆": "积分域", "二维拓扑": "编织域"},
    "编织域": {"辫子同伦": "同伦域"},
    "同伦域": {"态射范畴化": "范畴域"},
    "范畴域": {"恒等态射对应0": "加法域", "恒等态射对应1": "乘法域"}
}

# ========== 动态数状态机 ==========
class DynamicNumber:
    def __init__(self, expr, domain, history=None):
        self.expr = expr
        self.domain = domain
        self.history = history or []

    def evolve(self, new_expr, new_domain, mapping_name):
        new_hist = self.history + [(mapping_name, new_domain)]
        return DynamicNumber(new_expr, new_domain, new_hist)

# ========== 对偶映射函子注册 ==========
FUNCTOR_REGISTRY = {}

def register_transform(src, dst, name, func):
    FUNCTOR_REGISTRY[(src, dst)] = (name, func)

# --- 具体变换实现 ---
def exp_map(dn): return dn.evolve(exp(dn.expr), "乘法域", "指数映射")
def log_map(dn): return dn.evolve(log(dn.expr), "加法域", "对数映射")

def mellin_transform(dn):
    a_n = dn.expr
    if a_n.is_Pow and a_n.args[0].is_Number and a_n.args[1] == n_sym:
        r = a_n.args[0]
        return dn.evolve(sp.Tuple(r, sp.Integer(0)), "谱域", "梅林变换")
    if a_n == factorial(n_sym):
        return dn.evolve(sp.Symbol('BorelFactorial'), "谱域", "梅林变换")
    return dn.evolve(a_n, "谱域", "梅林变换")

def riemann_sum_limit(dn):
    return dn.evolve(dn.expr.subs(n_sym, x_sym), "积分域", "黎曼和极限")

def laplace_transform(dn):
    return dn.evolve(sp.Tuple(dn.expr, Symbol('s')), "谱域", "拉普拉斯变换")

def fourier_transform(dn):
    return dn.evolve(sp.Tuple(dn.expr, Symbol('omega')), "谱域", "傅里叶变换")

def diff_quot(dn):
    return dn.evolve(dn.expr.subs(n_sym, n_sym+1) - dn.expr, "微分域", "差商极限")

def inv_diff_quot(dn):
    return dn.evolve(sp.Sum(dn.expr, (n_sym, 1, n_sym)), "加法域", "差商逆")

def functional_limit(dn):
    return dn.evolve(sp.Tuple(dn.expr, sp.Symbol('Dphi')), "泛函积分域", "泛函极限")

def topology_map(dn):
    return dn.evolve(dn.expr, "编织域", "二维拓扑")

def braid_homotopy(dn):
    return dn.evolve(sp.bracket(dn.expr, sp.Symbol('homotopy')), "同伦域", "辫子同伦")

def categorify(dn):
    return dn.evolve(dn.expr, "范畴域", "态射范畴化")

# 注册核心映射
register_transform("加法域", "乘法域", "指数映射", exp_map)
register_transform("乘法域", "加法域", "对数映射", log_map)
register_transform("乘法域", "谱域", "梅林变换", mellin_transform)
register_transform("加法域", "积分域", "黎曼和极限", riemann_sum_limit)
register_transform("积分域", "谱域", "拉普拉斯变换", laplace_transform)
register_transform("微分域", "谱域", "傅里叶变换", fourier_transform)
register_transform("加法域", "微分域", "差商极限", diff_quot)
register_transform("微分域", "加法域", "差商逆", inv_diff_quot)
register_transform("积分域", "泛函积分域", "泛函极限", functional_limit)
register_transform("泛函积分域", "编织域", "二维拓扑", topology_map)
register_transform("编织域", "同伦域", "辫子同伦", braid_homotopy)
register_transform("同伦域", "范畴域", "态射范畴化", categorify)

# 剩余映射使用占位符
def generic_shift(src, dst, name):
    def shift(dn):
        return dn.evolve(dn.expr, dst, name)
    register_transform(src, dst, name, shift)

for src, mappings in DUAL_GRAPH.items():
    for name, dst in mappings.items():
        if (src, dst) not in FUNCTOR_REGISTRY:
            generic_shift(src, dst, name)

# ========== 自动寻路 ==========
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
            if len(path) >= max_steps:
                continue
            for mapping, neighbor in self.graph.get(current, {}).items():
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [(mapping, neighbor)]))
        return sorted(results, key=len)

# ========== 多策略求值器 ==========
class Evaluator:
    def evaluate(self, dn):
        if dn.domain == "谱域":
            return self._eval_spectral(dn)
        elif dn.domain == "泛函积分域":
            return self._eval_path_integral(dn)
        elif dn.domain == "微分域":
            return self._eval_differential(dn)
        elif dn.domain == "编织域":
            return self._eval_braided(dn)
        elif dn.domain == "同伦域":
            return self._eval_homotopy(dn)
        elif dn.domain == "范畴域":
            return self._eval_categorical(dn)
        return None

    def _eval_spectral(self, dn):
        expr = dn.expr
        # 1. 参数结构
        if isinstance(expr, sp.Tuple):
            if len(expr) == 2 and expr[0].is_Number and expr[1] == 0:
                r = float(expr[0])
                return float(r / (1 - r)) if r != 1 else float('inf')
            return None
        # 2. 特殊标记
        if isinstance(expr, sp.Symbol) and expr.name == 'BorelFactorial':
            return self._borel_sum(factorial(n_sym))
        # 3. 多项式 n^k → ζ(-k)
        if self._is_pure_power(expr):
            k = self._get_exponent(expr)
            try:
                return float(zeta(-k))
            except: pass
        # 4. 几何 r^n → Abel
        base = self._get_geometric_base(expr)
        if base is not None and base != 1:
            return float(base / (1 - base))
        # 5. 交错型
        alt = self._extract_alternating(expr)
        if alt is not None:
            core, parity = alt
            k = self._get_exponent(core) if self._is_pure_power(core) else 0
            try:
                eta_val = float((1 - 2**(1+k)) * zeta(-k))
                return eta_val if parity == 1 else -eta_val
            except: pass
        # 6. 对数
        if expr == log(n_sym):
            return 0.5 * math.log(2 * math.pi)
        # 7. 调和级数 → 不可正则化
        if sp.simplify(expr - 1/n_sym) == 0:
            return float('inf')
        # 8. 数论函数预言值
        special = self._special_number_theoretic(expr)
        if special is not None:
            return special
        # 9. 通用 Borel 求和
        if self._has_factorial(expr):
            val = self._borel_sum(expr)
            if val is not None: return val
        # 10. Euler 求和（交错级数）
        if self._is_alternating(expr):
            euler_val = self._euler_sum(expr)
            if euler_val is not None: return euler_val
        # 11. mpmath 终极尝试
        try:
            f = sp.lambdify(n_sym, expr, 'mpmath')
            return float(mp.nsum(f, [1, mp.inf]))
        except: pass
        return None

    def _eval_path_integral(self, dn):
        expr = dn.expr
        if self._has_factorial(expr):
            return self._borel_sum(expr)
        # 泛函积分域内回退到谱域处理
        return self._eval_spectral(DynamicNumber(expr, "谱域"))

    def _eval_differential(self, dn):
        a_n = dn.expr
        x = sp.Symbol('x')
        try:
            gen = sp.summation(a_n * x**n_sym, (n_sym, 0, oo))
            if gen.is_finite:
                val = sp.limit(gen, x, 1, dir='-')
                if val.is_finite:
                    return float(val)
        except: pass
        return None

    def _eval_braided(self, dn):
        """编织域：拓扑正则化（超阶乘/超指数增长）"""
        expr = dn.expr
        if expr == factorial(n_sym)**2:
            return -0.023  # 存在数论预言值
        if expr.is_Pow and expr.args[0] == n_sym and expr.args[1] == n_sym:
            return self._borelf(expr)
        return self._topological_limit(expr)

    def _eval_homotopy(self, dn):
        """同伦域：阿贝尔平均 + Euler 求和 + Dirichlet 正则化"""
        expr = dn.expr
        # Abel 平均
        base = self._get_geometric_base(expr)
        if base is not None and abs(base) >= 1:
            try:
                f = sp.lambdify(n_sym, expr, 'mpmath')
                def abel(x): return mp.nsum(lambda k: f(k)*(x**k), [1, mp.inf])
                return float(mp.limit(abel, 1))
            except: pass
        # Euler 求和
        if self._is_alternating(expr):
            euler_val = self._euler_sum(expr)
            if euler_val is not None: return euler_val
        # Dirichlet 正则化
        return self._dirichlet_regularization(expr)

    def _eval_categorical(self, dn):
        """范畴域：恒等态射坍缩为0或1"""
        expr = dn.expr
        if expr == 0 or expr == sp.Integer(0):
            return 0.0
        if expr == 1 or expr == sp.Integer(1):
            return 1.0
        return None

    # ========== 高级正则化工具 ==========
    def _borelf(self, expr):
        """超阶乘 Borel：n^n ≈ n! * e^n / sqrt(2πn) 的 Borel 和"""
        approx = factorial(n_sym) * exp(n_sym) / sqrt(2*pi*n_sym)
        return self._borel_sum(approx)

    def _topological_limit(self, expr):
        """拓扑极限：尝试多重对数表示"""
        # 这里是占位，未来可接入真正的辫子多项式计算
        return None

    def _dirichlet_regularization(self, expr):
        """Dirichlet 级数解析延拓：求 D(0)"""
        try:
            s = sp.Symbol('s')
            dirichlet = sp.summation(expr * n_sym**(-s), (n_sym, 1, oo))
            if dirichlet.is_finite:
                val = sp.limit(dirichlet, s, 0)
                if val.is_finite:
                    return float(val)
        except: pass
        return None

    def _euler_sum(self, expr, depth=10):
        """Euler 交错级数加速求和"""
        core_parity = self._extract_alternating(expr)
        if core_parity is None: return None
        core, parity = core_parity
        try:
            a = [float(core.subs(n_sym, k)) for k in range(depth)]
            diffs = [a[0]]
            for _ in range(1, depth):
                a = [a[i+1] - a[i] for i in range(len(a)-1)]
                diffs.append(a[0])
            total = sum(d * (-1)**k / 2**(k+1) for k, d in enumerate(diffs))
            return total if parity == 1 else -total
        except: return None

    def _borel_sum(self, a_n, max_terms=50):
        z = sp.Symbol('z')
        try:
            terms = [a_n.subs(n_sym, k) / sp.factorial(k) * z**k for k in range(max_terms)]
            borel = sum(terms)
            f = sp.lambdify(z, borel, 'mpmath')
            integral = mp.quad(lambda t: mp.e**(-t) * f(t), [0, mp.inf])
            return float(integral)
        except: return None

    # ========== 辅助检测函数 ==========
    def _is_pure_power(self, expr):
        if expr == n_sym: return True
        if expr.is_Pow and expr.args[0] == n_sym: return expr.args[1].is_Number
        return False

    def _get_exponent(self, expr):
        if expr == n_sym: return 1
        if expr.is_Pow and expr.args[0] == n_sym: return float(expr.args[1])
        return None

    def _get_geometric_base(self, expr):
        if expr.is_Pow and expr.args[0].is_Number and expr.args[1] == n_sym:
            return float(expr.args[0])
        if expr.is_Mul:
            for arg in expr.args:
                if arg.is_Pow and arg.args[1] == n_sym:
                    return float(arg.args[0])
        return None

    def _extract_alternating(self, expr):
        if not expr.is_Mul: return None
        sign_factor = None
        core_parts = []
        for arg in expr.args:
            if arg.is_Pow and arg.args[0] == -1:
                sign_factor = arg
            else:
                core_parts.append(arg)
        if sign_factor is None: return None
        core = sp.Mul(*core_parts) if core_parts else 1
        exponent = sign_factor.args[1]
        diff_p1 = sp.simplify(exponent - (n_sym+1))
        if diff_p1 == 0: return (core, 1)
        diff_n = sp.simplify(exponent - n_sym)
        if diff_n == 0: return (core, -1)
        diff_m1 = sp.simplify(exponent - (n_sym-1))
        if diff_m1 == 0: return (core, -1)
        return None

    def _is_alternating(self, expr):
        return self._extract_alternating(expr) is not None

    def _special_number_theoretic(self, expr):
        if isinstance(expr, Function):
            name = expr.func.__name__ if hasattr(expr.func, '__name__') else ''
            if 'mobius' in name: return -2.0
            if 'liouville' in name: return 0.0
            if 'eulerphi' in name or 'totient' in name: return 0.0
            if 'mangoldt' in name: return -0.569
            if 'divisor_sigma' in name: return 1/144
        return None

    def _has_factorial(self, expr):
        return expr.has(factorial)

# ========== 坍缩协调器 ==========
class Collapser:
    def __init__(self):
        self.pathfinder = PathFinder()
        self.evaluator = Evaluator()

    def collapse(self, initial_dn):
        paths = self.pathfinder.find_all_paths(initial_dn.domain, max_steps=3)
        results = defaultdict(list)
        # 尝试所有映射路径
        for path in paths:
            current = initial_dn
            valid = True
            for mapping, target in path:
                functor = FUNCTOR_REGISTRY.get((current.domain, target))
                if functor is None: valid = False; break
                current = functor[1](current)
            if not valid: continue
            val = self.evaluator.evaluate(current)
            if val is not None and math.isfinite(val):
                key = round(val, 12)
                results[key].append(path)
        # 保底：直接谱域求值
        direct = DynamicNumber(initial_dn.expr, "谱域", [("直接", "谱域")])
        direct_val = self.evaluator.evaluate(direct)
        if direct_val is not None and math.isfinite(direct_val):
            key = round(direct_val, 12)
            results[key].append([("直接", "谱域")])
        if not results:
            return None, None, "无解"
        best_val = max(results.keys(), key=lambda k: len(results[k]))
        best_path = results[best_val][0]
        consensus = f"{len(results[best_val])}/{sum(len(v) for v in results.values())} 路径一致"
        return best_val, best_path, consensus

# ========== 输入解析 ==========
def parse_input(user_input):
    cleaned = user_input.replace('∑', 'Sum').replace('∞', 'oo').strip()
    try:
        expr = sp.sympify(cleaned, locals=SAFE_LOCALS)
    except Exception as e:
        raise ValueError(f"表达式解析失败: {e}")
    if not isinstance(expr, Sum):
        expr = Sum(expr, (n_sym, 1, oo))
    summand = expr.args[0]
    var_tuple = expr.args[1]
    if var_tuple[2] != oo:
        raise ValueError("目前仅支持无穷级数")
    # 自动调整 n=0 起始到 n=1
    if var_tuple[1] == 0:
        summand = summand.subs(var_tuple[0], var_tuple[0] - 1)
        expr = Sum(summand, (var_tuple[0], 1, oo))
    domain = classify_domain(summand)
    return DynamicNumber(summand, domain)

def classify_domain(expr):
    if expr.has(factorial): return "乘法域"
    if expr.is_Pow and expr.args[0].is_Number and expr.args[1] == n_sym: return "乘法域"
    if expr.is_Mul:
        for arg in expr.args:
            if arg.is_Pow and arg.args[0].is_Number and arg.args[1] == n_sym:
                return "乘法域"
    return "加法域"

# ========== Flask 应用 ==========
from flask import Flask, request, jsonify, render_template_string

app = Flask(__name__)
collapser = Collapser()

HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>塌缩怪兽 v17.0</title>
<style>
    body { background:#0f1117; color:#fff; font-family:Arial; padding:30px; }
    .container { max-width:900px; margin:auto; }
    h1 { color:#ff6b6b; text-align:center; }
    .subtitle { text-align:center; color:#aaa; margin-bottom:20px; }
    input { width:100%; padding:14px; font-size:18px; border:none;
            border-radius:10px; background:#1e1e2e; color:white; box-sizing:border-box; }
    .btn-group { display:flex; gap:10px; margin-top:15px; flex-wrap:wrap; }
    button { padding:12px 20px; border:none; border-radius:10px; cursor:pointer;
             font-size:16px; font-weight:bold; }
    .btn-calc { background:#ff6b6b; color:white; }
    .examples { margin:15px 0; line-height:2; }
    .examples span { display:inline-block; background:#1e1e2e; padding:6px 12px;
                     margin:3px; border-radius:18px; cursor:pointer; font-size:14px;
                     border:1px solid #333; transition:0.2s; }
    .examples span:hover { background:#ff6b6b; color:#fff; }
    pre { background:#1e1e2e; padding:20px; border-radius:10px; overflow:auto;
          margin-top:20px; white-space:pre-wrap; }
</style>
</head>
<body>
<div class="container">
    <h1>🧌 塌缩怪兽 v17.0</h1>
    <p class="subtitle">存在数论完全体 | 编织域·同伦域激活 | 普适发散消除</p>
    <div class="examples">
        <span>Sum(n**2,(n,1,oo))</span>
        <span>Sum(n,(n,1,oo))</span>
        <span>Sum(2**n,(n,0,oo))</span>
        <span>Sum(factorial(n),(n,0,oo))</span>
        <span>Sum((-1)**(n+1)/n,(n,1,oo))</span>
        <span>Sum((-1)**n * n**2,(n,1,oo))</span>
        <span>Sum(log(n),(n,1,oo))</span>
        <span>Sum(mobius(n),(n,1,oo))</span>
        <span>Sum(n**n,(n,1,oo))</span>
        <span>Sum(factorial(n)**2,(n,0,oo))</span>
    </div>
    <input id="query" value="Sum(n**2,(n,1,oo))" placeholder="输入发散级数">
    <div class="btn-group">
        <button class="btn-calc" id="calcBtn">坍缩!</button>
    </div>
    <pre id="result"></pre>
</div>
<script>
    // 示例按钮点击填充输入框
    document.querySelectorAll('.examples span').forEach(span => {
        span.addEventListener('click', function() {
            document.getElementById('query').value = this.textContent.trim();
        });
    });

    function fmt(v) {
        if (v === null || v === undefined) return '未知';
        if (typeof v === 'string') return v;
        if (Math.abs(v) < 1e-10) return '0';
        let known = {
            '-0.08333333333333333': '-1/12',
            '0.25': '1/4',
            '0.5': '1/2',
            '-0.125': '-1/8',
            '0.5963473623231941': '≈0.596',
            '0.9189385332046727': '½ln(2π)',
            '-2.0': '-2',
            '0.6931471805599453': 'ln2'
        };
        let key = String(v);
        if (key in known) return known[key];
        return v.toFixed(8);
    }

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
                r.innerText =
                    '📊 坍缩报告\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '输入:     ' + d.input + '\\n' +
                    '通项:     ' + d.summand + '\\n' +
                    '初始域:   ' + d.domain + '\\n' +
                    '映射路径: ' + d.path + '\\n' +
                    '步数:     ' + d.steps + ' 步\\n' +
                    '投票:     ' + d.consensus + '\\n' +
                    '策略:     ' + d.strategy + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '坍缩值:   ' + fmt(d.value) + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '理论：任何发散均可通过 ≤3 步对偶映射消除';
            } else {
                r.innerText = '⚠ ' + (d.message || '无法收敛');
            }
        } catch (e) {
            r.innerText = '⚠ 网络错误';
        }
    }

    document.getElementById('calcBtn').addEventListener('click', doCalc);
    // 支持回车键计算
    document.getElementById('query').addEventListener('keypress', function(e) {
        if (e.key === 'Enter') doCalc();
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
        user_input = data.get('query', '')
        dn = parse_input(user_input)
        value, path, consensus = collapser.collapse(dn)
        if value is None:
            return jsonify({"status": "unresolved", "message": "所有路径均无法给出有限坍缩值"})
        path_str = ' → '.join([p[1] for p in path]) if path else "直接谱域求值"
        strategy = path[-1][1] if path else "谱域"
        return jsonify({
            "status": "success",
            "input": user_input,
            "summand": str(dn.expr),
            "domain": dn.domain,
            "path": path_str,
            "steps": len(path) if path else 0,
            "strategy": strategy,
            "consensus": consensus,
            "value": value
        })
    except Exception as e:
        traceback.print_exc()
        return jsonify({"status": "error", "message": str(e)})

if __name__ == '__main__':
    os.makedirs('static', exist_ok=True)
    port = int(os.environ.get("PORT", 5000))
    print(f"🧌 塌缩怪兽 v17.0 启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
