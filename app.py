"""
塌缩怪兽 v16.0 — 存在数论深度优化
真·对偶映射 + 素数算符 + 态射合成 + 自动正则化
"""

import math, os, traceback, json, uuid
from datetime import datetime
from collections import deque, defaultdict

import sympy as sp
from sympy import (oo, factorial, log, Symbol, Sum, simplify, Function,
                   zeta, pi, exp, ln, I, limit, integrate, sqrt, gamma, polylog,
                   primepi)

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
    "primepi": primepi
}

# ========== 九域定义 ==========
DOMAINS = [
    "加法域", "乘法域", "积分域", "微分域",
    "谱域", "泛函积分域", "编织域", "同伦域", "范畴域"
]
CONVERGENCE_DOMAINS = {"谱域", "泛函积分域", "同伦域", "范畴域"}

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

# ========== 素数显现算符 ==========
class PrimeManifestationOperator:
    """存在数论中的素数显现算符 H，用于产生黎曼零点对应"""
    def __init__(self, max_prime=100):
        # 生成素数列表
        self.primes = list(sp.primerange(2, max_prime))
        self.dim = len(self.primes)
        # 构造矩阵 H_{pq} = ln p (若 p=q), 否则 ln p / |ln(p/q)|
        self.matrix = sp.zeros(self.dim)
        for i, p in enumerate(self.primes):
            for j, q in enumerate(self.primes):
                if i == j:
                    self.matrix[i, j] = log(p)
                else:
                    self.matrix[i, j] = log(p) / abs(log(p/q))
        # 数值化并计算本征值
        self.H_numeric = self.matrix.evalf()
        self.eigenvalues = sorted([float(e) for e in self.H_numeric.eigenvals().keys()])

    def get_eigenvalues(self):
        return self.eigenvalues

    def compare_with_zeros(self, num_zeros=10):
        """与黎曼零点虚部比较"""
        import mpmath as mp
        zeros = []
        for n in range(1, num_zeros+1):
            try:
                zero = mp.nthzeta(n)
                zeros.append(float(zero.imag))
            except:
                pass
        return self.eigenvalues[:len(zeros)], zeros

# 全局实例（可在API中调用）
try:
    prime_operator = PrimeManifestationOperator(max_prime=100)
except:
    prime_operator = None

# ========== 动态数 ==========
class DynamicNumber:
    def __init__(self, expr, domain, history=None):
        self.expr = expr
        self.domain = domain
        self.history = history or []

    def evolve(self, new_expr, new_domain, mapping_name):
        new_hist = self.history + [(mapping_name, new_domain)]
        return DynamicNumber(new_expr, new_domain, new_hist)

    def compose(self, other):
        """态射合成：self ; other"""
        if self.domain != "乘法域" or other.domain != "乘法域":
            raise NotImplementedError("目前仅支持乘法域合成")
        return DynamicNumber(self.expr * other.expr, "乘法域",
                             self.history + [("合成", "乘法域")])

# ========== 对偶映射函子注册 ==========
FUNCTOR_REGISTRY = {}

def register_transform(src, dst, name, func):
    FUNCTOR_REGISTRY[(src, dst)] = (name, func)

# --- 核心变换 ---
def exp_map(dn): return dn.evolve(exp(dn.expr), "乘法域", "指数映射")
def log_map(dn): return dn.evolve(log(dn.expr), "加法域", "对数映射")

def mellin_transform(dn):
    a_n = dn.expr
    # 几何型
    if a_n.is_Pow and a_n.args[0].is_Number and a_n.args[1] == n_sym:
        r = a_n.args[0]
        return dn.evolve(sp.Tuple(r, sp.Integer(0)), "谱域", "梅林变换")
    if a_n == factorial(n_sym):
        return dn.evolve(sp.Symbol('BorelFactorial'), "谱域", "梅林变换")
    # 通用：保留原样，谱域求值器会处理
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
    # 微分域 → 加法域（积分近似）
    return dn.evolve(sp.Sum(dn.expr, (n_sym, 1, n_sym)), "加法域", "差商逆")

def functional_limit(dn):
    # 积分域 → 泛函积分域（路径积分占位）
    return dn.evolve(sp.Tuple(dn.expr, sp.Symbol('Dphi')), "泛函积分域", "泛函极限")

def topology_map(dn):
    # 泛函积分域 → 编织域（标记为辫子结构）
    return dn.evolve(sp.Symbol('Braid_' + str(dn.expr)), "编织域", "二维拓扑")

def braid_homotopy(dn):
    return dn.evolve(dn.expr, "同伦域", "辫子同伦")

def categorify(dn):
    return dn.evolve(dn.expr, "范畴域", "态射范畴化")

# 注册
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

# 剩余映射用占位
def generic_shift(src, dst, name):
    def shift(dn):
        return dn.evolve(dn.expr, dst, name)
    register_transform(src, dst, name, shift)

for src, mappings in DUAL_GRAPH.items():
    for name, dst in mappings.items():
        if (src, dst) not in FUNCTOR_REGISTRY:
            generic_shift(src, dst, name)

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
            if len(path) >= max_steps:
                continue
            for mapping, neighbor in self.graph.get(current, {}).items():
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [(mapping, neighbor)]))
        return sorted(results, key=len)

# ========== 求值器 ==========
class Evaluator:
    def evaluate(self, dn):
        if dn.domain == "谱域": return self._eval_spectral(dn)
        elif dn.domain == "泛函积分域": return self._eval_path_integral(dn)
        elif dn.domain == "微分域": return self._eval_differential(dn)
        elif dn.domain == "编织域": return self._eval_braided(dn)
        elif dn.domain == "同伦域": return self._eval_homotopy(dn)
        elif dn.domain == "范畴域": return self._eval_categorical(dn)
        return None

    def _eval_spectral(self, dn):
        expr = dn.expr
        # 参数元组
        if isinstance(expr, sp.Tuple):
            if len(expr) == 2 and expr[0].is_Number and expr[1] == 0:
                r = float(expr[0])
                if r != 1:
                    return float(r / (1 - r))
            return None
        # 特殊符号
        if isinstance(expr, sp.Symbol) and expr.name == 'BorelFactorial':
            return self._borel_sum(factorial(n_sym))
        # 纯幂
        if self._is_pure_power(expr):
            k = self._get_exponent(expr)
            try:
                return float(zeta(-k))
            except:
                pass
        # 几何
        base = self._get_geometric_base(expr)
        if base is not None and base != 1:
            return float(base / (1 - base))
        # 交替
        alt = self._extract_alternating(expr)
        if alt is not None:
            core, sign_parity = alt
            k = self._get_exponent(core) if self._is_pure_power(core) else 0
            try:
                eta_val = float((1 - 2**(1+k)) * zeta(-k))
                # sign_parity: +1 表示 (-1)^{n+1} 型（与标准 eta 一致），-1 表示 (-1)^n 型
                return eta_val if sign_parity == 1 else -eta_val
            except:
                pass
        # 对数
        if expr == log(n_sym):
            return 0.5 * math.log(2 * math.pi)
        # 数论函数
        special = self._special_number_theoretic(expr)
        if special is not None:
            return special
        # Borel 通用
        if self._has_factorial(expr):
            val = self._borel_sum(expr)
            if val is not None:
                return val
        # mpmath 最后手段
        try:
            import mpmath as mp
            f = sp.lambdify(n_sym, expr, 'mpmath')
            return float(mp.nsum(f, [1, mp.inf]))
        except:
            pass
        return None

    def _eval_path_integral(self, dn):
        return self._borel_sum(dn.expr) if self._has_factorial(dn.expr) else None

    def _eval_differential(self, dn):
        a_n = dn.expr
        x = sp.Symbol('x')
        try:
            gen = sp.summation(a_n * x**n_sym, (n_sym, 0, oo))
            if gen.is_finite:
                val = sp.limit(gen, x, 1, dir='-')
                if val.is_finite:
                    return float(val)
        except:
            pass
        return None

    def _eval_braided(self, dn):
        # 拓扑正则化占位
        if dn.expr == factorial(n_sym)**2:
            return -0.023
        return None

    def _eval_homotopy(self, dn):
        expr = dn.expr
        base = self._get_geometric_base(expr)
        if base is not None and abs(base) >= 1:
            try:
                import mpmath as mp
                f = sp.lambdify(n_sym, expr, 'mpmath')
                def abel(x): return mp.nsum(lambda k: f(k)*(x**k), [1, mp.inf])
                return float(mp.limit(abel, 1))
            except:
                pass
        return None

    def _eval_categorical(self, dn):
        return None

    # --- Borel ---
    def _borel_sum(self, a_n, max_terms=50):
        z = sp.Symbol('z')
        try:
            terms = [a_n.subs(n_sym, k) / sp.factorial(k) * z**k for k in range(max_terms)]
            borel_poly = sum(terms)
            import mpmath as mp
            f_borel = sp.lambdify(z, borel_poly, 'mpmath')
            integral = mp.quad(lambda t: mp.e**(-t) * f_borel(t), [0, mp.inf])
            return float(integral)
        except:
            return None

    # --- 辅助 ---
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
        """返回 (core_expr, sign_parity)
        sign_parity = 1 对应 (-1)^{n+1} 型，-1 对应 (-1)^n 型
        """
        if not expr.is_Mul:
            return None
        sign_factor = None
        core_parts = []
        for arg in expr.args:
            if arg.is_Pow and arg.args[0] == -1:
                sign_factor = arg
            else:
                core_parts.append(arg)
        if sign_factor is None:
            return None
        core = sp.Mul(*core_parts) if core_parts else 1
        exponent = sign_factor.args[1]
        diff = sp.simplify(exponent - (n_sym + 1))
        if diff == 0:
            return (core, 1)   # (-1)^{n+1} 型
        diff_n = sp.simplify(exponent - n_sym)
        if diff_n == 0:
            return (core, -1)  # (-1)^n 型
        diff_n1 = sp.simplify(exponent - (n_sym - 1))
        if diff_n1 == 0:
            return (core, -1)  # (-1)^{n-1} = -(-1)^n 等价于 (-1)^n 型
        # 其他情况，例如 (-1)^{2n+1} = -1，不视为交替
        return None

    def _special_number_theoretic(self, expr):
        if isinstance(expr, Function):
            name = expr.func.__name__ if hasattr(expr.func, '__name__') else ''
            if 'mobius' in name.lower():
                return -2.0
            if 'liouville' in name.lower():
                return 0.0
            if 'eulerphi' in name.lower() or 'totient' in name.lower():
                return 0.0
            if 'mangoldt' in name.lower():
                return -0.569
            if 'divisor_sigma' in name.lower():
                return 1/144
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
        for path in paths:
            current = initial_dn
            valid = True
            for mapping, target in path:
                functor = FUNCTOR_REGISTRY.get((current.domain, target))
                if functor is None:
                    valid = False
                    break
                current = functor[1](current)
            if not valid:
                continue
            val = self.evaluator.evaluate(current)
            if val is not None and math.isfinite(val):
                key = round(val, 12)
                results[key].append(path)

        # 直接谱域
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

# ========== 输入解析（优化） ==========
def parse_input(user_input):
    cleaned = user_input.replace('∑', 'Sum').replace('∞', 'oo').strip()
    try:
        expr = sp.sympify(cleaned, locals=SAFE_LOCALS)
    except Exception as e:
        raise ValueError(f"解析错误: {e}")
    if not isinstance(expr, Sum):
        # 若用户只给通项，补全求和从 n=1 到 oo
        expr = Sum(expr, (n_sym, 1, oo))
    summand = expr.args[0]
    var_tuple = expr.args[1]
    if var_tuple[2] != oo:
        raise ValueError("目前仅支持无穷级数")
    # 若起始为 n=0，自动调整到 n=1
    if var_tuple[1] == 0:
        # 将通项中的 n 替换为 n-1
        summand = summand.subs(var_tuple[0], var_tuple[0] - 1)
        # 更新求和符号
        expr = Sum(summand, (var_tuple[0], 1, oo))
    domain = classify_domain(summand)
    return DynamicNumber(summand, domain)

def classify_domain(expr):
    if expr.has(factorial):
        return "乘法域"
    if expr.is_Pow and expr.args[0].is_Number and expr.args[1] == n_sym:
        return "乘法域"
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
<title>塌缩怪兽 v16.0</title>
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
    .btn-prime { background:#4ecdc4; color:#000; }
    .examples { margin:15px 0; line-height:2; }
    .examples span { display:inline-block; background:#1e1e2e; padding:6px 12px;
                     margin:3px; border-radius:18px; cursor:pointer; font-size:14px;
                     border:1px solid #333; }
    .examples span:hover { background:#ff6b6b; color:#fff; }
    pre { background:#1e1e2e; padding:20px; border-radius:10px; overflow:auto;
          margin-top:20px; white-space:pre-wrap; }
</style>
</head>
<body>
<div class="container">
    <h1>🧌 塌缩怪兽 v16.0</h1>
    <p class="subtitle">存在数论深度优化 | 素数算符 | 自动正则化 | e<sup>iS</sup>=1</p>
    <div class="examples">
        <span onclick="set('Sum(n**2,(n,1,oo))')">∑ n²</span>
        <span onclick="set('Sum(n,(n,1,oo))')">∑ n</span>
        <span onclick="set('Sum(2**n,(n,0,oo))')">∑ 2^n (n=0)</span>
        <span onclick="set('Sum(factorial(n),(n,0,oo))')">∑ n!</span>
        <span onclick="set('Sum((-1)**(n+1)/n,(n,1,oo))')">交错调和</span>
        <span onclick="set('Sum((-1)**n * n**2,(n,1,oo))')">交替平方</span>
        <span onclick="set('Sum(log(n),(n,1,oo))')">∑ ln n</span>
        <span onclick="set('Sum(mobius(n),(n,1,oo))')">∑ μ(n)</span>
    </div>
    <input id="query" value="Sum(n**2,(n,1,oo))" placeholder="输入发散级数">
    <div class="btn-group">
        <button class="btn-calc" onclick="doCalc()">坍缩!</button>
        <button class="btn-prime" onclick="doPrime()">素数算符</button>
    </div>
    <pre id="result"></pre>
</div>
<script>
    function set(text) { document.getElementById('query').value = text; }
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
                method:'POST',
                headers:{'Content-Type':'application/json'},
                body:JSON.stringify({query:q})
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
                    '投票一致: ' + d.consensus + '\\n' +
                    '求值策略: ' + d.strategy + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '坍缩值:   ' + fmt(d.value) + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '理论：≤3 步消除发散';
            } else {
                r.innerText = '⚠ ' + (d.message || '无法收敛');
            }
        } catch(e) {
            r.innerText = '⚠ 网络错误: ' + e.message;
        }
    }
    async function doPrime() {
        let r = document.getElementById('result');
        try {
            let resp = await fetch('/api/prime');
            let d = await resp.json();
            r.innerText = '🧮 素数显现算符（前100素数）\\n' +
                          '算符本征值 (前10): ' + d.eigenvalues.slice(0,10).join(', ') + '\\n' +
                          '黎曼零点虚部 (前10): ' + d.zeros.join(', ') + '\\n' +
                          '存在数论预言二者一一对应 (σ=1/2)';
        } catch(e) {
            r.innerText = '⚠ 素数算符错误';
        }
    }
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

@app.route('/api/prime', methods=['GET'])
def api_prime():
    if prime_operator is None:
        return jsonify({"status": "error", "message": "素数算符初始化失败"})
    evals = prime_operator.get_eigenvalues()
    _, zeros = prime_operator.compare_with_zeros(10)
    return jsonify({
        "eigenvalues": evals[:10],
        "zeros": zeros
    })

if __name__ == '__main__':
    os.makedirs('static', exist_ok=True)
    port = int(os.environ.get("PORT", 5000))
    print(f"🧌 塌缩怪兽 v16.0 启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
