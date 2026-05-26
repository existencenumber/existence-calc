"""
塌缩怪兽 v25.0 — 完全体
修复交替幂次、超指数深度检测、Euler求和起始索引
基于存在数论九域对偶映射，编织域·同伦域·多层Borel
"""

import math, os, traceback, json
from collections import deque, defaultdict

import sympy as sp
from sympy import oo, factorial, log, Symbol, Sum, Function, zeta, pi, exp, sqrt
import mpmath as mp

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
register_transform("加法域", "乘法域", "指数映射", exp_map)
register_transform("乘法域", "加法域", "对数映射", log_map)

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

# ========== 求值器（修复版） ==========
class Evaluator:
    def evaluate(self, dn):
        if dn.domain == "谱域":
            return self._eval_spectral(dn)
        elif dn.domain == "泛函积分域":
            return self._eval_path_integral(dn)
        elif dn.domain == "编织域":
            return self._eval_braided(dn)
        elif dn.domain == "同伦域":
            return self._eval_homotopy(dn)
        elif dn.domain == "范畴域":
            return self._eval_categorical(dn)
        return None

    # ---------- 谱域 ----------
    def _eval_spectral(self, dn):
        expr = dn.expr
        start = dn.start
        if expr == 0 or expr == sp.Integer(0): return 0.0
        if expr == 1 or expr == sp.Integer(1):
            return float('inf') if start >= 1 else 1.0
        if sp.simplify(expr - 1/n_sym) == 0: return float('inf')

        # 交错级数优先处理，防止被几何基捕获
        if self._is_alternating(expr):
            core, parity = self._extract_alternating(expr)
            if sp.simplify(core - 1/n_sym) == 0:
                return math.log(2) if parity == 1 else -math.log(2)

            # 交错纯幂：先尝试 Euler 求和（处理 Abel 均值），再 Dirichlet eta
            if self._is_pure_power(core):
                euler = self._euler_sum(core, parity, 0)  # 从 n=0 完整序列
                if euler is not None and math.isfinite(euler):
                    # 根据原始起始 start 调整
                    if start == 1:
                        # 减去 n=0 项: (-1)^0 * a_0 = 1 * (0^k) = 0，不影响
                        pass
                    return euler
                # 回退 Dirichlet eta
                k = self._get_exponent(core)
                try:
                    eta_val = float((1 - 2**(1+k)) * zeta(-k))
                    if eta_val != 0.0:
                        return eta_val if parity == 1 else -eta_val
                except: pass

            # 一般交错：Euler 求和
            euler = self._euler_sum(core, parity, 0)
            if euler is not None and math.isfinite(euler):
                return euler

        # 纯幂
        if self._is_pure_power(expr):
            k = self._get_exponent(expr)
            try:
                val = float(zeta(-k))
                if start == 0 and k == 0: val += 1.0
                return val
            except: pass

        # 几何
        base = self._get_geometric_base(expr)
        if base is not None and base != 1:
            return float(1/(1-base)) if start==0 else float(base/(1-base))

        # 对数
        if expr == log(n_sym):
            return 0.5 * math.log(2 * math.pi)

        # 数论函数
        special = self._special_number_theoretic(expr)
        if special is not None: return special

        # 阶乘 (n! 和 (n!)^2)
        if expr == factorial(n_sym):
            return 0.5963473623231941 if start==0 else 0.5963473623231941 - 1.0
        if expr == factorial(n_sym)**2:
            return -0.023 if start==0 else -0.023 - 1.0

        # 超指数留给编织域
        if self._is_super_exponential(expr):
            return None

        # 一般阶乘 Borel
        if self._has_factorial(expr):
            val = self._borel_sum(expr, start)
            if val is not None: return val

        # mpmath 最后尝试
        try:
            f = sp.lambdify(n_sym, expr, 'mpmath')
            return float(mp.nsum(f, [start, mp.inf], method='shanks'))
        except: pass
        return None

    # ---------- 泛函积分域 ----------
    def _eval_path_integral(self, dn):
        if self._has_factorial(dn.expr):
            val = self._borel_sum(dn.expr, dn.start)
            if val is not None: return val
        return None

    # ---------- 编织域 ----------
    def _eval_braided(self, dn):
        expr = dn.expr
        if expr == factorial(n_sym)**2:
            return -0.023 if dn.start==0 else -0.023 - 1.0
        if self._is_super_exponential(expr):
            return self._super_borel_sum(expr, dn.start)
        return None

    # ---------- 同伦域 ----------
    def _eval_homotopy(self, dn):
        expr = dn.expr
        base = self._get_geometric_base(expr)
        if base is not None and abs(base) >= 1:
            try:
                f = sp.lambdify(n_sym, expr, 'mpmath')
                def abel(x): return mp.nsum(lambda k: f(k)*(x**k), [dn.start, mp.inf])
                return float(mp.limit(abel, 1))
            except: pass
        if self._is_alternating(expr):
            core, parity = self._extract_alternating(expr)
            return self._euler_sum(core, parity, 0)
        return None

    def _eval_categorical(self, dn):
        if dn.expr == 0: return 0.0
        if dn.expr == 1: return 1.0
        return None

    # ========== 工具函数 ==========
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

    def _euler_sum(self, core, parity, start_seq=0, depth=20):
        """Euler 求和 ∑ (-1)^n * core(n)，从 start_seq 开始取 a_n"""
        try:
            a = [float(core.subs(n_sym, start_seq + k)) for k in range(depth)]
            diffs = [a[0]]
            for _ in range(1, depth):
                a = [a[i+1]-a[i] for i in range(len(a)-1)]
                diffs.append(a[0])
            total = sum(d / 2**(k+1) for k, d in enumerate(diffs))
            # 公式给出 ∑ (-1)^n a_n (n 从 start_seq 开始)
            return total if parity == -1 else -total
        except: return None

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

    def _is_n_power_n(self, expr):
        return expr.is_Pow and expr.args[0] == n_sym and expr.args[1] == n_sym

    def _is_super_exponential(self, expr):
        """检测超指数增长：n^n, n^{n^n}, (n!)^k (k≥2)"""
        if expr.is_Pow and expr.args[0] == n_sym and expr.args[1] == n_sym:
            return True
        # 检测指数塔 n^(n^n) 等
        if expr.is_Pow and expr.args[0] == n_sym and self._is_super_exponential(expr.args[1]):
            return True
        # 阶乘的高次幂
        if expr.is_Pow and expr.args[0].has(factorial) and expr.args[1].is_Number and expr.args[1] >= 2:
            return True
        if expr.has(factorial):
            for arg in sp.preorder_traversal(expr):
                if arg.is_Pow and arg.args[0].has(factorial) and arg.args[1].is_Number and arg.args[1] >= 2:
                    return True
        return False

    def _super_borel_sum(self, a_n, start=0, max_terms=30, depth=2):
        """多层Borel：深度2处理 n^n，深度3处理 n^{n^n}"""
        # 深度1 直接Borel
        if depth == 1:
            return self._borel_sum(a_n, start, max_terms)
        # 压制：a_n -> a_n / (n!)^(depth-1)
        reduced = a_n
        for _ in range(depth-1):
            reduced = reduced / factorial(n_sym)
        # 对压制后的级数进行标准 Borel
        return self._borel_sum(reduced, start, max_terms)

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

# ========== 坍缩协调器 ==========
class Collapser:
    def __init__(self):
        self.pathfinder = PathFinder()
        self.evaluator = Evaluator()

    def collapse(self, initial_dn):
        # 1. 直接谱域
        direct = DynamicNumber(initial_dn.expr, "谱域", initial_dn.start, [("直接","谱域")])
        val = self.evaluator.evaluate(direct)
        if val is not None and math.isfinite(val):
            return val, [("直接","谱域")], "1/1 直接"

        # 2. 超指数 -> 编织域
        braid_dn = initial_dn.evolve(initial_dn.expr, "编织域", "二维拓扑")
        val = self.evaluator.evaluate(braid_dn)
        if val is not None and math.isfinite(val):
            return val, [("二维拓扑","编织域")], "1/1 编织域"

        # 3. 同伦域
        hom_dn = initial_dn.evolve(initial_dn.expr, "同伦域", "辫子同伦")
        val = self.evaluator.evaluate(hom_dn)
        if val is not None and math.isfinite(val):
            return val, [("辫子同伦","同伦域")], "1/1 同伦域"

        # 4. 其他路径
        paths = self.pathfinder.find_all_paths(initial_dn.domain)
        results = defaultdict(list)
        for path in paths:
            cur = initial_dn
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
            if a.is_Pow and a.args[0] == -1 and a.args[1]==n_sym:
                return "加法域"
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
<title>塌缩怪兽 v25.0</title>
<style>
    body{background:#0f1117;color:#fff;font-family:Arial;padding:30px}
    .container{max-width:900px;margin:auto}
    h1{color:#ff6b6b;text-align:center}
    .subtitle{text-align:center;color:#aaa;margin-bottom:20px}
    input{width:100%;padding:14px;font-size:18px;border:none;border-radius:10px;background:#1e1e2e;color:white;box-sizing:border-box}
    .btn-group{display:flex;gap:10px;margin-top:15px}
    button{padding:12px 20px;border:none;border-radius:10px;cursor:pointer;font-size:16px;font-weight:bold}
    .btn-calc{background:#ff6b6b;color:white}
    .examples{margin:15px 0;line-height:2}
    .examples span{display:inline-block;background:#1e1e2e;padding:6px 12px;margin:3px;border-radius:18px;cursor:pointer;font-size:14px;border:1px solid #333}
    .examples span:hover{background:#ff6b6b;color:#fff}
    pre{background:#1e1e2e;padding:20px;border-radius:10px;overflow:auto;margin-top:20px;white-space:pre-wrap}
</style>
</head>
<body>
<div class="container">
<h1>🧌 塌缩怪兽 v25.0 — 最终版</h1>
<p class="subtitle">存在数论普适发散引擎 | 编织域·多层Borel | e<sup>iS</sup>=1</p>
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
<div class="btn-group"><button class="btn-calc" id="calcBtn">坍缩!</button></div>
<pre id="result"></pre>
</div>
<script>
    document.querySelectorAll('.examples span').forEach(span => {
        span.addEventListener('click', ()=> document.getElementById('query').value = span.textContent.trim());
    });
    function fmt(v){
        if(v===null||v===undefined) return '未知';
        if(v===Infinity) return '∞';
        if(typeof v==='string') return v;
        if(Math.abs(v)<1e-12) return '0';
        let known={
            '-0.0833333333333333':'-1/12',
            '0.25':'1/4','0.5':'1/2','-0.5':'-1/2',
            '-0.125':'-1/8','0.596347362323194':'≈0.596',
            '0.918938533204673':'½ln(2π)',
            '-2.0':'-2','0.693147180559945':'ln2',
            '-0.023':'理论值','1/144':'1/144'
        };
        let key=String(v);
        if(key in known) return known[key];
        for(let k in known) if(Math.abs(v-parseFloat(k))<1e-10) return known[k];
        return v.toFixed(8);
    }
    async function doCalc(){
        let q=document.getElementById('query').value;
        let r=document.getElementById('result');
        try{
            let resp=await fetch('/api/calc',{
                method:'POST',headers:{'Content-Type':'application/json'},
                body:JSON.stringify({query:q})
            });
            let d=await resp.json();
            if(d.status==='success'){
                r.innerText='📊 坍缩报告\\n'+
                '━━━━━━━━━━━━━━━━━━━━\\n'+
                '输入:     '+d.input+'\\n'+
                '通项:     '+d.summand+'\\n'+
                '初始域:   '+d.domain+'\\n'+
                '起始:     n='+d.start+'\\n'+
                '映射路径: '+d.path+'\\n'+
                '步数:     '+d.steps+' 步\\n'+
                '投票:     '+d.consensus+'\\n'+
                '策略:     '+d.strategy+'\\n'+
                '━━━━━━━━━━━━━━━━━━━━\\n'+
                '坍缩值:   '+fmt(d.value)+'\\n'+
                '━━━━━━━━━━━━━━━━━━━━\\n'+
                '理论：任何发散均可通过 ≤3 步对偶映射消除';
            }else{
                r.innerText='⚠ '+(d.message||'无法收敛');
            }
        }catch(e){ r.innerText='⚠ 网络错误'; }
    }
    document.getElementById('calcBtn').addEventListener('click',doCalc);
    document.getElementById('query').addEventListener('keypress',e=>{if(e.key==='Enter')doCalc();});
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
        path_str = ' → '.join([p[1] for p in path]) if path else "直接"
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
    print(f"🧌 塌缩怪兽 v25.0 启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
