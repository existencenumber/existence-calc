"""
塌缩怪兽 v11.5 — 稳健交替级数检测
基于存在数论的非微扰计算工具
"""

import re, math, os, json, traceback, uuid
from datetime import datetime

import sympy as sp
from sympy import Sum, oo, factorial, log, Symbol, simplify, Mul, Pow, Integer

SAFE_LOCALS = {
    "Sum": Sum, "oo": oo, "factorial": factorial,
    "log": sp.log, "sin": sp.sin, "cos": sp.cos,
    "exp": sp.exp, "sqrt": sp.sqrt, "n": Symbol("n")
}

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from flask import Flask, request, jsonify, render_template_string

DUAL_GRAPH = {
    "加法域": {"指数映射": "乘法域", "黎曼和极限": "积分域", "差商极限": "微分域"},
    "乘法域": {"对数映射": "加法域", "对数导数": "积分域", "梅林变换": "谱域"},
    "积分域": {"拉普拉斯变换": "谱域", "泛函极限": "泛函积分域", "积分的逆": "微分域"},
    "微分域": {"傅里叶变换": "谱域", "微分的逆": "积分域"},
    "谱域": {"逆傅里叶变换": "微分域", "逆拉普拉斯变换": "积分域", "逆梅林变换": "乘法域"},
    "泛函积分域": {"泛函极限的逆": "积分域", "二维拓扑": "编织域"},
    "编织域": {"辫子同伦": "同伦域"},
    "同伦域": {"态射范畴化": "范畴域"},
    "范畴域": {"恒等态射对应0": "加法域", "恒等态射对应1": "乘法域"},
}
CONVERGENCE_DOMAINS = {"谱域", "泛函积分域", "同伦域", "范畴域"}

DOMAIN_POSITIONS = {
    "加法域": (0, 2), "乘法域": (2, 2), "谱域": (1, 0),
    "积分域": (0, 0), "微分域": (2, 0), "泛函积分域": (3, 1),
}
DOMAIN_COLORS = {
    "加法域": "#ff6b6b", "乘法域": "#4ecdc4", "谱域": "#ffe66d",
    "积分域": "#96ceb4", "微分域": "#45b7d1", "泛函积分域": "#dda0dd",
}

def compute_zeta(s):
    table = {0:-0.5, -1:-1/12, -2:0, -3:1/120, -4:0, -5:-1/252, -6:0, -7:-1/240, -8:0,
             1:float('inf'), 2:math.pi**2/6, 3:1.2020569031595942}
    if s in table: return table[s]
    if s < 0 and s % 2 == 0: return 0
    try: return float(sp.N(sp.zeta(s)))
    except: return None

def compute_eta(s):
    """eta(s) 或 eta((s, sign))"""
    sign = 1
    s_val = s
    if isinstance(s, tuple):
        s_val, sign = s
    eta_vals = {0:0.5, -1:0.25, -2:0, -3:-1/8, 1:math.log(2)}
    if s_val in eta_vals:
        return sign * eta_vals[s_val]
    z = compute_zeta(s_val)
    if z is None or z == float('inf'): return None
    return sign * (1 - 2**(1 - s_val)) * z

def compute_abel(r):
    if r == 1: return float('inf')
    try: return 1/(1-r)
    except: return None

def compute_euler(r): return 1/(1+abs(r))
def compute_borel(): return 0.5963473623231941
def compute_fibonacci(): return -1
def compute_mobius(): return -2
def compute_liouville(): return 0
def compute_euler_phi(): return 0
def compute_mangoldt(): return -0.569
def compute_divisor(): return 1/144
def compute_superfactorial(): return -0.082
def compute_primorial(): return -0.064
def compute_euler_char(): return -5/6
def compute_ramanujan_cf(): return (math.sqrt(5)+1)/2
def compute_gr_1loop(): return 0
def compute_gr_2loop(): return 0
def compute_prime_counting(): return -0.128
def compute_prime_zeta(s):
    m = {-1:-0.0942, -2:-0.023, -3:-0.008, -4:-0.003}
    return m.get(s, 0)
def compute_zeta_deriv(s):
    if s == 0: return 0.5*math.log(2*math.pi)
    return None

def analyze_term(expr, var='n'):
    sym = Symbol(var)
    expr_s = sp.simplify(expr)

    # ── 交替因子检测 ──
    # 寻找 (-1)**(a*n+b) 形式的因子
    alternating_sign = 0  # 0:未发现, 1:有交替, 指数中的n系数
    alternating_const = 0
    rest_expr = expr_s
    if isinstance(expr_s, Mul):
        factors = expr_s.args
    else:
        factors = [expr_s]
    new_factors = []
    for f in factors:
        if isinstance(f, Pow) and sp.simplify(f.base) == -1:
            exp = f.exp
            if exp.has(sym):
                # 提取系数
                coeff = exp.coeff(sym)
                const_term = sp.simplify(exp - coeff * sym)
                alternating_sign = int(coeff) if coeff.is_Integer else 0
                alternating_const = int(const_term) if const_term.is_Integer else 0
                continue  # 去掉这个因子
        new_factors.append(f)
    if alternating_sign != 0:
        rest_expr = Mul(*new_factors) if new_factors else 1
        rest_expr = sp.simplify(rest_expr)
        # 判断相对于标准 eta 定义 ∑ (-1)^(n-1) 的符号
        # 我们的因子是 (-1)^(k*n + b)
        # 标准: (-1)^(n-1) => k=1, b=-1
        # 对于 (-1)^(n+1): k=1, b=1, 符号 = (-1)^{(b+1)} = (-1)^2 = 1
        # 对于 (-1)^n: k=1, b=0, 符号 = (-1)^{(0+1)} = -1
        if alternating_sign == 1:
            sign = (-1)**(alternating_const + 1)
        elif alternating_sign == -1:
            sign = (-1)**(-alternating_const + 1)
        else:
            sign = 1
        if rest_expr.is_polynomial(sym):
            deg = sp.degree(rest_expr, gen=sym)
            return {"type": "alternating_power", "domain": "加法域",
                    "mapping": ["加法域", "谱域"],
                    "method": "eta", "param": (-deg, sign)}

    # 多项式 n^k
    if expr_s.is_polynomial(sym):
        deg = sp.degree(expr_s, gen=sym)
        if deg >= 0:
            return {"type": "polynomial", "domain": "加法域",
                    "mapping": ["加法域", "乘法域", "谱域"],
                    "method": "zeta", "param": -deg}

    # 几何 r^n
    if expr_s.is_Pow and expr_s.exp.has(sym):
        base = float(expr_s.base)
        return {"type": "geometric", "domain": "乘法域",
                "mapping": ["乘法域", "谱域"],
                "method": "abel", "param": base}
    for atom in expr_s.atoms():
        if atom.is_Pow and atom.exp == sym:
            base = float(atom.base)
            return {"type": "geometric", "domain": "乘法域",
                    "mapping": ["乘法域", "谱域"],
                    "method": "abel", "param": base}
    for atom in expr_s.atoms():
        if atom.is_Pow and atom.exp.has(sym):
            base = float(atom.base)
            return {"type": "geometric", "domain": "乘法域",
                    "mapping": ["乘法域", "谱域"],
                    "method": "abel", "param": base}

    # 阶乘 n!
    if expr_s.has(sp.factorial):
        return {"type": "factorial", "domain": "乘法域",
                "mapping": ["乘法域", "谱域"],
                "method": "borel", "param": None}

    # 调和 1/n
    if sp.simplify(expr_s - 1/sym) == 0:
        return {"type": "harmonic", "domain": "加法域",
                "mapping": ["加法域", "谱域"],
                "method": "zeta", "param": 1}

    # 对数 ln n
    if expr_s == sp.log(sym):
        return {"type": "logarithmic", "domain": "加法域",
                "mapping": ["加法域", "谱域"],
                "method": "zeta_deriv", "param": 0}

    # 特殊函数
    fn = str(expr_s.func) if hasattr(expr_s, 'func') else ''
    if 'mobius' in fn.lower():
        return {"type": "mobius", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "mobius", "param": None}
    if 'liouville' in fn.lower():
        return {"type": "liouville", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "liouville", "param": None}
    if 'eulerphi' in fn.lower() or 'totient' in fn.lower():
        return {"type": "euler_phi", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "euler_phi", "param": None}
    if 'divisor_sigma' in fn.lower():
        return {"type": "divisor", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "divisor", "param": None}
    if 'fibonacci' in fn.lower():
        return {"type": "fibonacci", "domain": "乘法域",
                "mapping": ["乘法域", "谱域"], "method": "fibonacci", "param": None}

    return {"type": "unknown", "domain": "未知",
            "mapping": ["未知"], "method": None, "param": None}

def compute_spectral(method, param):
    methods = {
        'zeta': compute_zeta, 'abel': compute_abel, 'euler': compute_euler,
        'borel': compute_borel, 'zeta_deriv': compute_zeta_deriv,
        'fibonacci': compute_fibonacci, 'mobius': compute_mobius,
        'liouville': compute_liouville, 'euler_phi': compute_euler_phi,
        'mangoldt': compute_mangoldt, 'divisor': compute_divisor,
        'superfactorial': compute_superfactorial, 'primorial': compute_primorial,
        'euler_char': compute_euler_char, 'ramanujan_cf': compute_ramanujan_cf,
        'gr_1loop': compute_gr_1loop, 'gr_2loop': compute_gr_2loop,
        'prime_counting': compute_prime_counting, 'prime_zeta': compute_prime_zeta,
    }
    if method == 'eta':
        return compute_eta(param)
    if method in methods:
        return methods[method](param) if param is not None else methods[method]()
    return None

def compute(user_input):
    if not user_input or not user_input.strip():
        return {"status": "error", "message": "输入不能为空"}
    try:
        cleaned = user_input.replace('∑', 'Sum').replace('∞', 'oo').strip()
        n = Symbol('n')
        expr = sp.sympify(cleaned, locals=SAFE_LOCALS)
        if not isinstance(expr, Sum):
            expr = Sum(expr, (n, 1, oo))
        summand = expr.args[0]
        var_tuple = expr.args[1]
        upper = var_tuple[2]
        if upper != oo:
            return {"status": "finite", "message": "有限级数可直接求和", "summand": str(summand)}
        analysis = analyze_term(summand)
        value = compute_spectral(analysis["method"], analysis["param"])
        result_value = value
        if value == float('inf'):
            result_value = "∞"
        return {
            "status": "success", "timestamp": datetime.utcnow().isoformat(),
            "input": user_input, "summand": str(summand),
            "divergence_type": analysis["type"], "domain": analysis["domain"],
            "mapping_path": " → ".join(analysis["mapping"]),
            "steps": len(analysis["mapping"]) - 1,
            "method": analysis["method"], "value": result_value
        }
    except Exception as e:
        traceback.print_exc()
        return {"status": "error", "message": str(e)}

def visualize_mapping(path):
    if not HAS_MPL: return None
    filename = f"graph_{uuid.uuid4().hex}.png"
    save_path = os.path.join("static", filename)
    os.makedirs("static", exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(-1, 4); ax.set_ylim(-1, 3); ax.axis("off")
    for domain in path:
        if domain in DOMAIN_POSITIONS:
            x, y = DOMAIN_POSITIONS[domain]
            c = plt.Circle((x, y), 0.35, color=DOMAIN_COLORS.get(domain, "#ccc"),
                           ec="black", linewidth=2, zorder=2)
            ax.add_patch(c)
            ax.text(x, y, domain, ha="center", va="center", fontsize=9, weight="bold")
    for i in range(len(path)-1):
        d1, d2 = path[i], path[i+1]
        if d1 in DOMAIN_POSITIONS and d2 in DOMAIN_POSITIONS:
            x1, y1 = DOMAIN_POSITIONS[d1]; x2, y2 = DOMAIN_POSITIONS[d2]
            ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle='->', color='red', lw=2.5, zorder=3))
    plt.title(" → ".join(path), fontsize=12)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    return f"/static/{filename}"

FEEDBACK_LOG = "feedback.log"

def log_unresolved(user_input, error_msg):
    entry = {"timestamp": datetime.utcnow().isoformat(), "input": user_input, "error": error_msg}
    with open(FEEDBACK_LOG, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

app = Flask(__name__)

HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>塌缩怪兽 v11.5</title>
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
    .btn-viz { background:#4ecdc4; color:#000; }
    .examples { margin:15px 0; line-height:2; }
    .examples span { display:inline-block; background:#1e1e2e; padding:6px 12px;
                     margin:3px; border-radius:18px; cursor:pointer; font-size:14px;
                     border:1px solid #333; }
    .examples span:hover { background:#ff6b6b; color:#fff; }
    pre { background:#1e1e2e; padding:20px; border-radius:10px; overflow:auto;
          margin-top:20px; white-space:pre-wrap; }
    img { width:100%; margin-top:20px; border-radius:10px; }
    .loading { text-align:center; padding:20px; color:#aaa; display:none; }
    .loading.show { display:block; }
</style>
</head>
<body>
<div class="container">
    <h1>🧌 塌缩怪兽 v11.5</h1>
    <p class="subtitle">存在数论非微扰计算器 | 九域对偶映射 | e<sup>iS</sup> = 1</p>
    <div class="examples">
        <span onclick="set('Sum(n**2,(n,1,oo))')">∑ n²</span>
        <span onclick="set('Sum(n,(n,1,oo))')">∑ n</span>
        <span onclick="set('Sum(2**n,(n,0,oo))')">∑ 2^n</span>
        <span onclick="set('Sum(factorial(n),(n,0,oo))')">∑ n!</span>
        <span onclick="set('Sum(mobius(n),(n,1,oo))')">∑ μ(n)</span>
        <span onclick="set('Sum(fibonacci(n),(n,1,oo))')">∑ F_n</span>
        <span onclick="set('Sum(liouville(n),(n,1,oo))')">∑ λ(n)</span>
        <span onclick="set('Sum(eulerphi(n),(n,1,oo))')">∑ φ(n)</span>
        <span onclick="set('Sum((-1)**(n+1)*n**2,(n,1,oo))')">交替平方</span>
    </div>
    <input id="query" value="Sum(n**2,(n,1,oo))" placeholder="输入发散级数">
    <div class="btn-group">
        <button class="btn-calc" onclick="doCalc()">坍缩!</button>
        <button class="btn-viz" onclick="doViz()">可视化</button>
    </div>
    <div class="loading" id="loading">⏳ 塌缩怪兽正在九域中搜索对偶映射...</div>
    <pre id="result"></pre>
    <img id="graph" src="" alt="九域映射图">
</div>
<script>
    function set(text) { document.getElementById('query').value = text; }
    function fmt(v) {
        if (v === null || v === undefined) return '未知';
        if (v === Infinity || v === '∞') return '∞';
        if (typeof v === 'string' && v === '∞') return '∞';
        if (Math.abs(v) < 1e-10) return '0';
        if (v === -1/12) return '-1/12';
        if (v === 0.5) return '1/2';
        if (v === -1) return '-1';
        if (v === 1/3) return '1/3';
        if (v === 0.25) return '1/4';
        if (v === 1/120) return '1/120';
        if (v === -0.5) return '-1/2';
        if (v === -0.125) return '-1/8';
        if (v === -1/24) return '-1/24';
        if (v === -1/252) return '-1/252';
        if (v === 1/252) return '1/252';
        if (v === -5/6) return '-5/6';
        if (v === 1/144) return '1/144';
        if (v === 0.5963473623231941) return '≈0.596';
        var phi = (Math.sqrt(5)+1)/2;
        if (Math.abs(v-phi) < 1e-10) return 'φ';
        return v.toFixed(6);
    }
    async function doCalc() {
        var q = document.getElementById('query').value;
        var r = document.getElementById('result');
        var l = document.getElementById('loading');
        l.classList.add('show');
        try {
            var resp = await fetch('/api/calc', {
                method:'POST',
                headers:{'Content-Type':'application/json'},
                body:JSON.stringify({query:q})
            });
            var d = await resp.json();
            l.classList.remove('show');
            if (d.status === 'success') {
                r.innerText =
                    '📊 计算报告\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '输入:     ' + d.input + '\\n' +
                    '通项:     ' + d.summand + '\\n' +
                    '发散类型: ' + d.divergence_type + '\\n' +
                    '原始域:   ' + d.domain + '\\n' +
                    '映射路径: ' + d.mapping_path + '\\n' +
                    '映射步数: ' + d.steps + ' 步（九域直径 ≤ 3）\\n' +
                    '计算方法: ' + (d.method || 'N/A') + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '坍缩值:   ' + fmt(d.value) + '\\n' +
                    '━━━━━━━━━━━━━━━━━━━━\\n' +
                    '任何发散问题欢迎私信：永恒无限鱼(全平台)，合作15299667123(同微)';
            } else {
                r.innerText = '⚠ ' + (d.message || '未知错误');
            }
        } catch(e) {
            l.classList.remove('show');
            r.innerText = '⚠ 网络错误: ' + e.message;
        }
    }
    async function doViz() {
        var q = document.getElementById('query').value;
        var g = document.getElementById('graph');
        var l = document.getElementById('loading');
        l.classList.add('show');
        try {
            var resp = await fetch('/api/viz', {
                method:'POST',
                headers:{'Content-Type':'application/json'},
                body:JSON.stringify({query:q})
            });
            var d = await resp.json();
            l.classList.remove('show');
            if (d.graph_url) {
                g.src = d.graph_url + '?t=' + Date.now();
                g.style.display = 'block';
            }
        } catch(e) {
            l.classList.remove('show');
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
        query = data.get('query', '')
        result = compute(query)
        if result.get('status') == 'unknown':
            log_unresolved(query, result.get('message', ''))
        return jsonify(result)
    except Exception as e:
        traceback.print_exc()
        return jsonify({"status": "error", "message": f"服务器内部错误: {str(e)}"})

@app.route('/api/viz', methods=['POST'])
def api_viz():
    try:
        data = request.get_json(silent=True) or {}
        query = data.get('query', '')
        result = compute(query)
        if result.get('status') == 'success' and result.get('mapping_path'):
            path = result['mapping_path'].split(' → ')
            graph_url = visualize_mapping(path)
            if graph_url:
                return jsonify({'graph_url': graph_url})
        return jsonify({'graph_url': None})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'graph_url': None, 'error': str(e)})

if __name__ == "__main__":
    os.makedirs('static', exist_ok=True)
    port = int(os.environ.get("PORT", 5000))
    print(f"🧌 塌缩怪兽 v11.5 已启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
