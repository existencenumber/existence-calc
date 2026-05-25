"""
存在数论塌缩怪兽 v9.0 — 云部署版
基于Flask的Web服务，支持自主推理与九域可视化
"""

import re, math, os, json
from collections import deque

# ── 可选依赖 ──
try:
    import sympy as sp
    from sympy import Sum, oo, factorial, log, Symbol, simplify
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from flask import Flask, request, jsonify, render_template_string

# ═══════════════════════════════════════════
# 九域对偶图
# ═══════════════════════════════════════════
DUAL_GRAPH = {
    "加法域": {"指数映射": "乘法域", "黎曼和极限": "积分域", "差商极限": "微分域"},
    "乘法域": {"对数映射": "加法域", "对数导数": "积分域", "梅林变换": "谱域"},
    "积分域": {"拉普拉斯变换": "谱域", "泛函极限": "泛函积分域", "积分的逆": "微分域"},
    "微分域": {"傅里叶变换": "谱域", "微分的逆": "积分域", "差商极限的逆": "加法域"},
    "谱域": {"逆傅里叶变换": "微分域", "逆拉普拉斯变换": "积分域", "逆梅林变换": "乘法域"},
    "泛函积分域": {"泛函极限的逆": "积分域", "二维拓扑": "编织域"},
    "编织域": {"辫子同伦": "同伦域", "二维拓扑的逆": "泛函积分域"},
    "同伦域": {"态射范畴化": "范畴域", "辫子同伦的逆": "编织域"},
    "范畴域": {"恒等态射对应0": "加法域", "恒等态射对应1": "乘法域", "态射范畴化的逆": "同伦域"},
}
CONVERGENCE_DOMAINS = {"谱域", "泛函积分域", "同伦域", "范畴域"}

DOMAIN_POSITIONS = {
    "加法域": (0, 1), "乘法域": (2, 1), "积分域": (1, 2),
    "微分域": (-1, 2), "谱域": (1, 3), "泛函积分域": (-1, 0),
    "编织域": (1, 0), "同伦域": (0, -1), "范畴域": (2, -1)
}
DOMAIN_COLORS = {
    "加法域": "#FF6B6B", "乘法域": "#4ECDC4", "积分域": "#45B7D1",
    "微分域": "#96CEB4", "谱域": "#FFEAA7", "泛函积分域": "#DDA0DD",
    "编织域": "#98D8C8", "同伦域": "#F7DC6F", "范畴域": "#BB8FCE"
}

# ═══════════════════════════════════════════
# 谱域计算函数
# ═══════════════════════════════════════════
def compute_zeta(s):
    zeta_vals = {0:-0.5, -1:-1/12, -2:0, -3:1/120, -4:-1/252, -5:1/252, -6:0, -7:-1/240, -8:0,
                 1:float('inf'), 2:math.pi**2/6, 3:1.2020569031595942}
    if s in zeta_vals: return zeta_vals[s]
    if -8 <= s <= -1: return 0
    return None

def compute_eta(s):
    return {0:0.5, -1:0.25, -2:0, -3:-1/8, 1:math.log(2)}.get(s)

def compute_abel(r): return 1/(1-r)
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

def compute_zeta_deriv(s):
    if s == 0: return 0.5*math.log(2*math.pi)
    return None

def compute_prime_zeta(s):
    m = {-1:-0.0942, -2:-0.023, -3:-0.008, -4:-0.003}
    if s in m: return m[s]
    if s > 0: return None
    return 0

def compute_simplex(k):
    return compute_zeta(-k) if k >= 2 else compute_zeta(-1)

# ═══════════════════════════════════════════
# BFS域映射
# ═══════════════════════════════════════════
def find_shortest_path(source_domain):
    if source_domain in CONVERGENCE_DOMAINS:
        return [source_domain], 0
    queue = deque([(source_domain, [source_domain])])
    visited = {source_domain}
    while queue:
        current, path = queue.popleft()
        for _, neighbor in DUAL_GRAPH.get(current, {}).items():
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = path + [neighbor]
                if neighbor in CONVERGENCE_DOMAINS:
                    return new_path, len(new_path) - 1
                queue.append((neighbor, new_path))
    return [source_domain, "谱域"], 1

# ═══════════════════════════════════════════
# 自主推理引擎
# ═══════════════════════════════════════════
def analyze_term(expr, var='n'):
    expr_s = sp.simplify(expr)
    if expr_s.is_polynomial(var):
        deg = sp.degree(expr_s, gen=var)
        if deg >= 0:
            return ('polynomial', deg, '加法域',
                    '加法域 → 指数映射 → 乘法域 → 梅林变换 → 谱域', 2, 'zeta', -deg)
    atoms = expr_s.atoms()
    for a in atoms:
        if a.is_Pow and a.exp == Symbol(var):
            base = a.base
            if base.is_Integer and base > 1:
                return ('geometric', int(base), '乘法域',
                        '乘法域 → 泛函积分域（阿贝尔求和）', 1, 'abel', int(base))
    if expr_s.has(sp.factorial):
        return ('factorial', None, '乘法域',
                '乘法域 → 泛函积分域（波雷尔求和）', 1, 'borel', None)
    if expr_s == 1/Symbol(var):
        return ('harmonic', None, '加法域',
                '加法域 → 梅林变换 → 谱域', 1, 'zeta', 1)
    if expr_s == sp.log(Symbol(var)):
        return ('logarithmic', None, '加法域',
                '加法域 → 梅林变换 → 谱域', 1, 'zeta_deriv', 0)
    fn = str(expr_s.func) if hasattr(expr_s, 'func') else ''
    if 'mobius' in fn.lower(): return ('mobius', None, '加法域', '加法域 → 梅林变换 → 谱域', 1, 'mobius', None)
    if 'liouville' in fn.lower(): return ('liouville', None, '加法域', '加法域 → 梅林变换 → 谱域', 1, 'liouville', None)
    if 'eulerphi' in fn.lower(): return ('euler_phi', None, '加法域', '加法域 → 梅林变换 → 谱域', 1, 'euler_phi', None)
    if 'divisor' in fn.lower(): return ('divisor', None, '加法域', '加法域 → 梅林变换 → 谱域', 1, 'divisor', None)
    if 'fibonacci' in fn.lower(): return ('fibonacci', None, '乘法域', '乘法域 → 泛函积分域', 1, 'fibonacci', None)
    return ('unknown', None, '未知', '无映射', 0, None, None)

def compute_spectral(method, param):
    methods = {
        'zeta': compute_zeta, 'eta': compute_eta, 'abel': compute_abel,
        'euler': compute_euler, 'borel': compute_borel, 'zeta_deriv': compute_zeta_deriv,
        'fibonacci': compute_fibonacci, 'mobius': compute_mobius, 'liouville': compute_liouville,
        'euler_phi': compute_euler_phi, 'mangoldt': compute_mangoldt, 'divisor': compute_divisor,
        'superfactorial': compute_superfactorial, 'primorial': compute_primorial,
        'euler_char': compute_euler_char, 'ramanujan_cf': compute_ramanujan_cf,
        'gr_1loop': compute_gr_1loop, 'gr_2loop': compute_gr_2loop,
        'prime_counting': compute_prime_counting, 'prime_zeta': compute_prime_zeta,
        'simplex': compute_simplex,
    }
    if method in methods:
        return methods[method](param) if param is not None else methods[method]()
    return None

def compute(user_input):
    if not HAS_SYMPY:
        return {"status": "error", "message": "SymPy未安装，请联系管理员。"}
    s = user_input.strip()
    try:
        expr_str = s.replace('∑', 'Sum').replace('∞', 'oo').replace(' ', '')
        n = Symbol('n')
        expr = sp.sympify(expr_str)
        if not isinstance(expr, Sum):
            expr = Sum(expr, (n, 1, oo))
        if isinstance(expr, Sum):
            summand = expr.args[0]
            var_tuple = expr.args[1]
            var_name = str(var_tuple[0])
            upper = var_tuple[2]
            if upper != oo:
                return {"status": "finite", "message": f"该级数上限为{str(upper)}，可直接计算。", "summand": str(summand)}
        else:
            return {"status": "error", "message": "无法识别为无穷级数。"}
        asym_type, param, domain, mapping, steps, method, method_param = analyze_term(summand, var_name)
        if asym_type == 'unknown':
            return {"status": "unknown", "input": s, "summand": str(summand),
                    "message": "无法识别通项的类型。请尝试已知数论函数名。"}
        value = compute_spectral(method, method_param)
        return {"input": s, "summand": str(summand), "domain": domain,
                "divergence": asym_type, "mapping_path": mapping, "steps": steps,
                "method": method, "param": method_param, "value": value,
                "status": "success" if value is not None else "unknown"}
    except Exception as e:
        return {"status": "error", "message": f"解析错误: {str(e)}。请尝试 Sum(n**2, (n,1,oo)) 格式。"}

# ═══════════════════════════════════════════
# 可视化
# ═══════════════════════════════════════════
def visualize_mapping(source_domain, mapping_path, save_path="static/graph.png"):
    if not HAS_MPL:
        return None
    os.makedirs('static', exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(-3, 4); ax.set_ylim(-2, 5); ax.set_aspect('equal'); ax.axis('off')
    for domain, (x, y) in DOMAIN_POSITIONS.items():
        color = DOMAIN_COLORS.get(domain, "#CCCCCC")
        c = plt.Circle((x, y), 0.4, color=color, ec='black', linewidth=2, zorder=2)
        ax.add_patch(c)
        ax.text(x, y, domain, ha='center', va='center', fontsize=9, weight='bold')
    for src, edges in DUAL_GRAPH.items():
        if src not in DOMAIN_POSITIONS: continue
        x1, y1 = DOMAIN_POSITIONS[src]
        for _, dst in edges.items():
            if dst not in DOMAIN_POSITIONS: continue
            x2, y2 = DOMAIN_POSITIONS[dst]
            ax.plot([x1, x2], [y1, y2], 'gray', linewidth=0.5, alpha=0.3, zorder=1)
    if len(mapping_path) >= 2:
        for i in range(len(mapping_path)-1):
            d1, d2 = mapping_path[i], mapping_path[i+1]
            if d1 in DOMAIN_POSITIONS and d2 in DOMAIN_POSITIONS:
                x1, y1 = DOMAIN_POSITIONS[d1]; x2, y2 = DOMAIN_POSITIONS[d2]
                ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                            arrowprops=dict(arrowstyle='->', color='red', lw=2.5, zorder=3))
    if mapping_path:
        start = mapping_path[0]
        if start in DOMAIN_POSITIONS:
            x, y = DOMAIN_POSITIONS[start]
            ax.add_patch(plt.Circle((x, y), 0.45, color='none', ec='red', linewidth=3, zorder=4))
        end = mapping_path[-1]
        if end in DOMAIN_POSITIONS:
            x, y = DOMAIN_POSITIONS[end]
            ax.add_patch(plt.Circle((x, y), 0.45, color='none', ec='green', linewidth=3, zorder=4))
    plt.title(f"九域对偶映射\n{source_domain} → {' → '.join(mapping_path)}", fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    return save_path

# ═══════════════════════════════════════════
# 辅助函数
# ═══════════════════════════════════════════
def format_value(value):
    if value is None: return "未知"
    if isinstance(value, float) and value == float('inf'): return "∞"
    if isinstance(value, float) and abs(value) < 1e-10: return "0"
    if isinstance(value, float):
        frac_map = {-1/12:"-1/12", 0.5:"1/2", -1:"-1", 1/3:"1/3", 0.25:"1/4",
                    1/120:"1/120", -0.5:"-1/2", -0.125:"-1/8", -1/24:"-1/24",
                    -1/252:"-1/252", 1/252:"1/252", -5/6:"-5/6", 1/144:"1/144"}
        for num, frac_str in frac_map.items():
            if abs(value-num) < 1e-10: return f"{value} = {frac_str}"
        phi = (math.sqrt(5)+1)/2
        if abs(value-phi) < 1e-10: return f"{value} = φ（黄金比例）"
        return f"{value:.6f}"
    return str(value)

# ═══════════════════════════════════════════
# Flask 应用
# ═══════════════════════════════════════════
app = Flask(__name__)

HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>存在数论塌缩怪兽 v9.0</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, sans-serif; background: #0a0a1a; color: #e0e0e0;
               display: flex; justify-content: center; padding: 30px 15px; }
        .container { max-width: 800px; width: 100%; }
        h1 { text-align: center; color: #ff6b6b; margin-bottom: 6px; font-size: 2em; }
        .subtitle { text-align: center; color: #aaa; margin-bottom: 25px; font-size: 0.9em; }
        .input-group { display: flex; gap: 10px; margin-bottom: 18px; flex-wrap: wrap; }
        input[type="text"] { flex: 1; min-width: 250px; padding: 14px; font-size: 1.05em;
               border: 2px solid #333; background: #1a1a2e; color: #fff; border-radius: 8px; outline: none; }
        input[type="text"]:focus { border-color: #ff6b6b; }
        button { padding: 14px 24px; font-size: 1em; border: none; border-radius: 8px;
                 cursor: pointer; font-weight: bold; transition: all 0.3s; }
        .btn-calc { background: #ff6b6b; color: #fff; }
        .btn-calc:hover { background: #ff4444; }
        .btn-viz { background: #4ecdc4; color: #000; }
        .btn-viz:hover { background: #3dbdb5; }
        .examples { margin-bottom: 18px; line-height: 2; }
        .examples span { display: inline-block; background: #1a1a2e; padding: 6px 12px;
                         margin: 3px; border-radius: 18px; cursor: pointer; font-size: 0.85em;
                         border: 1px solid #333; transition: all 0.2s; }
        .examples span:hover { background: #ff6b6b; color: #fff; border-color: #ff6b6b; }
        .result { background: #1a1a2e; border-radius: 12px; padding: 25px; margin-top: 20px;
                  border: 1px solid #333; display: none; }
        .result.show { display: block; }
        .result h3 { color: #ff6b6b; margin-bottom: 15px; }
        .result-row { display: flex; justify-content: space-between; padding: 8px 0;
                      border-bottom: 1px solid #222; flex-wrap: wrap; }
        .result-label { color: #aaa; }
        .result-value { color: #fff; font-weight: bold; }
        .result-highlight { font-size: 1.5em; color: #4ecdc4; text-align: center;
                            padding: 15px 0; border-bottom: 1px solid #222; }
        .graph-img { max-width: 100%; margin-top: 20px; border-radius: 8px; display: none; }
        .graph-img.show { display: block; }
        .loading { text-align: center; padding: 20px; color: #aaa; display: none; }
        .loading.show { display: block; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧌 塌缩怪兽 v9.0</h1>
        <p class="subtitle">存在数论非微扰计算器 | 九域对偶映射 | e<sup>iS</sup> = 1</p>
        <div class="examples">
            <span onclick="set('Sum(n**2, (n, 1, oo))')">∑ n²</span>
            <span onclick="set('Sum(n, (n, 1, oo))')">∑ n</span>
            <span onclick="set('Sum(2**n, (n, 0, oo))')">∑ 2^n</span>
            <span onclick="set('Sum(factorial(n), (n, 0, oo))')">∑ n!</span>
            <span onclick="set('Sum(mobius(n), (n, 1, oo))')">∑ μ(n)</span>
            <span onclick="set('Sum(fibonacci(n), (n, 1, oo))')">∑ F_n</span>
            <span onclick="set('Sum(liouville(n), (n, 1, oo))')">∑ λ(n)</span>
            <span onclick="set('Sum(eulerphi(n), (n, 1, oo))')">∑ φ(n)</span>
        </div>
        <div class="input-group">
            <input type="text" id="query" placeholder="输入发散级数，如 Sum(n**2, (n,1,oo)) 或 1+2+3+..."
                   value="Sum(n**2, (n, 1, oo))">
            <button class="btn-calc" onclick="doCalc()">坍缩!</button>
            <button class="btn-viz" onclick="doViz()">可视化</button>
        </div>
        <div class="loading" id="loading">⏳ 塌缩怪兽正在九域中搜索对偶映射...</div>
        <div class="result" id="result"></div>
        <img class="graph-img" id="graph" src="" alt="九域映射图">
    </div>
    <script>
        function set(text) { document.getElementById('query').value = text; }
        function fmt(v) {
            if (v === null || v === undefined) return '未知';
            if (v === Infinity || v === '∞') return '∞';
            if (Math.abs(v) < 1e-10) return '0';
            var fracMap = {};
            fracMap[-1/12] = '-1/12'; fracMap[0.5] = '1/2'; fracMap[-1] = '-1';
            fracMap[1/3] = '1/3'; fracMap[0.25] = '1/4'; fracMap[1/120] = '1/120';
            fracMap[-0.5] = '-1/2'; fracMap[-0.125] = '-1/8'; fracMap[-1/24] = '-1/24';
            fracMap[-1/252] = '-1/252'; fracMap[1/252] = '1/252'; fracMap[-5/6] = '-5/6';
            fracMap[1/144] = '1/144';
            for (var k in fracMap) {
                if (Math.abs(v - parseFloat(k)) < 1e-10) return v + ' = ' + fracMap[k];
            }
            var phi = (Math.sqrt(5)+1)/2;
            if (Math.abs(v - phi) < 1e-10) return v + ' = φ';
            return typeof v === 'number' ? v.toFixed(6) : String(v);
        }
        function buildHTML(d) {
            var h = '<h3>📊 计算报告</h3>';
            h += '<div class="result-row"><span class="result-label">输入</span><span class="result-value">' + (d.input || 'N/A') + '</span></div>';
            h += '<div class="result-row"><span class="result-label">通项</span><span class="result-value">' + (d.summand || 'N/A') + '</span></div>';
            h += '<div class="result-row"><span class="result-label">原始域</span><span class="result-value">' + (d.domain || 'N/A') + '（' + (d.divergence || 'N/A') + '）</span></div>';
            h += '<div class="result-row"><span class="result-label">映射路径</span><span class="result-value">' + (d.mapping_path || 'N/A') + '</span></div>';
            h += '<div class="result-row"><span class="result-label">映射步数</span><span class="result-value">' + (d.steps || 'N/A') + ' 步（九域直径 ≤ 3）</span></div>';
            h += '<div class="result-highlight">坍缩值: ' + fmt(d.value) + '</div>';
            h += '<p style="text-align:center;color:#aaa;margin-top:10px;">发散是表象，守恒是本质。e^{iS} = 1</p>';
            return h;
        }
        async function doCalc() {
            var q = document.getElementById('query').value;
            var r = document.getElementById('result');
            var l = document.getElementById('loading');
            var g = document.getElementById('graph');
            g.classList.remove('show');
            r.classList.remove('show');
            l.classList.add('show');
            try {
                var resp = await fetch('/api/calc', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({query: q})
                });
                var d = await resp.json();
                l.classList.remove('show');
                r.innerHTML = buildHTML(d);
                r.classList.add('show');
            } catch(e) {
                l.classList.remove('show');
                r.innerHTML = '<h3>错误</h3><p>网络错误: ' + e.message + '</p>';
                r.classList.add('show');
            }
        }
        async function doViz() {
            var q = document.getElementById('query').value;
            var g = document.getElementById('graph');
            var l = document.getElementById('loading');
            g.classList.remove('show');
            l.classList.add('show');
            try {
                var resp = await fetch('/api/viz', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({query: q})
                });
                var d = await resp.json();
                l.classList.remove('show');
                if (d.graph_url) {
                    g.src = d.graph_url;
                    g.classList.add('show');
                }
            } catch(e) {
                l.classList.remove('show');
                alert('可视化失败: ' + e.message);
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
    data = request.get_json()
    query = data.get('query', '')
    result = compute(query)
    return jsonify(result)

@app.route('/api/viz', methods=['POST'])
def api_viz():
    data = request.get_json()
    query = data.get('query', '')
    result = compute(query)
    if result.get('status') in ('success', 'unknown') and result.get('domain'):
        path_parts = result['mapping_path'].split(' → ')
        graph_path = visualize_mapping(result['domain'], path_parts, 'static/graph.png')
        if graph_path:
            return jsonify({'graph_url': '/static/graph.png'})
    return jsonify({'graph_url': None})

# ═══════════════════════════════════════════
# 启动
# ═══════════════════════════════════════════
if __name__ == "__main__":
    os.makedirs('static', exist_ok=True)
    port = int(os.environ.get("PORT", 5000))
    print(f"🧌 塌缩怪兽已启动: http://0.0.0.0:{port}")
    app.run(host="0.0.0.0", port=port, debug=False)
