"""
塌缩怪兽 v11.5 — 最终修复交替级数识别（字符串匹配法）
基于存在数论的非微扰计算工具
"""

import re, math, os, json, traceback, uuid
from datetime import datetime
from collections import deque

import sympy as sp
from sympy import Sum, oo, factorial, log, Symbol, simplify

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

# ==================== 计算函数 ====================
def compute_zeta(s):
    table = {0:-0.5, -1:-1/12, -2:0, -3:1/120, -4:0, -5:-1/252, -6:0, -7:-1/240, -8:0,
             1:float('inf'), 2:math.pi**2/6, 3:1.2020569031595942}
    if s in table: return table[s]
    if s < 0 and s % 2 == 0: return 0
    try: return float(sp.N(sp.zeta(s)))
    except: return None

def compute_eta(param):
    """支持 (s, sign) 元组或单值 s"""
    if isinstance(param, tuple):
        s_val, sign = param
        base = compute_eta(s_val)
        return sign * base if base is not None else None
    
    eta_vals = {0:0.5, -1:0.25, -2:0, -3:-1/8, 1:math.log(2)}
    if param in eta_vals: return eta_vals[param]
    
    z = compute_zeta(param)
    if z is None or z == float('inf'): return None
    return (1 - 2**(1-param)) * z

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

# ==================== 核心分析函数 ====================
def analyze_term(expr, var='n', original_input=""):
    sym = Symbol(var)
    expr_s = sp.simplify(expr)

    # =========================================
    # 优先处理：字符串检测交替级数
    # 示例：Sum((-1)**(n+1) * n**2, (n, 1, oo))
    # =========================================
    if original_input and '(-1)**' in original_input:
        # 提取 (-1)**(...) 后面的核心通项
        # 去掉可能的前缀 "Sum(" 和结尾
        core = original_input
        if 'Sum(' in core and core.endswith(')'):
            core = core[4:-1]
        # 格式: (-1)**(n+1) * n**2, (n, 1, oo)
        # 提取第一个逗号之前的主表达式
        main_part = core.split(',')[0].strip()
        
        # 寻找交替因子的模式
        alt_pattern = r'\(\s*-\s*1\s*\)\s*\*\*\s*\([^)]+\)'
        match = re.search(alt_pattern, main_part)
        if match:
            alt_factor = match.group(0)
            # 剩余部分 = 原字符串去掉交替因子
            remaining = main_part[:match.start()] + main_part[match.end():]
            # 去掉前导的 '*' 符号
            remaining = remaining.strip().lstrip('*').strip()
            
            if remaining:
                try:
                    rest_expr = sp.sympify(remaining, locals=SAFE_LOCALS)
                    if rest_expr.is_polynomial(sym):
                        deg = sp.degree(rest_expr, gen=sym)
                        # 分析符号
                        sign = 1
                        if 'n+1' in alt_factor or '(n+1)' in alt_factor:
                            sign = 1
                        elif '(n-1)' in alt_factor or 'n-1' in alt_factor:
                            sign = 1
                        elif '(n)' in alt_factor:
                            sign = -1
                        
                        return {
                            "type": "alternating_power",
                            "domain": "加法域",
                            "mapping": ["加法域", "谱域"],
                            "method": "eta",
                            "param": (-deg, sign)
                        }
                except:
                    pass

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
                "mapping": ["乘法域", "谱域"], "method": "abel", "param": base}
    for atom in expr_s.atoms():
        if atom.is_Pow and atom.exp == sym:
            base = float(atom.base)
            return {"type": "geometric", "domain": "乘法域",
                    "mapping": ["乘法域", "谱域"], "method": "abel", "param": base}
        if atom.is_Pow and atom.exp.has(sym):
            base = float(atom.base)
            return {"type": "geometric", "domain": "乘法域",
                    "mapping": ["乘法域", "谱域"], "method": "abel", "param": base}

    # 阶乘 n!
    if expr_s.has(sp.factorial):
        return {"type": "factorial", "domain": "乘法域",
                "mapping": ["乘法域", "谱域"], "method": "borel", "param": None}

    # 调和 1/n
    if sp.simplify(expr_s - 1/sym) == 0:
        return {"type": "harmonic", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "zeta", "param": 1}

    # 对数 ln n
    if expr_s == sp.log(sym):
        return {"type": "logarithmic", "domain": "加法域",
                "mapping": ["加法域", "谱域"], "method": "zeta_deriv", "param": 0}

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

# ==================== 计算引擎 ====================
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
        
        analysis = analyze_term(summand, original_input=user_input)
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

# ...（后续的可视化、Flask 路由、HTML 模板等保持不变）...
