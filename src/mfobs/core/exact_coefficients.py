"""Exact rational coefficients, optionally over declared conserved parameters."""
import ast
from fractions import Fraction
import re
import sympy as sp


def coefficient(value, parameters=()):
    if not parameters:
        f = Fraction(value)
        return sp.Rational(f.numerator, f.denominator)
    if any(not isinstance(name, str) or not re.fullmatch('[A-Za-z_][A-Za-z_0-9]*', name) for name in parameters):
        raise ValueError('Invalid exact coefficient parameter')
    symbols = {name: sp.Symbol(name) for name in parameters}
    tree = ast.parse(str(value), mode='eval')
    if sum(1 for _ in ast.walk(tree)) > 1000:
        raise ValueError('Coefficient expression too large')
    def read(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name) and node.id in symbols:
            return symbols[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return read(node.operand)*(1 if isinstance(node.op, ast.UAdd) else -1)
        if isinstance(node, ast.BinOp):
            a, b = read(node.left), read(node.right)
            if isinstance(node.op, ast.Add): return a+b
            if isinstance(node.op, ast.Sub): return a-b
            if isinstance(node.op, ast.Mult): return a*b
            if isinstance(node.op, ast.Div) and b != 0: return a/b
            if isinstance(node.op, ast.Pow) and b.is_Integer and -16 <= b <= 16 and not (a == 0 and b < 0): return a**b
        raise ValueError('Only exact rational arithmetic over declared parameters is permitted')
    return sp.cancel(read(tree.body))


def rational_vector(values):
    """A sparse vector with exact coefficients and zero entries removed."""
    return {key: Fraction(value) for key, value in values.items() if value}


def add_vectors(left, right, scale=1):
    """Add one scaled sparse rational vector; keys describe formal coordinates."""
    result = rational_vector(left)
    scale = Fraction(scale)
    for key, value in right.items():
        result[key] = result.get(key, Fraction()) + scale * Fraction(value)
        if not result[key]:
            del result[key]
    return result
