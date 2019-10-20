import tokenize
from io import StringIO
from pprint import pprint
from collections import namedtuple

import numpy as np

from patsy import demo_data
from patsy import dmatrix

from patsy.parse_formula import Operator
from patsy.parse_formula import _default_ops
from patsy.parse_formula import parse_formula
from patsy.parse_formula import _tokenize_formula
from patsy.parse_formula import _atomic_token_types

from patsy.infix_parser import infix_parse
from patsy.infix_parser import ParseNode

from patsy.desc import Evaluator
from patsy.desc import ModelDesc
from patsy.desc import IntermediateExpr
from patsy.desc import INTERCEPT

# this is our data structure for random effects terms
RandomEffectsTerm = namedtuple('RandomEffectsTerm', ['expr', 'factor'])


def eval_bar(evaluator, tree):
    """Evaluation function for the bar operator AST node."""

    assert len(tree.args) == 2
    expr_node, factor_node = tree.args

    # create model description for the expression left of the bar
    expr_node = ParseNode("~", None, [expr_node], expr_node.origin)
    expr_md = ModelDesc.from_formula(expr_node)

    # create model description for grouping factor right of the bar
    factor_node = ParseNode("~", None, [factor_node], factor_node.origin)
    factor_md = ModelDesc.from_formula(factor_node)
    factor_md.rhs_termlist.remove(INTERCEPT)

    # combine these in a random effects term
    ret = RandomEffectsTerm(expr=expr_md, factor=factor_md)

    # return corresponding intermediate expression
    return IntermediateExpr(False, None, False, [ret])


def buildzi(ret, data, env=0):
    """Build Z_i (corresponds to a single random effects term)."""
    X = dmatrix(ret.expr, data, env)
    J = dmatrix(ret.factor, data, env)
    n, _ = X.shape
    Zi = np.array([np.kron(J[i], X[i]) for i in range(n)])
    return Zi


def buildz(rets, data, env=0):
    Zis = [buildzi(ret, data, env=env) for ret in rets]
    return np.concatenate(Zis, axis=1).T


def evaluate_formula(formula):
    """Given mixed effects formula, return a model description."""
    # mixed effects specific operators
    extra_operators = [
        Operator('|', 2, 50),
        Operator('|| ', 2, 50)
    ]

    # construct a list of operator strings needed for tokenization
    operators = _default_ops + extra_operators
    operator_strings = [op.token_type for op in operators]

    tokens = list(_tokenize_formula(formula, operator_strings))
    node = infix_parse(tokens, operators, _atomic_token_types)

    e = Evaluator()

    # we can't handle double bar yet
    e.add_op('|', 2, eval_bar)

    model_description = e.eval(node, require_evalexpr=False)

    return model_description


def get_matrices(data, formula):
    """Given the data and a formula, build Z and X matrices."""
    model_description = evaluate_formula(formula)

    fixef_terms, randef_terms = [], []
    for term in model_description.rhs_termlist:
        if isinstance(term, RandomEffectsTerm):
            randef_terms.append(term)
        else:
            fixef_terms.append(term)

    Z = buildz(randef_terms, data)
    X = dmatrix(ModelDesc(model_description.lhs_termlist, fixef_terms), data)

    _, q = Z.shape

    L = np.eye(q)

    return X, Z, L
