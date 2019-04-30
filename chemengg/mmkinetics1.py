"""@package docstring
Documentation for this module 
This is an attempt at using python sympy to {\em derive} the Michaelis
Menten equation.
"""

from sympy import *

def mmstandard():
    """Documentation for a function 
    This is the standard michaelis menten mechanism.  The binding of enzyme
    to substrate to form the enzyme substrate complex is considered reversible 
    while the formation of the product from the enzyme substrate complex with the release
    of the enzyme is considered irreversible.
    \ce{E + S <=>T[\cf{k_{1}}][\cf{k_{-1}}] ES ->T[\cf{k_2}] E + P}
    """
    S = symbols('S')
    E = symbols('E')
    ES = symbols('ES')
    P = symbols('P')
    k1 = symbols('k_1')
    km1 = symbols('k_m1')
    k2 = symbols('k_2')
    E0 = symbols('E0')
    V = symbols('V')
    ES = symbols('ES')
    E = E0 - ES
    ES = solveset(k1*E*S-km1*ES-k2*ES,ES)
    LES = latex(ES)
    V = k2*ES
    LV = latex(V)
    return LV


def mmEPreversible():
    """Documentation for a function 
    This is the michaelis menten mechanism, with a twist.  The binding of enzyme
    to substrate to form the enzyme substrate complex is still considered reversible 
    but the formation of the product from the enzyme substrate complex with the release
    of the enzyme is also considered reversible.
    \ce{E + S <=>T[\cf{k_{1}}][\cf{k_{-1}}] ES <=>T[\cf{k_2}][\cf{k_{-2}}] E + P}
    """
    S = symbols('S')
    E = symbols('E')
    ES = symbols('ES')
    P = symbols('P')
    k1 = symbols('k_1')
    km1 = symbols('k_m1')
    k2 = symbols('k_2')
    km2 = symbols('k_m2')
    E0 = symbols('E0')
    V = symbols('V')
    ES = symbols('ES')
    E = E0 - ES
    ES = solveset(k1*E*S-km1*ES-k2*ES + km2*E*P,ES)
    LES = latex(ES)
    V = k2*ES
    LV = latex(V)
    return LV


