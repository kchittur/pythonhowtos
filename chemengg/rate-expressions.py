"""
Documentation for this module 
This is an attempt at using python sympy to {\em derive} the Michaelis
Menten equation.
Needs python sympy with solveset
'''

from sympy import *

def mmstandard():
    r'''This is the standard Michaelis Menten form  
    The binding of enzyme
    to substrate to form the enzyme substrate complex is considered reversible 
    while the formation of the product from the enzyme substrate 
    complex with the release
    of the enzyme is considered irreversible.
    .. math:: 
    \ce{E + S <=>T[\cf{k_{1}}][\cf{k_{-1}}] ES ->T[\cf{k_2}] E + P}
    
    ''' 
    S = symbols('S')
    E = symbols('E')
    ES = symbols('ES')
    P = symbols('P')
    k1 = symbols('k_1')
    km1 = symbols('k_{-1}')
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
    r'''Michaelis Menten with product formation step reversible 
    This is the michaelis menten mechanism, with a twist.  
    The binding of enzyme
    to substrate to form the enzyme substrate complex 
    is still considered reversible 
    but the formation of the product from the enzyme 
    substrate complex with the release
    of the enzyme is also considered reversible.
    .. math:: 
    \ce{E + S <=>T[\cf{k_{1}}][\cf{k_{-1}}] ES <=>T[\cf{k_2}][\cf{k_{-2}}] E + P}
    
    '''
    S = symbols('S')
    E = symbols('E')
    ES = symbols('ES')
    P = symbols('P')
    k1 = symbols('k_1')
    km1 = symbols('k_{-1}')
    k2 = symbols('k_2')
    km2 = symbols('k_{-2}')
    E0 = symbols('E0')
    V = symbols('V')
    ES = symbols('ES')
    E = E0 - ES
    ES = solveset(k1*E*S-km1*ES-k2*ES + km2*E*P,ES)
    LES = latex(ES)
    V = k2*ES - km2*E*P
    LV = latex(V)
    return LV

def mmEqlb():
    r''' derivation using the equilibrium assumption
    .. math::  
    \ce{E + S <=>T[\cf{k_{1}}][\cf{k_{-1}}] ES ->T[\cf{k_2}] E + P}
    where \frac{k_{1}{k_{-1}} = K_{m} = \frac{ES]{[E][S]}

    '''
    S = symbols('S')
    E = symbols('E')
    ES = symbols('ES')
    P = symbols('P')
    k1 = symbols('k_1')
    km1 = symbols('k_{-1}')
    k2 = symbols('k_2')
    E0 = symbols('E0')
    V = symbols('V')
    ES = symbols('ES')
    Km = symbols('K_m')
    E = E0 - ES
    Km = km1/k1
    ES = solveset(Km*ES - E*S, ES)
    LES = latex(ES)
    V = k2*ES
    LV = latex(V)
    return LV    


standard = mmstandard()
mmEPreversible = mmEPreversible()
mmE = mmEqlb()
