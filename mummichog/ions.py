'''
Temporary solution form mummichog version 2. Version 3 will use mass2chem lib. 
'''

import re

PROTON = 1.00727646677
electron = 0.000549

def parse_chemformula(x):
    '''This does not deal with nested groups in chemical formula.
    Formula from HMDB 3.5 are compatible.
    '''
    p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
    d = {}
    for pair in p: d[pair[0]] = int( pair[1] or 1 )
    return d

def check_sub(Fm1, Fm2):
    '''Check if Fm2 a subunit of Fm1; Dictionary representation of chem formula
    '''
    Fm2_in_Fm1 = True
    Elements1 = Fm1.keys()
    Elements2 = Fm2.keys()
    if [x for x in Elements2 if x not in Elements1]:    #element not in Fm1
        return False
    else:
        for e in Elements2:
            if Fm2[e] > Fm1[e]: Fm2_in_Fm1 = False
        return Fm2_in_Fm1

def compute_adducts(mw, cFormula, mode='pos'):
    '''
    function calculating isotopes and adducts.
    # A few rules here:
    # C13, S34, Cl37 pos only applicable if the atom is in formula
    # and others in `required subgroup`, as the last item in the tuples
    # Better way is to based on chemical structure - future direction.
    mode could be passed as 'pos_default'.
    Temporary solution form mummichog version 2. Version 3 will use mass2chem lib.
    '''
    if 'pos' in mode:
        addList = [(mw - electron, 'M[1+]', ''), 
            (mw + PROTON, 'M+H[1+]', ''),
            (mw + 21.9820 + PROTON, 'M+Na[1+]', ''), 
            (mw + 18.0106 + PROTON, 'M+H2O+H[1+]', ''), 
            (mw +1.0034 - electron, 'M(C13)[1+]', ''),
            (mw +1.0034 + PROTON, 'M(C13)+H[1+]', 'C'),
            (mw/2 + PROTON, 'M+2H[2+]', ''),
            (mw/3 + PROTON, 'M+3H[3+]', ''),

            (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', 'C'),
            (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', 'C'),
            (mw +1.9958 + PROTON, 'M(S34)+H[1+]', 'S'),
            (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', 'Cl'),
            
            (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]', ''),
            (mw + 37.9555 + PROTON, 'M+K[1+]', ''), 
            
            (mw - 18.0106 + PROTON, 'M-H2O+H[1+]', 'H2O'), 
            (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]', 'H4O2'),
            (mw - 17.0265 + PROTON, 'M-NH3+H[1+]', 'NH3'),
            (mw - 27.9950 + PROTON, 'M-CO+H[1+]', 'CO'),
            (mw - 43.9898 + PROTON, 'M-CO2+H[1+]', 'CO2'),
            (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]', 'H2CO2'),
            (mw + 67.9874 + PROTON, 'M+HCOONa[1+]', ''),
            (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]', 'HCO2Na'),
            (mw + 57.9586 + PROTON, 'M+NaCl[1+]', ''), 
            (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]', 'C3H4O2'),
            (mw + 83.9613 + PROTON, 'M+HCOOK[1+]', ''),
            (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]', 'HCO2K'),
            ]
    else:
        addList =  [(mw - PROTON, 'M-H[-]', ''),
            (mw + electron, 'M[-]', ''), 
            (mw - 18.0106 - PROTON, 'M-H2O-H[-]', 'H2O'),
            (mw + 34.9689, 'M+Cl[-]', ''),

            (mw + 1.0034 - PROTON, 'M(C13)-H[-]', 'C'),
            (mw + 36.9659, 'M+Cl37[-]', ''),
            (mw/2 - PROTON, 'M-2H[2-]', ''),
            
            (mw + 1.9958 - PROTON, 'M(S34)-H[-]', 'S'),
            (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]', 'Cl'),
            (mw + 21.9820 - 2*PROTON, 'M+Na-2H[-]', ''),
            (mw + 37.9555 - 2*PROTON, 'M+K-2H[-]', ''),
            
            (mw + 78.9183, 'M+Br[-]', ''),
            (mw + 80.9163, 'M+Br81[-]', ''),
            (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]', ''),
            (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]', ''),
            (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]', ''),
            (mw - PROTON + 15.99491, 'M-H+O[-]', ''),
            ]

    dict_cFormula = parse_chemformula(cFormula)
    mydict = {}
    for x in addList:
        if check_sub(dict_cFormula, parse_chemformula(x[2])):
            mydict[x[1]] = x[0]

    return mydict


def compute_isotopes_adducts(cpd, mode='pos'):
    '''
    cpd format: {'formula': 'C29H51N2O8PRS', 'mw': 0, 'name': 'linoelaidic acid ACP (all trans)', 'adducts': {}}
    return dictionary of isotopes and adducts
    '''
    mw = cpd.get('mw', 0)
    f = cpd.get('formula', '')
    if mw and f:
        return compute_adducts(mw, f, mode)
    else:
        return {}
