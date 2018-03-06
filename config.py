# Copyright (c) 2010-2018 Shuzhao Li.
# All rights reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


'''
configuration and utility functions of mummichog
@author: Shuzhao Li

'''

VERSION = '2.0.5-beta-20180306' 
RELEASE = False
USE_DEBUG = False


import sys, os, inspect
import numpy as np


SEARCH_STEPS = 4
MODULE_SIZE_LIMIT = 100
SIGNIFICANCE_CUTOFF = 0.05
MASS_RANGE = (50, 2000)


RETENTION_TIME_TOLERANCE_FRAC = 0.01    # fraction of total retention time


PROTON = 1.00727646677


# currency metaboites of ubiquitous presence
currency = ['C00001', 'C00080', 'C00007', 'C00006', 'C00005', 'C00003',
            'C00004', 'C00002', 'C00013', 'C00008', 'C00009', 'C00011', 
            'G11113', '',
            'H2O', 'H+', 'Oxygen', 'NADP+', 'NADPH', 'NAD+', 'NADH', 'ATP', 
            'Pyrophosphate', 'ADP', 'Orthophosphate', 'CO2',]

# testing purpose only
# currency = []


shapeslist = ['diamond', 'octagon', 'invtrapezium', 'triangle', 'box', 
        'trapezium', 'invtriangle', 'parallelogram', 'polygon','egg']

#zcolors = ['#0000FF', '#0066FF', '#00FFFF', '#66FFFF', '#99FFFF','#FFFFFF',  
#           '#FFFF99', '#FFCC00', '#FF6633', '#FF3300', '#CC0000', '#CC0000',]

zcolors = ['#99BBFF', '#99CCFF', '#99DDFF', '#99EEFF', '#99FFFF',
           '#FFFFFF',
           '#FFEE44', '#FFCC44', '#FFBB44', '#FFAA44', '#FF9944', 
           #'#FF8844',
           ]


def sigmoid(x, theta=0.25):
    '''
    A logistic function for weighting. Not used.
    '''
    return 1/(1 + np.e ** (-theta * x))


# accuracy of the MS instrument
def mz_tolerance(mz, mode):
    '''
    The functions supplied here for FTMS and Oribtrap 
    were based on instrumentation in DPJ lab.
    Please verify or replace them before using on your own data.
    It's better to use a fixed value, e.g. 10 ppm, in case of uncertainty.
    
    will use flat ppm, not using this function for ------
    
    '''
    try:
        mode = int(mode)
        return 0.000001 * mode * mz
    
    except ValueError:
        if mode == 'FTMS':
            return max(0.00001*mz, 3*2**(1.98816*np.log2(mz) - 26.1945))
        
        elif mode == 'ORBITRAP':
            # needs further calibration
            # return max(0.00001*mz, 2*2**(1.45325*np.log2(mz) - 20.8554))
            return 0.000010*mz
        
        else:
            return 0.000010*mz
  
primary_ions = ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]']

wanted_adduct_list = {
    'pos_default': ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],
    'dpj_positive': ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]',
                    'M(Cl37)+H[1+]', 'M(S34)+H[1+]', 'M+K[1+]', 'M+HCOOK[1+]', 
                    ],
    'generic_positive': ['M[1+]','M+H[1+]','M+2H[2+]','M+3H[3+]','M(C13)+H[1+]','M(C13)+2H[2+]',
                    'M(C13)+3H[3+]','M(S34)+H[1+]','M(Cl37)+H[1+]','M+Na[1+]','M+H+Na[2+]','M+K[1+]',
                    'M+H2O+H[1+]','M-H2O+H[1+]','M-H4O2+H[1+]','M-NH3+H[1+]','M-CO+H[1+]',
                    'M-CO2+H[1+]','M-HCOOH+H[1+]','M+HCOONa[1+]','M-HCOONa+H[1+]','M+NaCl[1+]',
                    'M-C3H4O2+H[1+]','M+HCOOK[1+]','M-HCOOK+H[1+]',
                    ],
    'negative': ['M-H[-]','M-2H[2-]','M(C13)-H[-]','M(S34)-H[-]','M(Cl37)-H[-]',
                    'M+Na-2H[-]','M+K-2H[-]','M-H2O-H[-]','M+Cl[-]','M+Cl37[-]',
                    'M+Br[-]','M+Br81[-]','M+ACN-H[-]','M+HCOO[-]','M+CH3COO[-]','M-H+O[-]'
    
                    ],

    # to add options

    }


def adduct_function(mw, mode):
    '''
    return a list of derivatives/adducts according to operation mode.
    The most frequent derivatives under positive mode are adopted from
        Brown et al. Analyst, 2009, 134, 1322-1332.
    'dpj_positive' is a customized version for DPJ lab.
    Negative mode is an empirical compilation.
    
    Some derivatives are not possible for some compounds, 
    subject to future upgrade.
    
    
    Paul Benton sent a list used in XCMSonline.
    
    This is to be replaced by pre-computed tables based on chemical formula.
    
    '''
    if mode == 'dpj_positive':
        return [(mw, 'M[1+]'), 
                (mw + PROTON, 'M+H[1+]'),
                (mw/2 + PROTON, 'M+2H[2+]'),
                (mw +1.0034 + PROTON, 'M(C13)+H[1+]'),
                (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]'),
                (mw +1.9958 + PROTON, 'M(S34)+H[1+]'),
                (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]'),
                (mw + 21.9820 + PROTON, 'M+Na[1+]'), 
                (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]'),
                (mw + 37.9555 + PROTON, 'M+K[1+]'), 
                (mw + 67.9874 + PROTON, 'M+HCOONa[1+]'),
                (mw + 83.9613 + PROTON, 'M+HCOOK[1+]'),
                ]
    
    elif mode == 'generic_positive':
        return [(mw, 'M[1+]'), 
                (mw + PROTON, 'M+H[1+]'),
                (mw/2 + PROTON, 'M+2H[2+]'),
                (mw/3 + PROTON, 'M+3H[3+]'),
                (mw +1.0034 + PROTON, 'M(C13)+H[1+]'),
                (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]'),
                (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]'),
                (mw +1.9958 + PROTON, 'M(S34)+H[1+]'),
                (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]'),
                (mw + 21.9820 + PROTON, 'M+Na[1+]'), 
                (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]'),
                (mw + 37.9555 + PROTON, 'M+K[1+]'), 
                (mw + 18.0106 + PROTON, 'M+H2O+H[1+]'), 
                (mw - 18.0106 + PROTON, 'M-H2O+H[1+]'), 
                (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]'),
                (mw - 17.0265 + PROTON, 'M-NH3+H[1+]'),
                (mw - 27.9950 + PROTON, 'M-CO+H[1+]'),
                (mw - 43.9898 + PROTON, 'M-CO2+H[1+]'),
                (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]'),
                (mw + 67.9874 + PROTON, 'M+HCOONa[1+]'),
                (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]'),
                (mw + 57.9586 + PROTON, 'M+NaCl[1+]'), 
                (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]'),
                (mw + 83.9613 + PROTON, 'M+HCOOK[1+]'),
                (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]'),
                ]
    
    elif mode == 'park':
        return [(mw + PROTON, 'M+H[1+]'),
                (mw + 21.9820 + PROTON, 'M+Na[1+]'), 
                #(mw + 18.0106 + PROTON, 'M+H2O+H[1+]'), 
                (mw - 18.0106 + PROTON, 'M-H2O+H[1+]'), 
                (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]'),
                ]
        
    elif mode == 'negative':
        return [(mw - PROTON, 'M-H[-]'),
               (mw/2 - PROTON, 'M-2H[2-]'),
               (mw + 1.0034 - PROTON, 'M(C13)-H[-]'),
               (mw + 1.9958 - PROTON, 'M(S34)-H[-]'),
               (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]'),
               (mw + 21.9820 - 2*PROTON, 'M+Na-2H[-]'),
               (mw + 37.9555 - 2*PROTON, 'M+K-2H[-]'),
               (mw - 18.0106 - PROTON, 'M-H2O-H[-]'),
               (mw + 34.9689, 'M+Cl[-]'),
               (mw + 36.9659, 'M+Cl37[-]'),
               (mw + 78.9183, 'M+Br[-]'),
               (mw + 80.9163, 'M+Br81[-]'),
               (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]'),
               (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]'),
               (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]'),
               (mw - PROTON + 15.99491, 'M-H+O[-]'),
               ]

    elif mode == 'neutral':
        print "Neutral mode of instrumentation is not supported."
        return []
    
    else:
        print "Unrecognized mode of instrumentation."
        return []

# weighting function of isotopic derivatives and adducts
dict_weight_adduct = {
            'M[1+]': 5, 
            'M+H[1+]': 5,
            'M+2H[2+]': 3,
            'M+3H[3+]': 1,
            'M(C13)+H[1+]': 2,
            'M(C13)+2H[2+]': 1,
            'M(C13)+3H[3+]': 1,
            'M(S34)+H[1+]': 1,
            'M(Cl37)+H[1+]': 1,
            'M+Na[1+]': 3, 
            'M+H+Na[2+]': 2,
            'M+K[1+]': 2, 
            'M+H2O+H[1+]': 1, 
            'M-H2O+H[1+]': 1, 
            'M-H4O2+H[1+]': 1,
            'M-NH3+H[1+]': 1,
            'M-CO+H[1+]': 1,
            'M-CO2+H[1+]': 1,
            'M-HCOOH+H[1+]': 1,
            'M+HCOONa[1+]': 1,
            'M-HCOONa+H[1+]': 1,
            'M+NaCl[1+]': 1, 
            'M-C3H4O2+H[1+]': 1,
            'M+HCOOK[1+]': 1,
            'M-HCOOK+H[1+]': 1,
            # negative
            'M-H[-]': 5,
            'M-2H[2-]': 3,
            'M(C13)-H[-]': 2,
            'M(S34)-H[-]': 1,
            'M(Cl37)-H[-]': 1,
            'M+Na-2H[-]': 2,
            'M+K-2H[-]': 1,
            'M-H2O-H[-]': 1,
            'M+Cl[-]': 1,
            'M+Cl37[-]': 1,
            'M+Br[-]': 1,
            'M+Br81[-]': 1,
            'M+ACN-H[-]': 1,
            'M+HCOO[-]': 1,
            'M+CH3COO[-]': 1,
            'M-H+O[-]': 1,
                    }




fishlogo = '''     
    --------------------------------------------
    
             oO                      ooooooooo
           oOO   OOOOO  ooooo       ooo oooo
     oOO   O       ooooo  oooooo ooooo
    oooO           oooooo         oooo ooooo
        Oooo   o      OOOOOO   oooo   oooooooo
            ooooo  oooo      
                 o
    
    --------------------------------------------
    '''










