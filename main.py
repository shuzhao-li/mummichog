#!/usr/bin/env python
# Copyright (c) 2010-2015 Shuzhao Li.
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
mummichog 2,
pathway and network analysis for untargeted and targeted metabolomics

@author: Shuzhao Li, Andrei Todor

online documentation:
http://code.google.com/p/atcg/wiki/mummichog_for_metabolomics
to be moved to clinicalmetabolomics.org



Working to version 2.
Add support of targeted data.
Enforce retention time match for featruers in the same Compound class.




Adding
-t, -- targeted: [False, True] to switch btw targeted vs untargeted data

Removing
-i, -r

'''

import time, getopt, logging
from config import *

def cli_options(opts):
    optdict = {'analysis': 'total',
               'cutoff': 0,
               'targeted': False,
               'network': 'human_mfn',
               'modeling': 'gamma',
               'evidence': 3,
               'mode': 'dpj',
               'instrument': 'unspecified',
               'force_primary_ion': False,
               'visualization': 2,
               'workdir': '',
               'input': '',
               'reference': '',
               'infile': '',
               'output': '',
               'permutation': 100,
               'outdir': 'mcgresult',
               }
    booleandict = {'T': True, 'F': False, 1: True, 0: False, 
                   'True': True, 'False': False, 'TRUE': True, 'FALSE': False, 'true': True, 'false': False,
                   }
    # update default from user argument
    for o, a in opts:
        if o in ("-a", "--analysis"): optdict['analysis'] = a
        elif o in ("-c", "--cutoff"): optdict['cutoff'] = float(a)
        elif o in ("-t", "--targeted"): optdict['targeted'] = booleandict.get(a, False)
        elif o in ("-n", "--network"): optdict['network'] = a
        elif o in ("-z", "--force_primary_ion"): optdict['force_primary_ion'] = booleandict.get(a, False)
        elif o in ("-d", "--modeling"): optdict['modeling'] = a
        elif o in ("-e", "--evidence"): optdict['evidence'] = int(a)
        elif o in ("-m", "--mode"): optdict['mode'] = a
        elif o in ("-u", "--instrument"): optdict['instrument'] = a
        elif o in ("-v", "--visualization"): optdict['visualization'] = int(a)
        elif o in ("-k", "--workdir"): optdict['workdir'] = a
        elif o in ("-i", "--input"): optdict['input'] = a
        elif o in ("-r", "--reference"): optdict['reference'] = a
        elif o in ("-f", "--infile"): optdict['infile'] = a
        elif o in ("-o", "--output"):
            optdict['output'] = a.replace('.csv', '')
            optdict['outdir'] = '.'.join([str(time.time()), a.replace('.csv', '')])
            
        elif o in ("-p", "--permutation"): optdict['permutation'] = int(a)
        else: print "Unsupported argument ", o
    
    return optdict



def dispatcher():
    '''
    Dispatch command line arguments to corresponding functions.
    No user supplied id is used in version 1.
    User supplied IDs, str_mz_rtime IDs and targeted metabolites will be supported in version 2.
    
    To add support for user cutoff of significant list
    
    
    
    '''
    helpstr = '''
    Usage example:
    python main.py -f mydata.txt -o myoutput
    
        -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score
        
        -t, --targeted: set to True if using targeted metabolomics data
        -n, --network: network model to use (default human_mfn), 
              [human, human_mfn, mouse, fly, yeast]
        
        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.
        
        -m, --mode: analytical mode of mass spec, [positive, negative, dpj].
              Default is dpj, a short version of positive.
        -u, --instrument: [5, 10, 25, FTMS, ORBITRAP].
              Any integer is treated as ppm. Default is 10. 
              Instrument specific functions may be implemented.
              
        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: M+H[+] (M-H[-] for negative mode) must be 
              present for a predicted metabolite, [True, False].
              Default is False.
        
        -c, --cutoff: optional cutoff p-value in user supplied statistics,
              used to select significant list of features. 
        -e, --evidence: cutoff score for metabolite to be in activity network.
              Default is 3.
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is gamma.
        '''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:t:d:e:m:n:u:z:v:k:i:r:f:o:p:", 
                            ["analysis=", "cutoff", "targeted=", "modeling=", "evidence=", "mode=", 
                             "network=", "instrument=", "force_primary_ion",
                             "visualization=", "workdir=", "input=", 
                             "reference=", "infile=", "output=", "permutation="])
        if not opts:
            print helpstr
            sys.exit(2)
        
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    return cli_options(opts)
    
#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    print fishlogo
    print "mummichog version %s \n" %VERSION
    optdict = dispatcher()
    
    # locate current path, for safe pydata import
    cmd_folder = os.path.split(inspect.getfile( inspect.currentframe() ))[0]
    if cmd_folder not in sys.path: sys.path.insert(0, cmd_folder)
    
    # load specific metabolic model
    if optdict['network'] in ['human', 'hsa', 'Human']:
        import pydata.human_model_humancyc as MetabolicModel
    
    elif optdict['network'] in ['human_mfn', 'hsa_mfn',]:
        import pydata.human_model_mfn as MetabolicModel
        
    elif optdict['network'] in ['mouse', 'Mouse',]:
        import pydata.mouse_model_biocyc as MetabolicModel
        
    elif optdict['network'] in ['fly', 'Fly',]:
        import pydata.fly_model_biocyc as MetabolicModel
        
    elif optdict['network'] in ['yeast', 'Yeast',]:
        import pydata.yeast_model_biocyc as MetabolicModel
        
    else:
        raise KeyError( "Unsupported species/model. Pls contact author." )
        
    
    # prepare output
    os.mkdir(os.path.join(optdict['workdir'], optdict['outdir']))
    
    logging.basicConfig(filename=os.path.join(optdict['workdir'], optdict['outdir'], 'mummichog.log'), 
                        format='%(message)s', 
                        level=logging.INFO)
    logging.info('\n'.join(["mummichog version: %s" %VERSION,
                             "pwd: %s" %os.getcwd(),
                             "user command: %s" %' '.join(sys.argv),
                             "\n",
                             ]))
    
    from functional_analysis import *

    print_and_loginfo("Started @ %s\n" %time.asctime())
    
    HNET = HsaNetwork(MetabolicModel, optdict)
    
    # if not optdict['targeted']:
    # Untargeted analysis
    AC = AnalysisCentral(HNET, HNET.input_mzlist)
    AC.create_dirs()
    AC.run_all_analysis()
    AC.export_csv_data()
    AC.export_sif_related()
    AC.web_export()
    
    # else
    # Targeted analysis
    # may not necessary to switch here - use same flow as above?
    
    
    
    print_and_loginfo("\nFinished @ %s\n" %time.asctime())

