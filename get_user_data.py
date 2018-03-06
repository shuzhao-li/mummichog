# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
Data input functions in mummichog

Read user input data from
    text files; xls files
    web form 
    future API support

Overall design in v2: separating user input from theoretical model.

@author: Shuzhao Li

'''
import time, getopt, urllib, base64
import logging
import StringIO

import matplotlib.pyplot as plt


from config import *
from models import *


logging.basicConfig(format='%(message)s', level=logging.INFO)

def print_and_loginfo(s):
    '''
    This function should retire?
    '''
    #print s
    logging.info(s)

#
# Functions to take command line input
#

def cli_options(opts):
    '''
    Ongoing work in version 2, making some options obsolete.
    
    obsolete parameters:
    'analysis': 'total',
    'targeted': False,
    'evidence': 3,
    'visualization': 2,
    
    '''
    time_stamp = str(time.time())
    
    optdict = {
               'cutoff': 0,
               
               'network': 'human_mfn',
               'modeling': None,
               
               'mode': 'pos_default',
               'instrument': 'unspecified',
               'force_primary_ion': True,
               
               'workdir': '',
               'input': '',
               'reference': '',
               'infile': '',
               'output': '',
               'permutation': 100,
               'outdir': 'mcgresult' + time_stamp,
               }
    booleandict = {'T': True, 'F': False, 1: True, 0: False, 
                   'True': True, 'False': False, 'TRUE': True, 'FALSE': False, 'true': True, 'false': False,
                    }
    modedict = {'default': 'pos_default', 'pos': 'pos_default', 'pos_default': 'pos_default',
                'dpj': 'dpj_positive', 'positive': 'generic_positive', 'Positive': 'generic_positive',
                'negative': 'negative', 'Negative': 'negative',
                    }
    # update default from user argument
    for o, a in opts:
        if o in ("-a", "--analysis"): optdict['analysis'] = a
        elif o in ("-c", "--cutoff"): optdict['cutoff'] = float(a)
        elif o in ("-t", "--targeted"): optdict['targeted'] = booleandict.get(a, False)
        elif o in ("-n", "--network"): optdict['network'] = a
        elif o in ("-z", "--force_primary_ion"): optdict['force_primary_ion'] = booleandict.get(a, True)
        elif o in ("-d", "--modeling"): optdict['modeling'] = a
        elif o in ("-e", "--evidence"): optdict['evidence'] = int(a)
        elif o in ("-m", "--mode"): optdict['mode'] = modedict.get(a, a)
        elif o in ("-u", "--instrument"): optdict['instrument'] = a
        elif o in ("-v", "--visualization"): optdict['visualization'] = int(a)
        elif o in ("-k", "--workdir"): optdict['workdir'] = a
        elif o in ("-i", "--input"): optdict['input'] = a
        elif o in ("-r", "--reference"): optdict['reference'] = a
        elif o in ("-f", "--infile"): optdict['infile'] = a
        elif o in ("-o", "--output"):
            optdict['output'] = a.replace('.csv', '')
            optdict['outdir'] = '.'.join([time_stamp, a.replace('.csv', '')])
            
        elif o in ("-p", "--permutation"): optdict['permutation'] = int(a)
        else: print "Unsupported argument ", o
    
    return optdict



def dispatcher():
    '''
    Dispatch command line arguments to corresponding functions.
    No user supplied id is used in version 1.
    User supplied IDs, str_mz_rtime IDs and targeted metabolites will be supported in version 2.
    

    '''
    helpstr = '''
    Usage example:
    python main.py -f mydata.txt -o myoutput
    
        -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score
        
        -n, --network: network model to use (default human_mfn), 
              [human, human_mfn, mouse, fly, yeast]
        
        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.
        
        -m, --mode: analytical mode of mass spec, [positive, negative, pos_defult].
              Default is pos_defult, a short version of positive.
        -u, --instrument: Any integer, treated as ppm of instrument accuracy. Default is 10. 
              
        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: one of primary ions, 
              ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]'],  
              must be present for a predicted metabolite, [True, False].
              Default is True.
        
        -c, --cutoff: optional cutoff p-value in user supplied statistics,
              used to select significant list of features. 
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is no.
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
# classes for data structure

class MassFeature:
    '''
    Data model, to store info per input feature
    row_number is used as unique ID. A string like "row23" is used instead of integer for two reasons:
    to enforce unique IDs in string not in number, and for clarity throughout the code;
    to have human friendly numbering, starting from 1 not 0.
    '''
    def __init__(self, row_number, mz, retention_time, p_value, statistic, CompoundID_from_user=''):
        self.row_number = row_number        # Unique ID
        self.mz = mz
        self.retention_time = retention_time
        self.p_value = p_value
        self.statistic = statistic
        self.CompoundID_from_user = CompoundID_from_user
        
        self.matched_Ions = []
        self.matched_Compounds = []
        self.matched_EmpiricalCompounds = []
        
        self.is_significant = False
        
        # for future use
        self.peak_quality = 0
        self.database_match = []

    def make_str_output(self):
        return '\t'.join( [str(x) for x in [self.row_number, self.mz, self.retention_time, 
                            self.p_value, self.statistic, self.CompoundID_from_user,
                            ]] ) 


class EmpiricalCompound:
    '''
    EmpiricalCompound is a computational unit to include 
    multiple ions that belong to the same metabolite,
    and isobaric/isomeric metabolites when not distinguished by the mass spec data.
    Thought to be a tentative metabolite. 
    Due to false matches, one Compound could have more EmpiricalCompounds
    
    In mummichog, this replaces the Mnode class in version 1;
    and is the compound presentation for Activity network and HTML report.
    
    This class serves as in between user-input MassFetaure and theoretical model Compound.
    
    
    '''
    def __init__(self, listOfFeatures):
        '''
        Initiation using 
        listOfFeatures = [[retention_time, row_number, ion, mass, compoundID], ...]
        This will be merged and split later to get final set of EmpCpds.
        '''
        self.listOfFeatures = listOfFeatures
        self.listOfFeatures.sort(key=lambda x: x[1])
        self.str_row_ion = self.__make_str_row_ion__()          # also a unique ID
        self.__unpack_listOfFeatures__()
        
        self.EID = ''
        self.chosen_compounds = []
        self.face_compound = ''
        
        self.evidence_score = 0
        self.primary_ion_present = False
        self.statistic = 0
    
    def __make_str_row_ion__(self):
        '''
        feature order is fixed now after sorting by row_number
        '''
        return ";".join([x[1]+'_'+x[2] for x in self.listOfFeatures])
    
    
    def __unpack_listOfFeatures__(self):
        
        self.compounds = list(set([x[4] for x in self.listOfFeatures]))
        self.massfeature_rows = [x[1] for x in self.listOfFeatures]
        self.ions = dict([x[2:4] for x in self.listOfFeatures])
        self.row_to_ion = dict( [x[1:3] for x in self.listOfFeatures] )
        
        
    def join(self, E):
        '''
        join another instance with identical ions.
        str_row_ion must be the same to be joined.
        '''
        for c in E.compounds:
            if c not in self.compounds:
                self.compounds.append(c)
        
        
    def evaluate(self):
        '''
        test if EmpCpds has any of primary_ions as defined in config.py, 
        ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]']
        
        evidence_score is combining weight scores from multiple adducts.
        
        No need to separate ionMode upfront ("positive", "negative")
        Bad ones should not be used for downstream analysis...
        
        '''
        
        if set(self.ions.keys()).intersection(primary_ions):
            self.primary_ion_present = True
        
        for x in self.ions.keys(): 
            self.evidence_score += dict_weight_adduct[x]
            
            
    def update_chosen_cpds(self, cpd):
        if cpd not in self.chosen_compounds:
            self.chosen_compounds.append(cpd)

    def designate_face_cpd(self):
        '''
        When there are more than one compounds suggested by pathway and module analysis,
        one is arbitrarily designated as "face compound".
        '''
        self.face_compound = self.chosen_compounds[-1]

    def get_mzFeature_of_highest_statistic(self, dict_mzFeature):
        '''
        Take highest abs(statistic) among all matched ions,
        which will give statistic value for downstream output.
        '''
        all = [dict_mzFeature[r] for r in self.massfeature_rows]
        all.sort(key=lambda m: abs(m.statistic))
        self.mzFeature_of_highest_statistic = all[-1]


class InputUserData:
    '''
    
    backward compatibility, 1 or 2-file input formats
    Per Joshua C., there'd be an option to test user designated L_sig, but user specified IDs are required
    
    return ListOfMassFeatures
    self.input_featurelist is "L_sig".
    '''
    
    def __init__(self, paradict):
        self.paradict = paradict
        self.header_fields = []
        self.ListOfMassFeatures = []
        self.input_featurelist = []
        
        self.read()
        self.determine_significant_list(self.ListOfMassFeatures)
        
        self.max_retention_time = max([M.retention_time for M in self.ListOfMassFeatures])
        self.max_mz = max([M.mz for M in self.ListOfMassFeatures])
        
        
    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user
        '''
        #
        lines = self.__check_redundant__( textValue.splitlines() )
        self.header_fields = lines[0].rstrip().split(delimiter)
        excluded_list = []
        for ii in range(len( lines )-1):
            y = lines[ii+1].split('\t')
            
            CompoundID_from_user = ''
            if len(y) > 4: CompoundID_from_user = y[4]
            [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
            
            # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
            if MASS_RANGE[0] < mz < MASS_RANGE[1]:
                # row # human-friendly, numbering from 1
                self.ListOfMassFeatures.append( 
                    MassFeature('row'+str(ii+1), mz, retention_time, p_value, statistic, CompoundID_from_user) )
            else:
                excluded_list.append( (ii, mz, retention_time) )
        
        if excluded_list:
            print_and_loginfo( "Excluding %d features out of m/z range %s." %(len(excluded_list), str(MASS_RANGE)) )

        
    def read_from_file(self, inputFile):
        return open(inputFile).read()
    
    def read_from_webform(self, t):
        return t

    def __check_redundant__(self, L):
        redundant = len(L) - len(set(L))
        if redundant > 0:
            print_and_loginfo( "Your input file contains %d redundant features." %(redundant) )
        return L

    def read(self):
        '''
        Read input feature lists to ListOfMassFeatures. 
        Row_numbers (rowii+1) are used as primary ID.
        # not using readlines() to avoid problem in processing some Mac files
        '''
        self.text_to_ListOfMassFeatures( 
                open(os.path.join(self.paradict['workdir'], self.paradict['infile'])).read() )
        print_and_loginfo("Read %d features as reference list." %len(self.ListOfMassFeatures))
    
    
    # more work?
    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot,
        in which case, paradict['cutoff'] is updated accordingly.
        
        '''
        if not self.paradict['cutoff']:
            # automated cutoff
            new = sorted(all_feature_list, key=lambda x: x.p_value)
            
            p_hotspots = [ 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 ]
            N_hotspots = [ len([x for x in all_feature_list if x.p_value < pp]) for pp in p_hotspots ]
            
            N_quantile = len(new) / 4
            N_optimum, N_minimum = 300, 30
            chosen = 9999
            for ii in range( len(N_hotspots) ):
                # will get the smallest p as ii increases
                if N_optimum < N_hotspots[ii] < N_quantile:
                    chosen = ii
            
            # if nothing was chosen
            if chosen > 100:
                for ii in range( len(N_hotspots) ):
                    if N_minimum < N_hotspots[ii] < N_quantile:
                        chosen = ii
            
            if chosen > 100:
                N_chosen = int(N_quantile)
                self.paradict['cutoff'] = new[N_chosen+1].p_value
            else:
                #N_chosen = N_hotspots[chosen]
                
                self.paradict['cutoff'] = p_hotspots[chosen]
        
            print_and_loginfo("Automatically choosing (p < %f) as significant cutoff."  %self.paradict['cutoff'])  
        
        # mark MassFeature significant
        for f in self.ListOfMassFeatures:
            if f.p_value < self.paradict['cutoff']:
                f.is_significant = True
        
        self.input_featurelist = [f.row_number for f in self.ListOfMassFeatures if f.is_significant]
        print_and_loginfo("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  


    def make_manhattan_plots(self, outfile='mcg_MWAS'):
        '''
        Manhattan plots of significance vs m/z, rtime.
        To use with reporting class
        '''
        
        # determine parameters
        figsize = (6, 3)
        CutoffLine = -np.log10(self.paradict['cutoff'])
        sigList = [ f for f in self.ListOfMassFeatures if f.p_value < self.paradict['cutoff'] ]
        restList = [ f for f in self.ListOfMassFeatures if f.p_value >= self.paradict['cutoff'] ]
        Y_label = "-log10 p-value"
        Y_black = [-np.log10(f.p_value) for f in restList]
        Y_green = [-np.log10(f.p_value) for f in sigList]
        X_label = ["m/z", "Retention time"]
        X_black = [ [f.mz for f in restList], [f.retention_time for f in restList] ]
        X_green = [ [f.mz for f in sigList], [f.retention_time for f in sigList] ]
        X_max = [self.max_mz, self.max_retention_time]
        
        # plot two panels, MWAS on m/z and rtime
        fig, myaxes = plt.subplots(figsize=figsize, nrows=1, ncols=2)
        for ii in range(2):
            
            myaxes[ii].scatter( X_black[ii], Y_black, s = 5, c='black', linewidths =0, alpha=0.8 )
            myaxes[ii].scatter( X_green[ii], Y_green, s=5, c='green', linewidths =0, alpha=0.8 )
            # lines
            myaxes[ii].plot([0,X_max[ii]], [CutoffLine, CutoffLine], 'g--')
        
            myaxes[ii].spines['right'].set_visible(True)
            myaxes[ii].spines['top'].set_visible(True)
            myaxes[ii].yaxis.set_ticks_position('both')
            myaxes[ii].xaxis.set_ticks_position('both')
            
            myaxes[ii].set_xlabel(X_label[ii])
            myaxes[ii].set_ylabel(Y_label)
            # rotate to avoid overlap xticklabels
            plt.setp(myaxes[ii].get_xticklabels(), rotation=30, horizontalalignment='right')
            
        #plt.title("Feature significance")
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')
        
        # get in-memory string for web use
        figdata = StringIO.StringIO()
        plt.savefig(figdata, format='png')
        figdata.seek(0)
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(figdata.buf))
        return '<img src = "%s"/>' % uri
        
        
        
        
        

# metabolicNetwork

class DataMeetModel:
    '''
    
    working on v2
    
    many to many matches:
    when a Compound matched to multiple MassFeatures, split by retention time to EmpiricalCompounds;
    when a Mass Feature matched to multiple Compounds, no need to do anything.
    
    
    Default primary ion is enforced, so that for an EmpiricalCompound, primary ion needs to exist before other ions.
    
    
    also here, compile cpd adduct lists, build cpd tree
    
    
    This returns the tracking map btw massFeatures - EmpiricalCompounds - Compounds
    Number of EmpiricalCompounds will be used to compute pathway enrichment, and control for module analysis.
    

    '''
    def __init__(self, theoreticalModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (M.row_number, EmpiricalCompounds, Cpd)
        
        '''
        self.model = theoreticalModel
        self.data = userData
        
        # retention time window for grouping, based on fraction of max rtime
        self.rtime_tolerance = self.data.max_retention_time * RETENTION_TIME_TOLERANCE_FRAC
        
        # major data structures
        self.IonCpdTree = self.__build_cpdindex__( self.data.paradict['mode'] )
        self.rowDict = self.__build_rowindex__( self.data.ListOfMassFeatures )
        self.ListOfEmpiricalCompounds = self.get_ListOfEmpiricalCompounds()
        
        # this is the reference list
        self.mzrows = [M.row_number for M in self.data.ListOfMassFeatures]
        
        self.rowindex_to_EmpiricalCompounds = self.__make_rowindex_to_EmpiricalCompounds__()
        self.Compounds_to_EmpiricalCompounds = self.__index_Compounds_to_EmpiricalCompounds__()
        
        # this is the sig list
        self.significant_features = self.data.input_featurelist
        self.TrioList = self.batch_rowindex_EmpCpd_Cpd( self.significant_features )
        
        
        

    def __build_cpdindex__(self, msmode):
        '''
        indexed Compound list, to speed up m/z matching.
        Limited to MASS_RANGE (default 50 ~ 2000 dalton).
        
        
        changing from adduct_function to wanted_adduct_list dictionary
        
        wanted_adduct_list['pos_default'] = ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],
        
        # 
        >>> metabolicModels['human_model_mfn']['Compounds'].items()[92]
        ('C00217', {'formula': '', 'mw': 147.0532, 'name': 'D-Glutamate; D-Glutamic acid; D-Glutaminic acid; D-2-Aminoglutaric acid',
         'adducts': {'M+2H[2+]': 74.53387646677, 'M+Br81[-]': 227.9695, 'M-H2O+H[1+]': 130.04987646677, 
         'M-C3H4O2+H[1+]': 76.03937646677, 'M-HCOOH+H[1+]': 102.05507646676999, 'M-HCOONa+H[1+]': 80.07307646677, 
         'M+K[1+]': 186.01597646677, 'M+Cl[-]': 182.0221, 'M+Na-2H[-]': 167.02064706646001, 'M-CO2+H[1+]': 104.07067646677, 
         'M+Na[1+]': 170.04247646677, 'M+Br[-]': 225.9715, 'M(S34)-H[-]': 148.04172353323, 'M+H[1+]': 148.06047646677, 
         'M-H4O2+H[1+]': 112.03927646677, 'M(C13)-H[-]': 147.04932353323, 'M(Cl37)-H[-]': 148.04312353323, 'M+HCOONa[1+]': 216.04787646677, 'M(C13)+2H[2+]': 75.03557646677, 'M+HCOOK[1+]': 232.02177646677, 'M-CO+H[1+]': 120.06547646677, 'M+HCOO[-]': 192.050845, 'M(C13)+3H[3+]': 50.359409800103336, 'M(Cl37)+H[1+]': 150.05767646677, 'M-H[-]': 146.04592353323, 'M+ACN-H[-]': 187.07246853323, 'M+Cl37[-]': 184.0191, 'M-H2O-H[-]': 128.03532353322998, 'M(S34)+H[1+]': 150.05627646677002, 'M-HCOOK+H[1+]': 64.09917646677, 'M+3H[3+]': 50.025009800103334, 'M+CH3COO[-]': 206.066495, 'M(C13)+H[1+]': 149.06387646677, 'M[1+]': 147.0532, 'M-NH3+H[1+]': 131.03397646677, 'M+NaCl[1+]': 206.01907646677, 'M+H+Na[2+]': 85.52487646677, 'M+H2O+H[1+]': 166.07107646677002, 'M-H+O[-]': 162.04083353323, 'M+K-2H[-]': 182.99414706646002, 'M-2H[2-]': 72.51932353323001}})
        >>> len(metabolicModels['human_model_mfn']['Compounds'])
        3560
        '''
        wanted_ions = wanted_adduct_list[msmode]
        IonCpdTree = []
        
        for ii in range(MASS_RANGE[1]+1): 
            IonCpdTree.append([])       #empty lists for anything below MASS_RANGE
            
        # iteritems vs items is contention of efficiency, but there's change btw Python 2 and Python 3...
        for c,d in self.model.Compounds.items():
            if d['mw']:                 #sanity check; bypass mistake in adducts type
                for ion,mass in d['adducts'].items():
                    if ion in wanted_ions and MASS_RANGE[0] < mass < MASS_RANGE[1]:
                        IonCpdTree[ int(mass) ].append( (c, ion, mass) )
                
        # tree: (compoundID, ion, mass), ion=match form; mass is theoretical
        return IonCpdTree


    def __build_rowindex__(self, ListOfMassFeatures):
        '''
        Index list of MassFeatures by row# in input data
        '''
        rowDict = {}
        for M in ListOfMassFeatures: rowDict[M.row_number] = M
        return rowDict


    def __match_all_to_all__(self):
        '''
        
        Major change of data structure here in version 2.
        In ver 1, matched m/z is stored in each Compound instance.
        Here, we produce mapping dictionaries for
            * mzFeatures to theoretical ions
            * Compounds to mzFeatures
        Then, 
            * EmpiricalCompounds are determined within Compound matched mzFeatures, considering retention time.
        
        
        '''
        self.__match_to_mzFeatures__()
        self.cpd2mzFeatures = self.index_Compounds_to_mzFeatures()
        return self.compound_to_EmpiricalCompounds()
        

    def __match_to_mzFeatures__(self):
        '''
        Fill mzFeatures with matched ions and compounds
        '''
        for M in self.data.ListOfMassFeatures:
            M.matched_Ions = self.__match_mz_ion__(M.mz, self.IonCpdTree)
        
        
    def index_Compounds_to_mzFeatures(self):
        '''
        compound ID - mzFeatures
        run after self.__match_to_mzFeatures__()
        L: (compoundID, ion, mass)
        cpd2mzFeatures[compoundID] = [(ion, mass, mzFeature), ...]
        '''
        cpd2mzFeatures = {}
        for M in self.data.ListOfMassFeatures:
            for L in M.matched_Ions:
                if cpd2mzFeatures.has_key(L[0]):
                    cpd2mzFeatures[L[0]].append( (L[1], L[2], M) )
                else:
                    cpd2mzFeatures[L[0]] = [(L[1], L[2], M)]
        
        print ("Got %d cpd2mzFeatures" %len(cpd2mzFeatures))
        return cpd2mzFeatures
        
        
    def __match_mz_ion__(self, mz, IonCpdTree):
        '''
        L: (compoundID, ion, mass)
        return ions matched to m/z
        '''
        floor = int(mz)
        matched = []
        mztol = mz_tolerance(mz, self.data.paradict['mode'])
        for ii in [floor-1, floor, floor+1]:
            for L in IonCpdTree[ii]:
                if abs(L[2]-mz) < mztol:
                    matched.append( L )
                    
        return matched

    def compound_to_EmpiricalCompounds(self):
        '''
        EmpiricalCompounds are constructed in this function.
        First splitting features matching to same compound by retention time;
        then merging those matched to same m/z features.
        run after self.index_Compounds_to_mzFeatures()
        '''
        ListOfEmpiricalCompounds = []
        for k,v in self.cpd2mzFeatures.items():
            ListOfEmpiricalCompounds += self.__split_Compound__(k, v, self.rtime_tolerance)      # getting inital instances of EmpiricalCompound
            
        print ("Got %d ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        
        # merge compounds that are not distinguished by analytical platform, e.g. isobaric
        return self.__merge_EmpiricalCompounds__( ListOfEmpiricalCompounds )
        
        
    def __split_Compound__(self, compoundID, list_match_mzFeatures, rtime_tolerance):
        '''
        Determine EmpiricalCompounds among the ions matched to a Compound;
        return list of EmpiricalCompounds (not final, but initiated here).
        
        The retention time is grouped by tolerance value; 
        This method should be updated in the future.
        
        input data format:
        cpd2mzFeatures[compoundID] = list_match_mzFeatures = [(ion, mass, mzFeature), ...]
        
        '''
        # unpacked format: [retention_time, row_number, ion, mass, compoundID]
        all_mzFeatures = [(L[2].retention_time, L[2].row_number, L[0], L[1], compoundID) for L in list_match_mzFeatures]
        all_mzFeatures.sort()
        ECompounds = []
        tmp = [ all_mzFeatures[0] ]
        for ii in range(len(all_mzFeatures)-1):
            if all_mzFeatures[ii+1][0]-all_mzFeatures[ii][0] < rtime_tolerance:
                tmp.append(
                            all_mzFeatures[ii+1] )
            else:
                ECompounds.append( EmpiricalCompound( tmp ) )
                tmp = [ all_mzFeatures[ii+1] ]
        
        ECompounds.append( EmpiricalCompound( tmp ) )
        return ECompounds
        
    
    def __merge_EmpiricalCompounds__(self, ListOfEmpiricalCompounds):
        '''
        If ion/mzFeatures are the same, merge EmpiricalCompounds
        EmpiricalCompounds.join() adds Compounds
        
        Because EmpiricalCompounds.str_row_ion uses mzFeatures sorted by row_number, this is 
        '''
        mydict = {}
        for L in ListOfEmpiricalCompounds:
            if mydict.has_key(L.str_row_ion):
                mydict[ L.str_row_ion ].join(L)
            else:
                mydict[ L.str_row_ion ]= L
        
        print ("Got %d merged ListOfEmpiricalCompounds" %len(mydict))
        return mydict.values()

    def __make_rowindex_to_EmpiricalCompounds__(self):
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.massfeature_rows:
                if mydict.has_key(m):
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict

    def __index_Compounds_to_EmpiricalCompounds__(self):
        '''
        Make dict cpd - EmpiricalCompounds
        '''
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.compounds:
                if mydict.has_key(m):
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict
        

    def batch_rowindex_EmpCpd_Cpd(self, list_features):
        '''
        Batch matching from row feature to Ecpds; Use trio data structure, (M.row_number, EmpiricalCompounds, Cpd).
        Will be used to map for both sig list and permutation lists.
        '''
        new = []
        for f in list_features:
            for E in self.rowindex_to_EmpiricalCompounds.get(f, []):
                for cpd in E.compounds:
                    new.append((f, E, cpd))
            
        return new

            
    def get_ListOfEmpiricalCompounds(self):
        '''
        Collect EmpiricalCompounds.
        Initiate EmpCpd attributes here.
        '''
        ListOfEmpiricalCompounds, ii = [], 1
        for EmpCpd in self.__match_all_to_all__():
            EmpCpd.evaluate()
            EmpCpd.EID = 'E' + str(ii)
            EmpCpd.get_mzFeature_of_highest_statistic( self.rowDict )
            ii += 1
            if self.data.paradict['force_primary_ion']:
                if EmpCpd.primary_ion_present:
                    ListOfEmpiricalCompounds.append(EmpCpd)
            else:
                ListOfEmpiricalCompounds.append(EmpCpd)
        
        print ("Got %d final ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        return ListOfEmpiricalCompounds

        