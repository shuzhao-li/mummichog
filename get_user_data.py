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

import logging
from config import *
from models import *


logging.basicConfig(format='%(message)s', level=logging.INFO)

def print_and_loginfo(s):
    print s
    logging.info(s)




class MassFeature:
    '''
    Data model, to store info per input feature
    
    '''
    def __init__(self, row_number, mz, retention_time, p_value, statistic, CompoundID_from_user=''):
        self.row_number = row_number
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
        return '\t'.join( [str(x) for x in [elf.row_number, self.mz, self.retention_time, 
                            self.p_value, self.statistic, self.CompoundID_from_user,
                            ','.join(self.matched_compounds),
                            ]] ) 

    def match_compounds(self):
        return 0





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
        
        compound IDs, row # for massFeatures
        
        listOfFeatures = [[retention_time, row_number, ion, mass, compoundID], ...]
        
        
        '''
        self.listOfFeatures = listOfFeatures
        self.str_row_ion = self.make_str_row_ion()
        self.unpack_listOfFeatures()
                      
        self.mostLikelyCompounds = []
        # mostLikelyCompound will be designated by module/pathway analysis
        self.evidence_score = 0
        self.chosen = None
        self.primary_ion_present = False
        #self.evaluate()
    
    
    def make_str_row_ion(self):
        return ";".join([x[1]+x[2] for x in self.listOfFeatures])
    
    
    
    def unpack_listOfFeatures(self):
        
        self.compounds = list(set([x[4] for x in self.listOfFeatures]))
        self.massfeature_rows = [x[1] for x in self.listOfFeatures]
        self.ions = dict([x[3:5] for x in self.listOfFeatures])
        
        
        
    def join(self, E):
        '''
        join another instance with identical ions
        '''
        for c in E.compounds:
            if c not in self.compounds:
                self.compounds.append(c)
        
        
    def evaluate(self):
        '''
        Compute evidence scores using a weight dictionary. The logic goes as -
        If 'M+H[1+]' present, 
        any of ('M+2H[2+]', 'M(C13)+H[1+]', 'M+Na[1+]', 'M+K[1+]', 'M+2Na[2+]') confirms.
        If 'M+2H[2+]' present,
        'M(C13)+2H[2+]' confirms.
        Compound.hitlist = [(input mz, match form, diff), ...] --- Changed!!!
        [( MassFeature instance, match form, diff), ...]
        
        
        
        to redo
        '''
        matchforms = set([x[1] for x in self.hitlist])
        if 'M+H[1+]' in matchforms or 'M-H[-]' in matchforms:
            self.primary_ion_present = True
        for x in matchforms: 
            self.evidence_score += dict_weight_adduct[x]
        





class InputUserData:
    '''
    
    backward compatibility, 1 or 2-file input formats
    Per Joshua, there'd be an option to test user designated L_sig, but user specified IDs are required
    
    return ListOfMassFeatures
    self.input_featurelist is "L_sig".
    
    Targeted data should use a separate workflow.
    
    '''
    
    def __init__(self, paradict):
        self.paradict = paradict
        self.header_fields = []
        self.ListOfMassFeatures = []
        self.read()
        self.determine_significant_list(self.ListOfMassFeatures)
        
        
    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user
        '''
        #
        lines = self.check_redundant( textValue.splitlines() )
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

    def check_redundant(self, L):
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
    
    
    
    # in work
    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot.
        
        
        '''
        if self.paradict['cutoff']:
            # use user specified cutoff
            self.input_featurelist = [x.row_number for x in all_feature_list if x.p_value < self.paradict['cutoff']]
        else:
            # automated cutoff
            new = [(x.p_value, x) for x in all_feature_list]
            new.sort()
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
                print_and_loginfo("Using %d top quantile of features as significant list." %N_chosen)
            else:
                N_chosen = N_hotspots[chosen]
                print_and_loginfo("Using %d features (p < %f) as significant list." %(N_chosen, p_hotspots[chosen]))
        
            self.input_featurelist = [x[1].row_number for x in new[:N_chosen]]
        



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
    Number of EmpiricalCompounds will be used to compute enrichment.
    

    '''
    def __init__(self, theoreticalModel, userData):
        '''
        working on v2
        # from ver 1 to ver 2, major change in .match()
        
        Important data indices:
        self.
        rowDict
        cpd2mzFeatures
        ListOfEmpiricalCompounds
        '''
        
        self.model = theoreticalModel
        self.data = userData
        self.build_cpdindex( self.data.paradict['mode'] )
        self.build_rowindex( self.data.ListOfMassFeatures )
        self.mzrows = [M.row_number for M in self.data.ListOfMassFeatures]
        
        self.match_all_to_all()
        self.rowindex_to_EmpiricalCompounds = self.make_rowindex_to_EmpiricalCompounds()
        self.Compounds_to_EmpiricalCompounds = self.index_Compounds_to_EmpiricalCompounds()
        self.rowindex_to_Compounds = self.make_rowindex_to_Compounds()
        
        self.significant_features = self.data.input_featurelist
        self.significant_EmpiricalCompounds = self.compile_significant_EmpiricalCompounds(self.significant_features)
        
        self.hit_EmpiricalCompounds = []



    def build_cpdindex(self, msmode):
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
        self.IonCpdTree = IonCpdTree


    def build_rowindex(self, ListOfMassFeatures):
        '''
        Index list of MassFeatures by row# in input data
        '''
        self.rowDict = {}
        for M in ListOfMassFeatures:
            self.rowDict[M.row_number] = M


    def match_all_to_all(self):
        '''
        
        Major change of data structure here in version 2.
        In ver 1, matched m/z is stored in each Compound instance.
        Here, we produce mapping dictionaries for
            * mzFeatures to theoretical ions
            * Compounds to mzFeatures
        Then, 
            * EmpiricalCompounds are determined within Compound matched mzFeatures, considering retention time.
        
        
        '''
        self.match_to_mzFeatures()
        self.index_Compounds_to_mzFeatures()
        self.compound_to_EmpiricalCompounds()
        

    def match_to_mzFeatures(self):
        '''
        Fill mzFeatures with matched ions and compounds
        '''
        for M in self.data.ListOfMassFeatures:
            M.matched_Ions = self.match_mz_ion(M.mz, self.IonCpdTree)
        
        
    def index_Compounds_to_mzFeatures(self):
        '''
        compound ID - mzFeatures
        run after self.match_to_mzFeatures()
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
        
        self.cpd2mzFeatures = cpd2mzFeatures
        print ("Got %d cpd2mzFeatures" %len(cpd2mzFeatures))
        
    def match_mz_ion(self, mz, IonCpdTree):
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

    def compound_to_EmpiricalCompounds(self, rtime_tolerance = 60):
        '''
        
        run after self.index_Compounds_to_mzFeatures()
        '''
        ListOfEmpiricalCompounds = []
        for k,v in self.cpd2mzFeatures.items():
            ListOfEmpiricalCompounds += self.split_Compound(k, v, rtime_tolerance)
            
        print ("Got %d ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        # merge compounds that are not distinguished by analytical platform, e.g. isobaric
        self.merge_EmpiricalCompounds( ListOfEmpiricalCompounds )
        print ("Got %d merged ListOfEmpiricalCompounds" %len(self.ListOfEmpiricalCompounds))
        
    def split_Compound(self, compoundID, list_match_mzFeatures, rtime_tolerance):
        '''
        Determine EmpiricalCompounds among the ions matched to a Compound;
        return list of EmpiricalCompounds
        
        The retention time is grouped by tolerance value; 
        This method should be updated in the future
        
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
        
    
    def merge_EmpiricalCompounds(self, ListOfEmpiricalCompounds):
        '''
        If ion/mzFeatures are the same, merge EmpiricalCompounds
        '''
        mydict = {}
        for L in ListOfEmpiricalCompounds:
            if mydict.has_key(L.str_row_ion):
                mydict[ L.str_row_ion ].join(L)
            else:
                mydict[ L.str_row_ion ]= L
            
        self.ListOfEmpiricalCompounds = mydict.values()


    def make_rowindex_to_EmpiricalCompounds(self):
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.massfeature_rows:
                if mydict.has_key(m):
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict


    def make_rowindex_to_Compounds(self):
        '''
        mapping through ListOfEmpiricalCompounds. This is a shortcut. Use with caution.
        returned dict may be redundant.
        '''
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.massfeature_rows:
                if mydict.has_key(m):
                    mydict[m] += E.compounds
                else:
                    mydict[m] = E.compounds
                    
        return mydict


    def index_Compounds_to_EmpiricalCompounds(self):
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
        

    def compile_significant_EmpiricalCompounds(self, significant_features):
        '''
        This function will be used to map for permutation lists too.
        '''
        EmCpds = []
        for f in significant_features:
            EmCpds += self.rowindex_to_EmpiricalCompounds.get(f, [])
        return set(EmCpds)
        

    def export_userData(self, outfile):
        '''
        to do
        
        Should add this function to write out user data with row_numbers,
        to help users to track data.
        
        '''
        with open(outfile, 'w') as O:
            O.write('\n'.join([ 
                    ]))














