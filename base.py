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
Basic data classes in mummichog

@author: Shuzhao Li, Andrei Todor
'''

import logging, csv
import networkx as nx

from config import *


logging.basicConfig(format='%(message)s', level=logging.INFO)

def print_and_loginfo(s):
    print s
    logging.info(s)


class Compound(object):
    '''
    This is a theoretical compound.
    A Compound instance is created each time of analysis,
    as the instrument operating mode is set per data set,
    and the instance also stores hits from m/z feature table.
    
    ver 2
    Adding retention time to consideration in mapping
    
    
    The two compound classes can incorporate more knowledge in the future.
    
    '''
    def __init__(self, id_str, molecular_weight, ms_mode, MetabolicModel):
        '''
        ms_mode designates the mode of MS operation.
        '''
        self.id = id_str
        self.name = MetabolicModel.dict_cpds_def.get(id_str, '')
        self.mw = molecular_weight
        self.mzlist = self.make_mzlist(mode=ms_mode)
        self.hitlist = []
        
    def make_mzlist(self, mode='dpj_positive'):
        '''
        return list of possible derivatives/adducts.
        mzlist is produced as [(mw, 'M[1+]'), (mw + PROTON, 'M+H[1+]'), ...]
        
        To do: 
        control adduct list by chemical formula
        
        '''
        if mode in ('dpj_positive', 'dpj', 'DPJ'):
            return adduct_function(self.mw, 'dpj_positive')
        
        elif mode in ('positive', 'generic_positive', 'generic', 'POS'):
            return adduct_function(self.mw, 'generic_positive')
            
        elif mode in ('neg', 'negative', 'NEG'):
            return adduct_function(self.mw, 'negative')

        else:
            return adduct_function(self.mw, mode)

    def match(self, mass_feature, mz_tolerance):
        '''
        record matched m/z to self.hitlist as
        [(input mz, match form, diff, retention_time), ...]
        
        use MassFeature instances;
        
        This still uses mz as primary ID; will change to row_number
        
        
        
        Change of hitlist data structure!!!
        
        
        '''
        for x in self.mzlist:
            diff = mass_feature.mz - x[0]
            if abs(diff) < tolerance:
                #self.hitlist.append( (mz, x[1], round(diff, 4), rtime) )
                self.hitlist.append( (mass_feature, x[1], round(diff, 4)) )
                
    def split(self, rtime_tolerance = 60):
        '''
        group and split hits by retention time, to EmoiricalCompound instances
        
        '''
        all_rtimes = [x[0].retention_time for x in self.hitlist]
        # now split the hits by retention_time
        
        # to finish
        
        
        
        # return [EmoiricalCompound instances, ...]
        return []
        
        
        
        



class EmpiricalCompound:
    '''
    
    EmpiricalCompound replacing the Mnode class -
    
    An Mnode is a node unit in the inspected network, 
    extension to Compound class to include evaluation of m/z patterns.
    
    
    # below could be dated
    It's equivalent to a tentative metabolite.
    It is a collection of m/z features that may belong to the same metabolite,
    reflecting isotopic forms and derivatives,
    and suggesting the most probable metabolite candidate.
    This is a unit for computational analysis and modeling, with room for error, 
    especially in case of metabolites of the same molecular weight,
    where a mixture is very probable.
    
    1. one instance per theoretical MW
    2. aggregate evidence of derivatives
       'M[1+]' AND 'M+2H[2+]' = True
       'M[1+]' AND 'M+Na[1+]' = True
       ....
    3. if cpds use same m/z values, choose cpd with stronger evidence from 2.
       if not resolved, use cpd with greatest degree in inspected network
    '''
    def __init__(self, cpd_id, cpd_mw):
        '''
        
        
        
        '''
        #self.cpd = CPD
        self.id = self.cpd_id
        self.mwstr = str(round(cpd_mw, 4))
        self.hitlist = []                   # [( MassFeature instance, match form, diff), ...]
        
        self.evidence_score = 0
        self.chosen = None
        self.primary_ion_present = False
        self.evaluate()
        
    def evaluate(self):
        '''
        Compute evidence scores using a weight dictionary. The logic goes as -
        If 'M+H[1+]' present, 
        any of ('M+2H[2+]', 'M(C13)+H[1+]', 'M+Na[1+]', 'M+K[1+]', 'M+2Na[2+]') confirms.
        If 'M+2H[2+]' present,
        'M(C13)+2H[2+]' confirms.
        Compound.hitlist = [(input mz, match form, diff), ...] --- Changed!!!
        [( MassFeature instance, match form, diff), ...]
        '''
        matchforms = set([x[1] for x in self.hitlist])
        if 'M+H[1+]' in matchforms or 'M-H[-]' in matchforms:
            self.primary_ion_present = True
            
        for x in matchforms: 
            self.evidence_score += dict_weight_adduct[x]
        




class MassFeature():
    '''
    Short hand for data model, to store info per input feature
    
    
    
    
    '''
    def __init__(self, row_number, mz, retention_time, p_value, statistic, CompoundID_from_user):
        self.row_number = row_number
        self.mz = mz
        self.retention_time = retention_time
        self.p_value = p_value
        self.statistic = statistic
        self.CompoundID_from_user = CompoundID_from_user
        
        self.is_significant = False
        
        # for future use
        self.peak_quality = 0
        self.database_match = []

    def make_str_output(self):
        pass


class Pathway:
    def __init__(self):
        self.id = ''
        self.name = ''
        self.rxns = []
        self.ecs = []
        self.ec_num = 0
        self.cpds = []
        self.cpd_num = 0
        
        self.selected_features = []
        self.combined_pvalue = 0    # -log10()
        
    def str_import(self, s):
        '''
        function to import pathway from plain text;
        excluding currency metabolites.
        Not used for now.
        '''
        a = s.rstrip().split('\t')
        self.id = a[0]
        self.name = a[1]
        self.rxns = a[2].split(';')
        self.ecs = a[3].split(';')
        self.ec_num = int(a[4])
        cpds = a[5].split(';')
        self.cpds = [x for x in cpds if x not in currency]
        self.cpd_num = len(self.cpds)
        
    def json_import(self, j):
        '''
        function to import pathway from JSON format.
        '''
        self.id = j['id']
        self.name = j['name']
        self.rxns = j['rxns']
        self.ecs = j['ecs']
        self.ec_num = len(j['ecs'])
        self.cpds = j['cpds']
        self.cpd_num = len(j['cpds'])



class HsaNetwork:
    '''
    Central data structure, containing
    1) a metabolites-centric metabolic model (not limited to Hsa; legacy class name)
    2) utility of processing user data
    3) record of m/z to cpd matching, which will be used by subsequent functions
    
    As the class is initiated for each new analysis, 
    analysis-wide user data are attached to its instance.
    
    
    Phased out 2-file input format.
    Support both targeted data (5th column as user supplied CompoundID)
    and untargeted data
    
    Will change class name to MetabolicNetwork
    
    Will change primary ID from m/z to row_number
    
    '''
    def __init__(self, MetabolicModel, paradict):
        '''
        Initiation of metabolic network model.
        Building Compound index.
        Parsing input files.
        Matching m/z - Compound.
        '''
        print_and_loginfo( "Loading metabolic network %s..." %MetabolicModel.version ) # version from metabolic model
        self.paradict = paradict
        self.MetabolicModel = MetabolicModel
        self.network = self.build_network(MetabolicModel.cpd_edges)
        self.total_cpd_list = self.network.nodes()
        
        self.cpd_dict = {}                  # id to Compound
        self.build_cpdindex(paradict['mode'])
                                            # to make CpdList and CpdTree
        self.input_mzlist = []
        self.ref_mzlist = []
        self.input_mzfcdict = {}            # fc can be any statistic from input file
        self.read()
        
        self.total_matched_cpds = []
        self.match_dict = {}                # mz to List of Compounds
        
        #
        # from ver 1 to ver 2, major change in .match()
        #
        self.match()
        
        self.significant_cpdlist = []
        
        
    def build_network(self, edges):
        return nx.from_edgelist( edges )
        
    def build_cpdindex(self, msmode):
        '''
        indexed Compound list, to speed up m/z matching.
        Limited to 50 ~ 2000 dalton.
        '''
        CpdList = []
        for k,v in self.MetabolicModel.dict_cpds_mass.items():
            if v > 0 and k in self.total_cpd_list:
                CpdList.append( Compound(k, v, msmode, self.MetabolicModel) )
                
        print_and_loginfo("cpds with MW: %d" %len(CpdList))
        
        CpdTree = {}
        for ii in range(49, 2001): CpdTree[ii] = []
            
        for C in CpdList:
            self.cpd_dict[C.id] = C
            neighborhood = [int(x[0]) for x in C.mzlist]
            for n in neighborhood:
                if 50 < n < 2001:
                    CpdTree[n].append(C)

        self.CpdList = CpdList
        self.CpdTree = CpdTree


    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot.
        
        
        '''
        if self.paradict['cutoff']:
            # use user specified cutoff
            self.input_mzlist = [x.mz for x in all_feature_list if x.p_value < self.paradict['cutoff']]
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
        
            self.input_mzlist = [x[1].mz for x in new[:N_chosen]]
        


    def read(self):
        '''
        Read input feature lists. 
        The significant list, input_mzlist, should be a subset of ref_mzlist.
        Checkpoint is added to limited m/z range (50, 2000).
        Made mz - fold change dictionary.
        Fold change is overwritten if redundant m/z entries are encountered.
        
        Adding new input format as single file...
        which contains a header, columns of m/z, rtime, p-value, statistic, [CompoundID].
        
        Cutoff for significant list is automatically determined in this function,
        by choosing a hotspot of p-value, or top quantile.
        
        This populates self.input_mzlist, self.ref_mzlist.
        The `statistic` column is used for `fold change` coloring.
        
        
        Primary ID for analysis will be row_number from input file;
        if user supplies CompoundID, CompoundID will be used for targeted analysis.
        
        
        
        If targeted analysis, limit to features with user supplied CompoundID.
        
        
        
        
        if user specifies p-value cutoff:
        
        
        error_message = "\nNot all significant features are in your reference list"

        
        if self.paradict['infile']:
            self.read_single_file_format()
        else:
            # self.read_two_file_format()
            print "Version 2 only supports one-file input."
            break
        
        # check m/z range (50, 2000); filter if necessary
        # This can be moved to read() as ver 2 now uses a single input file
        IL, RL = len(self.input_mzlist), len(self.ref_mzlist)
        self.input_mzlist = set([x for x in self.input_mzlist if 50 < x < 2001])
        self.ref_mzlist = set([x for x in self.ref_mzlist if 50 < x < 2001])
        IL2, RL2 = len(self.input_mzlist), len(self.ref_mzlist)
        if IL != IL2 or RL != RL2:
            print_and_loginfo( range_warning )
            
        print_and_loginfo("Got %d significant features from %d references" %(IL2, RL2))
        
        if IL2 > 2000 or RL2 > 20000:
            print size_warning
            logging.warning(size_warning)
        
        '''
        
        # not using readlines() to avoid problem in processing some Mac files
        rows = open(os.path.join(self.paradict['workdir'], 
                   self.paradict['infile'])).read().splitlines()    
        
        all_feature_list, excluded_list = [], []
        for ii in range(1, len(rows)):
            # skipping first row
            if rows[ii][0] != '#':
                y = rows[ii].rstrip().split('\t')
                CompoundID_from_user = ''
                if len(y) > 4: CompoundID_from_user = y[4]
                [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
                
                # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
                if 50 < mz < 2000:
                    all_feature_list.append( 
                        MassFeature(ii, mz, retention_time, p_value, statistic, CompoundID_from_user) )
                else:
                    excluded_list.append( (ii, mz, retention_time) )
        
        if excluded_list:
            print_and_loginfo("Excluding %d features out of m/z range (50, 2000)." %len(excluded_list))
        
        # if untargeted
        if not self.paradict['targeted']:
            # determining and populate significant list of features to use
            self.determine_significant_list( all_feature_list )
            
            for x in all_feature_list:
                self.ref_mzlist.append( x.mz )
                self.input_mzfcdict[ x.mz ] = x.statistic
            
        # if targeted, use CompoundID_from_user as primary ID
        # only supports KEGG style ID for now, as in mfn_default model
        # user has to supply cutoff p_value
        #
        #
        # Yet to see if this can be combined with pathway and module analysis
        #
        #
        else:
            all_feature_list = [x for x in all_feature_list if x.CompoundID_from_user]
            self.input_mzlist = [x.CompoundID_from_user for x in all_feature_list 
                        if x.p_value < self.paradict['cutoff']]
            self.ref_mzlist = [x.CompoundID_from_user for x in all_feature_list]
            for x in all_feature_list: self.input_mzfcdict[ x.CompoundID_from_user ] = x.statistic
            
        # will phase out when primary ID changes to row_number
        self.input_mzlist, self.ref_mzlist = set(self.input_mzlist), set(self.ref_mzlist)
        self.all_feature_list = all_feature_list
        

    def match(self):
        '''
        
        Major change of data structure here
        
        
        
        In ver 1, matched m/z is stored in each Compound instance.
        
        to split Compound to EmpiricalCompound,
        by retention time, ...
        
        so EmpiricalCompound is used for match_dict
        
        '''
        for x in self.all_feature_list:
            self.match_dict[x.mz] = []
            mztol = mz_tolerance(x.mz, self.paradict['mode'])
            floor = int(x.mz)
            for ii in [floor-1, floor, floor+1]:
                for C in self.CpdTree[ii]:
                    C.match(x, mztol)
                    
        for C in self.CpdList:
            # here to split Compound to EmpiricalCompound
            split_EmpiricalCompounds = C.split()
            
            for CC in split_EmpiricalCompounds:
                for x in CC.hitlist:
                    self.match_dict[x].append(CC)            # mz to EmpiricalCompound
                    self.total_matched_cpds.append(CC.id)
            
            
        # this uses cpd ID
        self.total_matched_cpds = set(self.total_matched_cpds)

















class Mmodule:
    '''
    Metabolites by their connection in metabolic network.
    A module is a subgraph, while modularity is calculated in 
    the background of reference hsanet.
    '''
    def __init__(self, hsanet, subgraph, seed_cpds, TF):
        self.network = hsanet.network
        self.num_ref_edges = hsanet.network.number_of_edges()
        self.num_ref_nodes = hsanet.network.number_of_nodes()
        self.graph = subgraph
        self.shave(seed_cpds)
        
        self.nodestr = self.make_nodestr()
        self.A = self.activity_score(seed_cpds, TF)
    
    def activity_score(self, seed_cpds, TF):
        '''
        A * (Ns/Nm)
        A = Newman-Girvan modularity score
        Ns = number of input cpds in module M
        Nm = number of total cpds in M
        Ns/Nm can be corrected as (Ns/total input size)/(Nm/network size), however,
        this normalization factor holds the same in permutations and can be left out.
        
        To reduce bias towards larger modules in Q:
        np.sqrt(self.num_ref_nodes/Nm) * 
        '''
        Ns = [x for x in self.graph.nodes() if x in seed_cpds]
        Ns = min(len(Ns), TF.count_cpd2mz(Ns))
        Nm = float(self.graph.number_of_nodes())
        self.compute_modularity()
        return np.sqrt(len(TF.input_cpdlist)/Nm) *self.Q * (Ns/Nm)
        
        
    def compute_modularity(self):
        '''
        To compute Newman-Girvan modularity for a single module,
        in reference to the whole network.
        '''
        m = self.num_ref_edges
        Nodes = self.graph.nodes()
        expected = 0
        for ii in Nodes:
            for jj in Nodes:
                if ii != jj:
                    expected += self.network.degree(ii) * self.network.degree(jj)
                    
        expected /= (4.0 * m)
        self.Q = (self.graph.number_of_edges() - expected) / m
    
    def test_compute_modularity(self):
        '''
        Alternative modularity measure as 
        edges in module over all edges on the same nodes
        '''
        m = float(self.graph.number_of_edges())
        expected = 0
        for ii in self.graph.nodes(): expected += self.network.degree(ii)
        self.Q = 2 * m * (np.sqrt(self.graph.number_of_nodes())) / expected
        
        
    def shave(self, seed_cpds):
        '''
        shave off nodes that do not connect seeds, i.e.
        any node with degree = 1 and is not a seed, iteratively.
        '''
        nonseeds = [x for x in self.graph.nodes() if x not in seed_cpds]
        excessive = [x for x in nonseeds if self.graph.degree(x)==1]
        while excessive:
            for x in excessive: self.graph.remove_node(x)
            nonseeds = [x for x in self.graph.nodes() if x not in seed_cpds]
            excessive = [x for x in nonseeds if self.graph.degree(x)==1]


    def make_nodestr(self):
        '''
        create an identifier using nodes in sorted order
        '''
        Nodes = self.graph.nodes()
        Nodes.sort()
        return ''.join(Nodes)

    def make_sif_str(self):
        s = ''
        for e in self.graph.edges(): s += e[0] + ' interact ' + e[1] + '\n'
        return s
        
    def export_sif(self, filename="mummichog_exportx.sif"):
        out = open(filename, 'w')
        out.write(self.make_sif_str())
        out.close()

    def export_network_txt(self, met_model, filename):
        '''
        To use .txt for Cytoscape 3, no need for .sif any more.
        '''
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in self.graph.edges():
            s += e[0] + '\t' + e[1] + '\t' + met_model.edge2enzyme.get(e, '') + '\n'
        
        out = open(filename, 'w')
        out.write(s)
        out.close()



class TableFeatures:
    '''
    Container to look up input mz list and make corresponding cpd list.
    Light weight class for use in permutation.
    The extended AnalysisCentral class organize much of the work.
    
    
    in ver 2,
    moving to single match step upfront
    
    
    '''
    def __init__(self, hsanet, input_mzlist):
        '''
        Hookup to reference network, input parameters.
        self.significant_cpdlist keeps track of cpds that bring significance 
        to pathway/modules.
        
        '''
        self.network = hsanet
        self.input_mzlist = set(input_mzlist)
        self.input_cpdlist = []
        self.make_cpdlist()
        
    def make_cpdlist(self):
        for x in self.input_mzlist:
            # get matched Compound/EmpiricalCompound instances from network.match_dict
            self.input_cpdlist += [C.id for C in 
                        self.network.match_dict.get(x, []) if C.id not in currency]
            
        self.input_cpdlist = set(self.input_cpdlist)

    def count_cpd2mz(self, cpdlist):
        '''
        return the number of m/z features responsible for cpdlist,
        as one m/z could match multiple cpds, skewing the statistics.
        cpdlist from matched only.
        '''
        result = []
        for c in cpdlist: result += [x[0] for x in self.network.cpd_dict[c].hitlist]
        return len(set(result).intersection(self.input_mzlist))



