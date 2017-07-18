# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
Data models in mummichog

Import from serialized Python objects, 
    or a central database in future versions

metabolicModels = {
    human_model_default: [
        ListOfCompounds: [Compound: {ID, name, mw, formula, AdductsAndDerivatives,
                                    }, ...],
        Pathway
        MetabolicNetwork
    
        ],
    
    ...,
}


CpdIndex can be built here, but better precomputed.




>>> import human_model_mfn as model
>>> dir()
['__builtins__', '__doc__', '__name__', '__package__', 'model']
>>> dir(model)
['__builtins__', '__doc__', '__file__', '__name__', '__package__', 
'cpd2pathways', 'cpd_edges', 'dict_cpds_def', 'dict_cpds_mass', 'edge2enzyme', 'edge2rxn', 
'metabolic_pathways', 'metabolic_rxns', 'pathdotdict', 'version']



PROTON = 1.00727646677

def compute_adducts(mw, PROTON):
    aList = [(mw, 'M[1+]'), 
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
                ] + [(mw - PROTON, 'M-H[-]'),
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
    mydict = {}
    for x in aList: mydict[x[1]] = x[0]
    return mydict

# below is a temporary hack; will update
# Main issue is chemical formula, which is not included now, and the adducts should be dependent on chemical formula

>>> from JSON_metabolicModels import metabolicModels
>>> metabolicModels.keys()
['human_model_mfn']
>>> metabolicModels['human_model_mfn'].keys()
['metabolic_rxns', 'metabolic_pathways', 'Compounds', 'cpd2pathways', 'edge2enzyme', 'dict_cpds_mass', 'cpd_edges', 'dict_cpds_def', 'pathdotdict', 'edge2rxn', 'version']
>>> 
>>> metabolicModels['human_model_mfn']['Compounds'].items()[:2]
[('lneldcACP', {'formula': '', 'mw': 0, 'name': 'linoelaidic acid ACP (all trans)', 'adducts': []}), ('G00160', {'formula': '', 'mw': 0, 'name': '(Gal)2 (GalNAc)2 (GlcA)2 (Xyl)1 (Ser)1; Glycoprotein; Glycosaminoglycan', 'adducts': []})]
>>> metabolicModels['human_model_mfn']['Compounds'].items()[92]
('C00217', {'formula': '', 'mw': 147.0532, 'name': 'D-Glutamate; D-Glutamic acid; D-Glutaminic acid; D-2-Aminoglutaric acid',
 'adducts': {'M+2H[2+]': 74.53387646677, 'M+Br81[-]': 227.9695, 'M-H2O+H[1+]': 130.04987646677, 
 'M-C3H4O2+H[1+]': 76.03937646677, 'M-HCOOH+H[1+]': 102.05507646676999, 'M-HCOONa+H[1+]': 80.07307646677, 
 'M+K[1+]': 186.01597646677, 'M+Cl[-]': 182.0221, 'M+Na-2H[-]': 167.02064706646001, 'M-CO2+H[1+]': 104.07067646677, 
 'M+Na[1+]': 170.04247646677, 'M+Br[-]': 225.9715, 'M(S34)-H[-]': 148.04172353323, 'M+H[1+]': 148.06047646677, 
 'M-H4O2+H[1+]': 112.03927646677, 'M(C13)-H[-]': 147.04932353323, 'M(Cl37)-H[-]': 148.04312353323, 'M+HCOONa[1+]': 216.04787646677, 'M(C13)+2H[2+]': 75.03557646677, 'M+HCOOK[1+]': 232.02177646677, 'M-CO+H[1+]': 120.06547646677, 'M+HCOO[-]': 192.050845, 'M(C13)+3H[3+]': 50.359409800103336, 'M(Cl37)+H[1+]': 150.05767646677, 'M-H[-]': 146.04592353323, 'M+ACN-H[-]': 187.07246853323, 'M+Cl37[-]': 184.0191, 'M-H2O-H[-]': 128.03532353322998, 'M(S34)+H[1+]': 150.05627646677002, 'M-HCOOK+H[1+]': 64.09917646677, 'M+3H[3+]': 50.025009800103334, 'M+CH3COO[-]': 206.066495, 'M(C13)+H[1+]': 149.06387646677, 'M[1+]': 147.0532, 'M-NH3+H[1+]': 131.03397646677, 'M+NaCl[1+]': 206.01907646677, 'M+H+Na[2+]': 85.52487646677, 'M+H2O+H[1+]': 166.07107646677002, 'M-H+O[-]': 162.04083353323, 'M+K-2H[-]': 182.99414706646002, 'M-2H[2-]': 72.51932353323001}})
>>> len(metabolicModels['human_model_mfn']['Compounds'])
3560
>>> len(metabolicModels['human_model_mfn']['metabolic_rxns'])
4204


@author: Shuzhao Li
'''


# metabolicModels from .py or from database
# empty adducts should be {}, not [] as in current test version

from JSON_metabolicModels import metabolicModels
import numpy as np
import networkx as nx

# yet testing
class metabolicPathway:
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





class metabolicNetwork:
    '''
    Metabolite-centric metabolic model 
    Theoretical model, no longer containing user data
    
    
    '''
    def __init__(self, MetabolicModel):
        '''
        Initiation of metabolic network model.
        Building Compound index.
        Parsing input files.
        Matching m/z - Compound.
        
        MetabolicModel['Compounds'] are subset of cpds in network/pathways with mw.
        Not all in total_cpd_list has mw.
        
        MetabolicModel needs to have the dicts as below
        
        '''
        #print_and_loginfo( "Loading metabolic network %s..." %MetabolicModel.version ) # version from metabolic model
        
        self.MetabolicModel = MetabolicModel
        self.network = self.build_network(MetabolicModel['cpd_edges'])
        
        self.version = MetabolicModel['version']
        self.Compounds = MetabolicModel['Compounds']
        self.metabolic_pathways = MetabolicModel['metabolic_pathways']
        self.dict_cpds_def = MetabolicModel['dict_cpds_def']
        self.cpd2pathways = MetabolicModel['cpd2pathways']
        self.edge2enzyme = MetabolicModel['edge2enzyme']
        self.total_cpd_list = self.network.nodes()
        
        
    def build_network(self, edges):
        return nx.from_edgelist( edges )
        

    def get_pathways(self):
        pass
 

class Mmodule:
    '''
    Metabolites by their connection in metabolic network.
    A module is a subgraph, while modularity is calculated in 
    the background of reference hsanet.
    
    
    need to record sig EmpCpds
    
    '''
    def __init__(self, network, subgraph, TrioList):
        '''
        TrioList (seeds) format: [(M.row_number, EmpiricalCompounds, Cpd), ...]
        to keep tracking of where the EmpCpd came from (mzFeature).
        
        network is the total parent metabolic network
        '''
        self.network = network
        self.num_ref_edges = self.network.number_of_edges()
        self.num_ref_nodes = self.network.number_of_nodes()
        self.graph = subgraph
        
        seed_cpds = [x[2] for x in TrioList]
        self.shave(seed_cpds)
        self.nodestr = self.make_nodestr()
        self.N_seeds = len(seed_cpds)
        self.A = self.activity_score(seed_cpds, self.get_num_EmpCpd(TrioList))
    
    def activity_score(self, seed_cpds, num_EmpCpd):
        '''
        A * (Ns/Nm)
        A = Newman-Girvan modularity score
        Ns = number of input cpds in module M
        Nm = number of total cpds in M
        Ns/Nm can be corrected as (Ns/total input size)/(Nm/network size), however,
        this normalization factor holds the same in permutations. 
        Use 100 here for network size/total input size.
        
        To reduce bias towards larger modules in Q:
        np.sqrt(len(seed_cpds)/Nm) * 
        
        Ns is now controlled by number of empiricalCompounds
        '''
            
        #Ns = len([x for x in self.graph.nodes() if x in seed_cpds])
        Ns = num_EmpCpd
        Nm = float(self.graph.number_of_nodes())
        if Nm > 0:
            self.compute_modularity()
            return np.sqrt(self.N_seeds/Nm) *self.Q * (Ns/Nm) * 100
        else:
            return 0
        
        
    def get_num_EmpCpd(self, TrioList):
        new = []
        subgraph_nodes = self.graph.nodes()
        for x in TrioList:
            if x[2] in subgraph_nodes:
                new.append(x[1])
                
        return len(set(new))
        
        
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

    def export_network_txt(self, met_model, filename):
        '''
        To use .txt for Cytoscape 3, no need for .sif any more.
        Edges are strings now as switching to JSON compatible.
        '''
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in self.graph.edges():
            s += e[0] + '\t' + e[1] + '\t' + met_model.edge2enzyme.get(','.join(sorted(e)), '') + '\n'
        
        out = open(filename, 'w')
        out.write(s)
        out.close()











