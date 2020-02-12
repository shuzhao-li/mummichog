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

Import from serialized Python objects.
(A server version will use or a central database)

metabolicModels = {
    human_model_default: [
        ListOfCompounds: [Compound: {ID, name, mw, formula, AdductsAndDerivatives,
                                    }, ...],
        Pathway
        MetabolicNetwork
    
        ],
    
    ...,
}

@author: Shuzhao Li
'''
from .JSON_metabolicModels import metabolicModels
import numpy as np
import networkx as nx


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
    Theoretical model, not containing user data
    '''
    def __init__(self, MetabolicModel):
        '''
        Initiation of metabolic network model.
        Building Compound index.
        Parsing input files.
        Matching m/z - Compound.
        
        MetabolicModel['Compounds'] are subset of cpds in network/pathways with mw.
        Not all in total_cpd_list has mw.
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
