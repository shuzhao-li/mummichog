# Copyright (c) 2010-2017 Shuzhao Li, Andrei Todor
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
Pathway, module analysis in mummichog;
then compute activity network.
Output includes HTML report, result.html, metabolite data and visualization files for Cytoscape 3.

Major change from version 1 to version 2: using EmpiricalCompound in place of cpd.

Separating I/O, to be used for both web apps and desktop apps


@author: Shuzhao Li, Andrei Todor
'''

import logging, random, itertools
from scipy import stats
import ng_modularity as NGM

from get_user_data import *


logging.basicConfig(format='%(message)s', level=logging.INFO)


# --------------------------------------------------------
#
# pathway analysis
#

class PathwayAnalysis:
    '''
    From matched features to pathway enrichment analysis.
    Using mfn human pathways for now.
    p-value is from Fisher exact test, 
    adjusted by resampling method in 
    GF Berriz, OD King, B Bryant, C Sander & FP Roth. 
    Characterizing gene sets with FuncAssociate. 
    Bioinformatics 19(18):2502-2504 (2003)
    
    "Adjusted_p" is not an accurate term. It's rather an empirical p-value.
    
    Note pathway_size is not different from version 1.
    
    version 2 moved everything into EmpiricalCompound space.
    
    
    
    '''
    def __init__(self, pathways, mixedNetwork):
        '''
        mixedNetwork contains both user input data, metabolic model,
        and mapping btw (mzFeature, EmpiricalCompound, cpd)
        
        '''
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        self.paradict = mixedNetwork.data.paradict
        
        self.pathways = self.get_pathways(pathways)
        self.resultListOfPathways = []          # will store result of pathway analysis
        
        # to help track wehre sig cpd comes from
        self.TrioList = self.mixedNetwork.TrioList
        self.significant_EmpiricalCompounds = set([x[1] for x in self.TrioList])
        
        self.ListOfEmpiricalCompounds = mixedNetwork.ListOfEmpiricalCompounds
        self.total_number_EmpiricalCompounds = len(self.ListOfEmpiricalCompounds)

        print_and_loginfo("\nPathway Analysis...")
        
        
    def get_pathways(self, pathways):
        '''
        convert pathways in JSON formats (import from .py) to list of Pathway class.
        Adding list of EmpiricalCompounds per pathway, which reflects the measured pathway coverage.
        '''
        new = []
        for j in pathways:
            P = metabolicPathway()
            P.json_import(j)
            P.adjusted_p = ''
            P.EmpiricalCompounds = self.__get_empiricalCompounds_by_cpds__(P.cpds)
            new.append(P)
        return new
        

    def __get_empiricalCompounds_by_cpds__(self, cpds):
        '''
        Mapping cpds to empirical_cpds. Also used for counting EmpCpds for each Pathway.
        '''
        cpds_empirical = []
        for c in cpds: cpds_empirical += self.mixedNetwork.Compounds_to_EmpiricalCompounds.get(c, [])
        return set(cpds_empirical)
        
        
    def do_permutations(self, pathways, num_perm):
        '''
        Modified from Berriz et al 2003 method.
        After collecting p-values from resampling, do a Gamma fit.
        
        Permutation is simplified in version 2; no more new TableFeatures instances.
        
        
        May consider fitting Gamma at log scale, to be more accurate --
        
        '''
        self.permutation_record = []
        print_and_loginfo("Resampling, %d permutations to estimate background ..." 
                          %num_perm)
        
        # this is feature number not cpd number
        N = len(self.mixedNetwork.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii + 1))
            sys.stdout.flush()
            random_Trios = self.mixedNetwork.batch_rowindex_EmpCpd_Cpd( random.sample(self.mixedNetwork.mzrows, N) )
            query_EmpiricalCompounds = set([x[1] for x in random_Trios])
            self.permutation_record += (self.__calculate_p_ermutation_value__(query_EmpiricalCompounds, pathways))
        
        print_and_loginfo("\nPathway background is estimated on %d random pathway values" 
                          %len(self.permutation_record))
        


    def __calculate_p_ermutation_value__(self, query_EmpiricalCompounds, pathways):
        '''
        calculate the FET p-value for all pathways.
        But not save anything to Pathway instances.
        '''
        p_of_pathways = [ ]
        query_set_size = len(query_EmpiricalCompounds)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        for P in pathways:
            overlap_features = query_EmpiricalCompounds.intersection(P.EmpiricalCompounds)
            overlap_size = len(overlap_features)
            ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                p_val = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                       [ecpd_num - overlap_size, negneg]], 'greater')[1]
                p_of_pathways.append(p_val)
            else: 
                p_of_pathways.append(1)
                
        return p_of_pathways


    def get_adjust_p_by_permutations(self, pathways):
        '''
        EASE score is used as a basis for adjusted p-values,
        as mummichog encourages bias towards more hits/pathway.
        pathways were already updated by first round of Fisher exact test,
        to avoid redundant calculations
        '''
        self.do_permutations(pathways, self.paradict['permutation'])
        
        if self.paradict['modeling'] == 'gamma':
            #vector_to_fit = [-np.log10(x) for x in self.permutation_record if x < 1]
            vector_to_fit = -np.log10(np.array(self.permutation_record))
            self.gamma = stats.gamma.fit(vector_to_fit)
            a, loc, scale = self.gamma
            
            for P in pathways: 
                P.adjusted_p = self.__calculate_gamma_p__(a, loc, scale, P.p_EASE)
        else:
            for P in pathways: P.adjusted_p = self.__calculate_p__(P.p_EASE, self.permutation_record)
        return pathways
        

    def __calculate_p__(self, x, record):
        '''
        calculate p-value based on the rank in record of permutation p-values
        '''
        total_scores = [x] + record
        total_scores.sort()
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
    
    def __calculate_gamma_p__(self, a, loc, scale, x):
        '''
        Use -log10 scale for model fitting
        '''
        return 1 - stats.gamma.cdf(-np.log10(x), a, loc, scale)
    
    
    def cpd_enrich_test(self):
        '''
        Fisher Exact Test in cpd space, after correction of detected cpds.
        Fisher exact test is using scipy.stats.fisher_exact
        for right-tail p-value:
        >>> stats.fisher_exact([[12, 5], [29, 2]], 'greater')[1]
        0.99452520602188932
        
        query size is now counted by EmpiricalCompounds.
        adjusted_p should be model p-value, not fdr.
        This returns a list of Pathway instances, with p-values.
        
                        P.p_EASE = stats.fisher_exact([[max(0, overlap_size - 1), query_set_size - overlap_size],
                                   [ecpd_num - overlap_size + 1, negneg]], 'greater')[1]
        '''
        FET_tested_pathways = []
        qset = self.significant_EmpiricalCompounds
        query_set_size = len(qset)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        print_and_loginfo("Query number of significant compounds = %d compounds" %query_set_size)
        
        for P in self.pathways:
            # use the measured pathway size
            P.overlap_EmpiricalCompounds = P.overlap_features = qset.intersection(P.EmpiricalCompounds)
            P.overlap_size = overlap_size = len(P.overlap_EmpiricalCompounds)
            P.EmpSize = ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                # Fisher's exact test
                P.p_FET = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                   [ecpd_num - overlap_size, negneg]], 'greater')[1]
                # EASE score as in Hosack et al 2003
                # taking out EASE, as the new approach of EmpiricalCompound is quite stringent already
                P.p_EASE = P.p_FET
                

            else:
                P.p_FET = P.p_EASE = 1
                
            FET_tested_pathways.append(P)
            #  (enrich_pvalue, overlap_size, overlap_features, P) 
            
        result = [(P.adjusted_p, P) for P in 
                                        self.get_adjust_p_by_permutations(FET_tested_pathways)]
        result.sort()
        self.resultListOfPathways = [x[1] for x in result]

    
    def collect_hit_Trios(self):
        '''
        get [(mzFeature, EmpiricalCompound, cpd),...] for sig pathways.
        Nominate top cpd for EmpCpd here, i.e.
        in an EmpCpd promoted by a significant massFeature, the cpd candidate is chosen from a significant pathway.
        If more than one cpds are chosen, keep multiple.
        
        '''
        overlap_EmpiricalCompounds = set([])
        for P in self.resultListOfPathways:
            if P.adjusted_p < SIGNIFICANCE_CUTOFF:
                # print(P.adjusted_p, P.name)
                overlap_EmpiricalCompounds = overlap_EmpiricalCompounds.union(P.overlap_EmpiricalCompounds)
        
        new = []
        for T in self.TrioList:
            # [(mzFeature, EmpiricalCompound, cpd),...]
            if T[1] in overlap_EmpiricalCompounds and T[0] in self.mixedNetwork.significant_features:
                # this does not apply to all sig EmpCpd
                T[1].update_chosen_cpds(T[2])
                T[1].designate_face_cpd()
                new.append(T)
        
        return new
                    
    
    def plot_model_pvalues(self, outfile='mcg_pathway_modeling'):
        '''
        Plot self.permutation_record
        P.p_EASE for P in self.resultListOfPathways
        Use -log10 scale, to show upward trend, consistent with other plots
        '''
        self.permutation_record.sort()
        Y_data = [-np.log10(x) for x in self.permutation_record]
        fig = plt.figure(figsize=(5,4))
        plt.plot(range(len(Y_data)), Y_data, 'b.')
        for P in self.resultListOfPathways[:10]:
            YY = -np.log10(P.p_EASE)
            plt.plot([0, 0.1*len(Y_data)], [YY, YY], 'r-')
        
        plt.ylabel("-log10 (FET p-value)")
        plt.xlabel("Number of permutation")
        plt.title("Modeling pathway significance")
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')

    
    def plot_bars_top_pathways(self, outfile='mcg_pathway_barplot'):
        '''
        Horizontal barplot of pathways.
        Also returnin-memory string for web use
        '''
        use_pathways = [P for P in self.resultListOfPathways if P.adjusted_p < SIGNIFICANCE_CUTOFF]
        if len(use_pathways) < 6:
            use_pathways = self.resultListOfPathways[:6]
        #plot use_pathways
        fig, ax = plt.subplots()
        ylabels = [P.name for P in use_pathways]
        data = [-np.log10(P.adjusted_p) for P in use_pathways]
        NN = len(data)
        ax.barh( range(NN), data, height=0.5, align='center', color="purple", alpha=0.4 )
        ax.set_yticks(range(NN))
        ax.set_yticklabels(ylabels)
        ax.set_xlabel('-log10 p-value')
        
        ax.plot([1.301, 1.301], [-0.5, NN], 'g--')    # NN is inverted too
        #ax.set_ylim(-0.5, NN)
        
        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')
        
        # get in-memory string for web use
        figdata = StringIO.StringIO()
        plt.savefig(figdata, format='png')
        figdata.seek(0)
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(figdata.buf))
        return '<img src = "%s"/>' % uri
        
    
    

# --------------------------------------------------------
#
# module analysis
#

class ModularAnalysis:
    '''
    1) Find modules from input list by connecting paths < 4;
    compute activity score that combines modularity and enrichment.
    2) Permutations by randomly selecting features from ref_mzlist;
    compute p-values based on permutations.
    
    Working on version 2:
    Module analysis will still be in the compound space, as network model is defined by compound edges.
    
    
    Need tracking the mapping btw compound and EmpiricalCompounds
    
    
    Tested to add a generator from EmpiricalCompounds to a bunch of combinations, 
    i.e. only one cpd from Ecpd is used at a time towards module analysis
    But it's too slow to be practical.
    
    
    '''
    def __init__(self, mixedNetwork):
        '''
        mapping btw (mzfeature, cpd) has to be via ListOfEmpiricalCompounds, 
        so that cpd can be tracked back to EmpiricalCompounds
        
        
        '''
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        self.paradict = mixedNetwork.data.paradict
        
        # both using row_numbers
        self.ref_featurelist = self.mixedNetwork.mzrows
        self.significant_features = self.mixedNetwork.significant_features
        self.significant_Trios = self.mixedNetwork.TrioList
        

    def dispatch(self):
        '''
        Only real modules are saved in total. 
        Permutated modules are not saved but their scores are recorded. 
        
        
        print("what found? -", self.modules_from_significant_features[3].nodestr)
        
        
        '''
        s = "\nModular Analysis, using %d permutations ..." %self.paradict['permutation']
        #print s
        logging.info(s)
        self.modules_from_significant_features = self.run_analysis_real()
        self.permuation_mscores = self.do_permutations(self.paradict['permutation'])

        self.rank_significance()        
        #for M in self.top_modules: print(M, M.A, nx.average_node_connectivity(M.graph))


    def run_analysis_real(self):
        return self.find_modules( self.significant_Trios )
        

    def do_permutations(self, num_perm):
        '''
        Run num_perm permutations on ref featurelist;
        populate activity scores from random modules in self.permuation_mscores
        
        '''
        permuation_mscores = []
        N = len(self.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii+1))
            sys.stdout.flush()
            random_trios = self.mixedNetwork.batch_rowindex_EmpCpd_Cpd( 
                                            random.sample(self.ref_featurelist, N) )
            permuation_mscores += [x.A for x in self.find_modules(random_trios)] or [0]
            
        return permuation_mscores
            
    def __generator_EmpiricalCompounds_cpds__(self, Ecpds):
        '''
        return the combinations of one cpd drawn from each EmpiricalCompound, N = len(Ecpds)
        as iterator
        This function is not used now, because too many combinations make software too slow.
        '''
        return itertools.product(*[ E.compounds for E in Ecpds ])


    def find_modules(self, TrioList):
        '''
        get connected nodes in up to 4 steps.
        modules are set of connected subgraphs plus split moduels within.
        A shaving step is applied to remove excessive nodes that do not 
        connect seeds (thus Mmodule initiation may reduce graph size). 
        A module is only counted if it contains more than one seeds.
        
        TrioList format: [(M.row_number, EmpiricalCompounds, Cpd), ...]
        '''
        global SEARCH_STEPS, MODULE_SIZE_LIMIT
        seeds = [x[2] for x in TrioList]
        modules, modules2, module_nodes_list = [], [], []

        for ii in range(SEARCH_STEPS):
            edges = nx.edges(self.network, seeds)
            if ii == 0:
                # step 0, counting edges connecting seeds
                edges = [x for x in edges if x[0] in seeds and x[1] in seeds]
                new_network = nx.from_edgelist(edges)
                
            else:
                # step 1, 2, 3, ... growing to include extra steps/connections
                new_network = nx.from_edgelist(edges)
                seeds = new_network.nodes()
            
            for sub in nx.connected_component_subgraphs(new_network):
                if 3 < sub.number_of_nodes() < MODULE_SIZE_LIMIT:
                    M = Mmodule(self.network, sub, TrioList)
                    modules.append(M)
                
        # add modules split from modules
        if USE_DEBUG:
            logging.info( '# initialized network size = %d' %len(seeds) )
            # need export modules for comparison to heinz
            self.__export_debug_modules__( modules )
            
        for sub in modules:
            if sub.graph.number_of_nodes() > 5:
                modules2 += [Mmodule(self.network, x, TrioList)
                             for x in self.__split_modules__(sub.graph)]
        
        new = []
        for M in modules + modules2:
            if M.graph.number_of_nodes() > 3 and M.nodestr not in module_nodes_list:
                new.append(M)
                module_nodes_list.append(M.nodestr)
                if USE_DEBUG: logging.info( str(M.graph.number_of_nodes()) + ', ' + str(M.A) )

        return new


    def __export_debug_modules__(self, modules):
        '''
        write out initial modules, to be split by alternative algorithm
        '''
        s = ''
        for M in modules: s += M.make_sif_str()
        out = open( os.path.join(self.modules_dir, 'debug_modules.sif'), 'a' )
        out.write(s + '#\n')
        out.close()
        
    def __split_modules__(self, g):
        '''
        return nx.graph instance after splitting the input graph
        by Newman's spectral split method
        Only modules more than 3 nodes are considered as good small modules 
        should have been generated in 1st connecting step.
        '''
        net = NGM.network()
        net.copy_from_graph(g)
        return [nx.subgraph(g, x) for x in net.specsplit() if len(x) > 3]

    # test alternative algorithm
    def __split_modules_nemo__(self, g):
        '''
        Alternative function using NeMo algorithm for module finding.
        Not used for now.
        '''
        net = NGM.nemo_network(g)
        return [nx.subgraph(g, x) for x in net.find_modules() if len(x) > 3]


    def rank_significance(self):
        '''
        compute p-values of modules. Either model based:
        scores of random modules are fitted to a Gamma distribution,
        p-value is calculated from CDF.
        Or rank based.
        '''
        print_and_loginfo("\nNull distribution is estimated on %d random modules" 
                          %len(self.permuation_mscores))
        print_and_loginfo("User data yield %d network modules" 
                          %len(self.modules_from_significant_features))
        
        if self.paradict['modeling'] == 'gamma':
            a, loc, scale = stats.gamma.fit(self.permuation_mscores)
            if USE_DEBUG:
                logging.info( 'Gamma fit parameters a, loc, scale = ' + ', '.join([str(x) for x in [a, loc, scale]]) )
            
            for M in self.modules_from_significant_features:
                M.p_value = 1 - stats.gamma.cdf(M.A, a, loc, scale)
            
        else:
            for M in self.modules_from_significant_features:
                M.p_value = self.__calculate_p__(M.A, self.permuation_mscores)
                
        top_modules = [M for M in self.modules_from_significant_features if M.p_value < SIGNIFICANCE_CUTOFF]
        self.top_modules = sorted(top_modules, key=lambda M: M.p_value)
        
        
        


    def __calculate_p__(self, x, record):
        '''
        calculate p-value based on the rank in record of scores
        '''
        total_scores = [x] + record
        total_scores.sort(reverse=True)
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
  
    def collect_hit_Trios(self):
        '''
        get [(mzFeature, EmpiricalCompound, cpd),...] for top_modules.
        Update EmpCpd chosen compounds.
        '''
        overlap_Cpds = []
        for M in self.top_modules:
            overlap_Cpds += M.graph.nodes()
        
        overlap_Cpds = set(overlap_Cpds)
        new = []
        for T in self.significant_Trios:
            if T[2] in overlap_Cpds:
                T[1].update_chosen_cpds(T[2])
                T[1].designate_face_cpd()
                new.append(T)
        
        return new



    def plot_model_pvalues(self, outfile='mcg_module_modeling.pdf'):
        '''
        Plot module activity against self.permuation_mscores
        
        '''
        self.permuation_mscores.sort(reverse=True)
        NN = len(self.permuation_mscores)
        fig = plt.figure(figsize=(5,4))
        plt.plot(range(NN), self.permuation_mscores, 'bo')
        for M in self.modules_from_significant_features:
            plt.plot([0, 0.1*NN], [M.A, M.A], 'r-')
        
        plt.ylabel("Activity score")
        plt.xlabel("Number of permutation")
        plt.title("Modeling module significance")
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')
        
    
    def draw_top_modules(self, outfile_prefix='mcg_module_.pdf'):
        for M in self.top_modules:
            #draw it
            #
            pass
    
    
# --------------------------------------------------------
#
# activity network analysis
#

class ActivityNetwork:
    '''
    Tally cpds responsible for significant pathways and modules,
    and build an acitvity network to represent top story in the data.
    Remove singletons in network. Try no more than 3 steps. AN can get too big or too small.
    
    '''
    def __init__(self, mixedNetwork, hit_Trios):
        '''
        Build a consensus network for hit_Trios,
        [(mzFeature, EmpiricalCompound, cpd),...] for top_modules and sig pathways.
        hit_Trios = set(PA.collect_hit_Trios() + MA.collect_hit_Trios())
        '''
        # also update to mixedNetwork
        mixedNetwork.hit_Trios = hit_Trios
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        
        nodes = [x[2] for x in hit_Trios]
        
        self.activity_network = self.build_activity_network(nodes)
        

    def build_activity_network(self, nodes, cutoff_ave_conn = 0.5, expected_size = 10):
        '''
        Get a network with good connections in no more than 3 steps.
        No modularity requirement for 1 step connected nodes.
        '''
        
        an = nx.subgraph(self.mixedNetwork.model.network, nodes)
        sub1 = self.__get_largest_subgraph__(an)
        
        if sub1.number_of_nodes() > expected_size:
            print_and_loginfo("\nActivity network was connected in 1 step.")
            return sub1
        
        else:   # expand 1 or 2 steps
            edges = nx.edges(self.mixedNetwork.model.network, nodes)
            new_network = self.__get_largest_subgraph__( nx.from_edgelist(edges) )
            conn = self.__get_ave_connections__(new_network)
            if an.number_of_nodes() > MODULE_SIZE_LIMIT or conn > cutoff_ave_conn:
                print_and_loginfo("\nActivity network was connected in 2 steps.")
                return new_network
            else:
                edges = nx.edges(self.mixedNetwork.model.network, new_network.nodes())
                new_network = self.__get_largest_subgraph__( nx.from_edgelist(edges) )
                conn = self.__get_ave_connections__(new_network)
                if conn > cutoff_ave_conn:
                    print_and_loginfo("\nActivity network was connected in 3 steps.")
                    return new_network
                else:
                    return an
    
    def export_network_txt(self, met_model, filename):
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in self.activity_network.edges():
            s += e[0] + '\t' + e[1] + '\t' + met_model.edge2enzyme.get(','.join(sorted(e)), '') + '\n'
        
        out = open(filename, 'w')
        out.write(s)
        out.close()
        
    def __get_largest_subgraph__(self, an):
        '''
        connected_component_subgraphs likely to return sorted subgraphs. Just to be sure here.
        '''
        connected = [(len(x),x) for x in nx.connected_component_subgraphs(an)]
        connected.sort(reverse=True)
        return connected[0][1]
        
    def __get_ave_connections__(self, N):
        '''
        nx.average_node_connectivity(G) is too slow; use self.__get_ave_connections__()
        '''
        try:                #Avoid ZeroDivisionError
            return N.number_of_edges()/float(N.number_of_nodes())
        except:
            return 0