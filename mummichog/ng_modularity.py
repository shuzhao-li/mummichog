"""
ng_modularity.py,
implementation of the module-finding algorithms from

A) Mark E.J. Newman (2006), PNAS 103:8577-8582
B) Rivera, Corban G., Rachit Vakil, and Joel S. Bader. BMC bioinformatics 11.Suppl 1 (2010): S61.

This was from Ca. 2008... 
Future work should improve its efficiency or use a more efficient module-finding algorithm.
Last modified to add fielding numpy.linalg.linalg.linAlgError, Shuzhao Li, 2017-05-15
"""

USE_TEST_MODE = False

from numpy import *

if USE_TEST_MODE:
    import re, itertools
    import networkx as nx
    from scipy.cluster import hierarchy


class module:
    """
    A mudule is a list of nodes within a reference network.
    self.divide(nobj) returns delta_Q and [group1, group2].
    If the division is not successful, returned groups are empty.
    The reference network is passed in via argument nobj.
    """
    def __init__(self, nodes):
        """
        Keep track of nodes by their names (ids)
        because the indexes change from network matrix to subnet.
        node names -> node index via reference_network.crd
        """
        self.nodes = nodes
        self.num_nodes = len(nodes)
        # modularity matrix for this module/subnet, B^(g)
        self.modularity_matrix = []
        # s_i = 1 if node i belongs to group1, -1 if group2
        self.s = [0]*self.num_nodes
        self.delta_Q = 0
        self.to_be_divide = True

    def quick_divide(self, nobj):
        g1, g2 = self.eigen_divide(nobj)
        self.to_be_divide = False
        return self.delta_Q, [g1, g2]
        

    def eigen_divide(self, nobj):
        """
        divide module in two by leading eigenvector.
        This function reproduced the karate club data as given in 
        http://www-personal.umich.edu/~mejn/courses/2004/cscs535/problems5.pdf
        
        Occassionally running into "numpy.linalg.linalg.linAlgError: Eigenvalues did not converge"
        
        """
        self.fetch_modularity_matrix(nobj)
        group1, group2 = [], []
        try:
            [eigenvalues, eigenvectors] = linalg.eig(self.modularity_matrix)
            if eigenvalues.max() > 0:
                # vector corresponding to max eigenvalue, normalized
                lead_vector = eigenvectors[:, eigenvalues.argmax()]
                for ii in range(self.num_nodes):
                    if lead_vector[ii] < 0:
                        group1.append(self.nodes[ii])
                    else:
                        group2.append(self.nodes[ii])
                self.calculate_s(group1, group2)
                self.delta_Q = self.compute_delta_Q(nobj)
        except numpy.linalg.linalg.linAlgError:
            pass
        
        return group1, group2

    def first_Q(self, nobj):
        """
        This is implementation of Eq.[2]. Not really used because 
        it is the same  as self.compute_delta_Q at the first network division.
        """
        s = array(self.s)
        return dot(dot(s.transpose(), nobj.modularity_matrix), s
                                            )/(4.0*nobj.num_edges)

    def compute_delta_Q(self, nobj):
        """
        compute gain in modularity, delta_Q, by Eq.[5].
        In Newman's paper, this step is suggested to be bypassed by Eq.[7].
        However, Eq.[5] can be written as (1/4m) * SUM[(s_i*s_j - 1)*B_(i,j)],
        and using B and B.sum() is more convinient in Numpy.
        I have not evaluated how my approach affects speed of the algorithm,
        though matrix operations in Numpy are respectable.
        """
        dQ = 0
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                dQ += (self.s[ii]*self.s[jj] - 1) * self.modularity_matrix[ii, jj]
        return dQ/(4.0*nobj.num_edges)


    def calculate_s(self, g1, g2):
        for ii in range(self.num_nodes):
            if self.nodes[ii] in g1:
                self.s[ii] = 1
            else:
                self.s[ii] = -1

    def fetch_modularity_matrix(self, nobj):
        am = zeros( (self.num_nodes, self.num_nodes) )
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                am[ii][jj] = nobj.modularity_matrix[nobj.crd[self.nodes[ii]],
                                                    nobj.crd[self.nodes[jj]]]
        self.modularity_matrix = am


class network:
    """
    A network is represented as a list of nodes and a list of edges,
    which can then be translated into an adjacency matrix.
    To stick to numpy consistency, all data are in "array" type.
    Please also note that Python indexing starts from 0.
    To use this class, self.nodes and self.edges should be obtained first,
    either from self.read_gml_network(infile) or elsewhere. 
    Then self.prime_network() fills in degrees and matrices etc.
    Finally, self.module_analyze() finds all modules in this network.
    """
    def __init__(self):
        """
        Network attributes, predefined for code clarity.
        Matrix in 'array' type (see class docstr above).
        node = vertex.
        An edge from a to b is written as (a, b).
        """
        self.desc = ''
        self.adjacency_matrix = []
        self.modularity_matrix = []
        self.nodes = []
        self.edges = []
        self.num_nodes = 0
        self.num_edges = 0
        self.degrees = []   # same sequence as in self.nodes
        self.modules = []   # module instances
        self.Q = 0          # modularity value, 0 for undivided network
        self.crd = {}       # coordinate for nodes


    def specsplit(self):
        """
        produce modules by eigenvector method alone
        """
        self.prime_network()
        m = module(self.nodes)
        self.modules = [m]
        divisible = [x for x in self.modules if x.to_be_divide]
        while divisible:
            for mod in divisible:
                #print "working on module with", mod.num_nodes
                delta_Q, newgroups = mod.quick_divide(self)
                # check if division successful
                if delta_Q > 0 and newgroups[0]:
                    self.Q += delta_Q
                    self.modules.remove(mod)
                    for g in newgroups:
                        self.modules.append(module(g))
                    
                # else mod is final
            divisible = [x for x in self.modules if x.to_be_divide]
        
        return [m.nodes for m in self.modules]


    #
    # house-keeping functions
    #

    def prime_network(self):
        """
        prepare network for analysis, starting from nodes and edges
        """
        if self.nodes and self.edges:
            self.num_edges = len(self.edges)
            self.num_nodes = len(self.nodes)
            self.make_node_index()
            self.make_adjacency_matrix()
            self.compute_degrees()
            self.make_modularity_matrix()
        else:
            raise Exception, "Network nodes or edges not properly defined!"

    def make_modularity_matrix(self):
        """
        Modularity matrix is defined in Eq.[3],
        B_(i,j) = A_(i,j) - (k_i*k_j)/2m.
        """
        am = zeros( (self.num_nodes, self.num_nodes) )
        for ii in range(self.num_nodes):
            for jj in range(self.num_nodes):
                am[ii, jj] = self.adjacency_matrix[ii, jj] - \
                         self.degrees[ii]*self.degrees[jj]/(2.0*self.num_edges)
        self.modularity_matrix = am
        
    def compute_degrees(self):
        """
        Compute degree for each node from adjacency matrix
        """
        if not self.degrees:
            for ii in range(self.num_nodes):
                self.degrees.append(self.adjacency_matrix[ii].sum())
        else:
            print "Degrees already exist!"

    def make_adjacency_matrix(self):
        """
        Make adjacency matrix from nodes and edges.
        Edges are treated as not directional, thus matrix is symmetrical.
        Multiple edges between two nodes are allowed.
        """
        am = zeros( (self.num_nodes, self.num_nodes) )
        for edge in self.edges:
            am[self.nodes.index(edge[0]), self.nodes.index(edge[1])] += 1
            # take out the next line will make a directional network
            am[self.nodes.index(edge[1]), self.nodes.index(edge[0])] += 1
        self.adjacency_matrix = am

    def make_node_index(self):
        """
        mapping node name/id to original index
        """
        for ii in range(self.num_nodes):
            self.crd[self.nodes[ii]] = ii

    #
    # The read_gml functions read network from GML format.
    # Alternatively, self.nodes and self.edges can be directly defined.
    #

    def copy_from_graph(self, g):
        '''
        g as a graph instance from nx.graph
        '''
        self.nodes, self.edges = g.nodes(), g.edges()

    def read_gml_network(self, infile):
        (self.nodes, self.edges) = self.read_gml(infile)
        
    def read_sif_network(self, infile):
        (self.nodes, self.edges) = self.read_sif(infile)

    def read_gml(self, infile):
        """
        read a GML file 
        (http://www.infosun.fim.uni-passau.de/Graphlet/GML/gml-tr.html),
        such as the karate.gml used in this demo,
        and return lists of nodes and edges. Based on regular expressions.
        This is not a bullet-proof parser of GML. But if you have looked at 
        your file and will check the result, this function is likely to get 
        it right. Only ids are used; labels are ignored.
        An alternative Python GML parser can be found 
        in NetworkX (networkx.lanl.gov).
        """
        def scan_node(str_block):
            # input example: 'node   [     id 5   ]'
            a = str_block.split()
            #return the id after 'id'
            return a[a.index('id') + 1]
        def scan_edge(str_block):
            # input example: 'edge   [     source 4     target 2   ]'
            a = str_block.split()
            sid, tid = a[a.index('source') + 1], a[a.index('target') + 1]
            return (sid, tid)
        s = open(infile).read()
        # remove line breaks to ease re.findall
        s = " ".join(s.splitlines())
        nodes = re.findall('node.*?\\[.*?id.*?\\]', s)
        edges = re.findall('edge.*?\\[.*?source.*?\\]', s)
        print "Found %d nodes and %d edges." %(len(nodes), len(edges))
        nodes = [scan_node(x) for x in nodes]
        edges = [scan_edge(x) for x in edges]
        return (nodes, edges)

    def read_sif(self, infile):
        """
        read sif file as defined in CytoScape; 
        simple hack that doesn't comply all sif definitinos.
        """
        w = open(infile).readlines()
        nodes, edges = [], []
        for line in w:
            a = line.strip().split()
            nodes.append(a[0])
            nodes.append(a[-1])
            edges.append((a[0], a[-1]))
        edges = list( set(edges) )
        nodes = list( set(nodes) )
        return (nodes, edges)
            

#
# NeMo algorithm for comparison
#

class nemo_network:
    """
    Implementation of 
    Rivera, Corban G., Rachit Vakil, and Joel S. Bader. "NeMo: network module identification in Cytoscape." 
    BMC bioinformatics 11.Suppl 1 (2010): S61.
    This class depends on Networkx.
    """
    def __init__(self, network):
        """
        network as nx.graph
        """
        self.network = network
        self.e = float( network.number_of_edges() )
        self.e_total = self.compute_e_total()
        self.sorted_pairs = []
        self.good_node_list = []        # nodes with r > 0
        self.distance_matrix = []
        self.modules = []

    def compute_e_total(self):
        '''
        part of lamda_expected equation
        '''
        et = 0
        for n in self.network.nodes():
            degree = self.network.degree(n)
            et += degree * (degree - 1)
            
        return et

    def compute_lamda_expected(self, n_a, n_b):
        """
        lamda_expected according to Equation (7)
        """
        return 0.25 * (self.e_total - (n_a*(n_a - 1) + n_b*(n_b - 1))) * n_a * n_b / (self.e)**2
        
    def compute_r(self, s_ab, lamda_expected):
        if s_ab * lamda_expected == 0:
            return 0
        else:
            return s_ab * log(s_ab / lamda_expected) - abs(s_ab - lamda_expected)
        
    def find_modules(self):
        '''
        return a list of modules, which are lists of nodes, 
        to be consistent with network.specsplit output
        '''
        self.compute_r_list()
        self.make_distance_matrix()
        self.hcl()
        return self.modules
    
    def compute_r_list(self):
        '''
        store list of node pairs
        by r in descending order
        '''
        mydict = {}
        for a, b in itertools.combinations( self.network.nodes(), 2):
            n_a = self.network.degree(a)
            n_b = self.network.degree(b)
            s_ab = len( set(self.network.neighbors(a)).intersection(set(self.network.neighbors(b))) )
            lamda_expected = self.compute_lamda_expected(n_a, n_b)
            mydict[ (a, b) ] = (self.compute_r(s_ab, lamda_expected), n_a, n_b, s_ab, lamda_expected)
            
        new = [(v[0], k) for k,v in mydict.iteritems() if v[0] > 0]
        new.sort(reverse=True)
        # [((a, b), r), ...]
        self.sorted_pairs = [ (x[1], x[0]) for x in new ]
        
        
    def make_distance_matrix(self):
        nodes = []
        for x in self.sorted_pairs: nodes += list(x[0])
        nodes = list(set(nodes))
        self.good_node_list = nodes
        self.N = N = len(nodes)
        nodedict = {}
        for ii in range(N): nodedict[ nodes[ii] ] = ii
        
        m = zeros( (N, N) )
        for x in self.sorted_pairs:
            m[nodedict[x[0][0]], nodedict[x[0][1]]] = x[1]
            m[nodedict[x[0][1]], nodedict[x[0][0]]] = x[1]
        
        self.distance_matrix = m
        
        
    def hcl(self):
        '''
        agglomerative clustering tree,
        represented as a linkage matrix.
        Then parsed by Adrian Rosebrock's method.
        Remove repeats.
        '''
        clusters = hierarchy.linkage(self.distance_matrix, "complete", "euclidean")
        new = []
        for num in range(2, self.N/2):
            new += self.parseClusters(clusters, num).values()
        
        new = [list(x) for x in new if len(x) > 2]
        for x in new: x.sort()
        
        modules = []
        for x in new:
            modules.append( tuple([self.good_node_list[ii] for ii in x]) )
            
        modules = list(set(modules))
        self.modules = [list(x) for x in modules]
        
    
    def parseClusters(self, clusters, numClusters):
        '''
        Courtesy of Adrian Rosebrock
        '''
        # initialize the parsed clusters dictionary and the current
        # cluster index
        parsed = {}
        #j = len(self.index)
        j = self.N

        # loop over the clusters
        for i in range(0, self.N - numClusters):
            # grab the index of the two clusters that were merged
            # together
            idA = int(clusters[i][0])
            idB = int(clusters[i][1])
    
            # try to grab the set of nodes already in the cluster
            # for each ID
            nodesA = parsed.get(idA, set())
            nodesB = parsed.get(idB, set())
    
            # if the first node set is larger than zero, delete the
            # cluster ID
            if len(nodesA) > 0:
                del parsed[idA]
    
            # if the second node set is larger than zero, delete the
            # cluster ID
            if len(nodesB) > 0:
                del parsed[idB]

            # if the first ID is less than the length of the index, then
            # add the ID to the nodes
            if idA < self.N:
                nodesA.add(idA)
    
            # if the second ID is less than the length of the index, then
            # add the ID to the nodes
            if idB < self.N:
                nodesB.add(idB)
    
            # take the union of the set of nodes, update the parsed clusters
            # dictionary, and then increment the next available cluster ID
            nodesA = nodesA.union(nodesB)
            parsed[j] = nodesA
            j += 1

        # return the parsed clusters dictionary
        return parsed




