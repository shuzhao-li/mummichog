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
Reporting functions in mummichog;
to write output results and do visualization.

@author: Shuzhao Li

'''


import csv, xlsxwriter, logging
import numpy as np
from websnippets import *

try:
    import pygraphviz as pgv
except ImportError:
    print("Pygraphviz is not found. Skipping...")



class HtmlExport:
    def __init__(self):
        self.elements = []
        self.jsdata = ''
        
        self.HTML_HEAD = HTML_HEAD
        self.HTML_END = HTML_END
        self.javascript_HEAD = javascript_HEAD
        self.javascript_END = javascript_END
        
        
    def write_tag(self, s, tag, classname='', htmlid=''):
        cn = ''
        if htmlid: cn += ' id="' + htmlid + '"'
        if classname: cn += ' class="' + classname + '"'
        
        return '<' + tag + cn + '>' + s + '</' + tag + '>'


    def add_element(self, s, tag, classname='', htmlid=''):
        self.elements.append( self.write_tag(s, tag, classname, htmlid) + '\n')

    def make_select_menu(self, N):
        '''
        Showing
        <select id="nmodule_select" onchange="updateDrawing();">
          <option value="an">Activity network</option>
          <option value="module1">Module 1</option>
          <option value="module2">Module 2</option>
        </select>
        Node size
        <select id="nsize_select"  onchange="updateDrawing();">
          <option value=16>Medium</option>
          <option value=8>Small</option>
          <option value=32>Large</option>
        </select>
        '''
        s = '''Showing
        <select id="nmodule_select" onchange="updateDrawing();">
        '''
        for m in ["Module " + str(ii) for ii in range(1, N)]:
            s += "<option>%s</option>" %m
            
        return s + '''  <option>Activity network</option>
                   </select>
                    Node size <select id="nsize_select"  onchange="updateDrawing();">
                      <option value=16>Medium</option>
                      <option value=8>Small</option>
                      <option value=32>Large</option>
                    </select>
                    Visual style <select id="vstyle_select"  onchange="updateDrawing();">
                      <option value=0>Drag</option>
                      <option value=1>Force</option>
                    </select>
                    '''


    def make_js_data(self, networks, cpdnamedict, dict_cpd_foldchange):
        '''
        Make data for both d3.js and cytoscape.js visualization.
        
        d3js visualization from multiple networks.
            var nodes = [
                [ {name:"Myriel",group:1},
                {name:"Napoleon",group:1},
                {name:"Mlle.Baptistine",group:1},
                {name:"Mme.Magloire",group:1},
                {name:"CountessdeLo",group:1},], ...];
                
                var links = [
               [ {source:1,target:0},
                {source:2,target:0},
                {source:3,target:0,value:10},
                {source:5,target:0,value:1},], ...];
        
        cytoscape.js visualization:
        	var cytonodes = [ [
				{ data: { id: "C01227", weight: 2.11 } },
				{ data: { id: "C05475", weight: 3.65 } },
				{ data: { id: "C05487", weight: 3.65 } },
				{ data: { id: "C03205", weight: -1.8 } },
					  ], ... ]
	  
					var cytoedges = [ [
				{ data: { id: "C01227C05498", weight: 1, source: "C01227", target: "C05498" } },
				{ data: { id: "C01227C05487", weight: 1, source: "C01227", target: "C05487" } },
				{ data: { id: "C01227C03205", weight: 1, source: "C01227", target: "C03205" } }, ... ]]
        
        Color coded by group. Thru additional argument colordict.
        The visualization can get a lot more complicated. Using only names for cytoscape.js for now. 
        When implementing other features, both id and name will be needed.
        '''
        cpdcolordict = self.rescale_color(dict_cpd_foldchange)
        
        nodestr, edgestr, cynodestr, cyedgestr = '', '', '', ''
        for network in networks:
            nodestr += '['
            cynodestr += '['
            ordered_nodes = network.nodes()
            nodesdict, ii = {}, -1
            for n in ordered_nodes:
                ii += 1
                nodesdict[n] = ii
                nodestr += '{name:%s, group:%d}, ' %(quote(cpdnamedict.get(n, n).split(';')[0]), cpdcolordict.get(n, 10))
                cynodestr += '{ data: ' + '{id:%s, group:%d} }, ' %(quote(cpdnamedict.get(n, n).split(';')[0]), cpdcolordict.get(n, 10))
                
            edgestr += '['
            cyedgestr += '['
            for e in network.edges():
                edgestr += '{source:%d, target:%d}, ' %(nodesdict[e[0]], nodesdict[e[1]])
                cyedgestr += '{ data: { id: "%s", weight: 1, source: %s, target: %s } }, ' %('-'.join(e), 
                				quote(cpdnamedict.get(e[0], e[0]).split(';')[0]), quote(cpdnamedict.get(e[1], e[1]).split(';')[0]) )

            nodestr += '], '
            cynodestr += '], '
            edgestr += '], '
            cyedgestr += '], '
        
        total_d3_data = 'var nodes = [ ' + nodestr + '];\n\n        var links = [' + edgestr + '];\n\n'
        total_cytoscapejs_data = '        var cytonodes = [ ' + cynodestr + '];\n\n        var cytoedges = [' + cyedgestr + '];\n\n'
        self.jsdata = total_d3_data + total_cytoscapejs_data
        

    def rescale_color(self, dict_cpd_foldchange):
        '''
        rescale color to 0~20
        '''
        max_fc = max([abs(x) for x in dict_cpd_foldchange.values()])
        newdict = {}
        for k,v in dict_cpd_foldchange.items(): newdict[k] = int(10*v/max_fc) + 9
        return newdict

    def export_text(self):
        s = self.HTML_HEAD
        for element in self.elements:
            s += element
            
        return s + self.javascript_HEAD + self.jsdata + self.javascript_END + self.HTML_END







# --------------------------------------------------------
# Visualization using DOT and pygraphviz
#
# a few function for visualization scheme
#

def quote(s):
    return '"'+s+'"'

def flatten(inlist):
    flattened = []
    for t in inlist:
        flattened.append(t[0])
        flattened.append(t[1])
    return flattened

def draw_path(dotfile, pngfile):
    G=pgv.AGraph(dotfile)
    G.draw(pngfile, prog='dot')

def write_xlsx(data, filename, sheetname=''):
    '''
    data in regular lists, [ [], [], ... ]
    '''
    workbook = xlsxwriter.Workbook(filename)
    ws = workbook.add_worksheet(sheetname)
    row = 0
    for line in data:
        col = 0
        for item in line:
            ws.write(row, col, item)
            col += 1
        row += 1
    
    workbook.close()


def makedict_cpd_foldchange(TF):
    '''
    {cpd: max_foldchange}
    '''
    def max_fold_change(L):
        L.sort()
        if abs(L[0]) > L[-1]: return L[0]
        else: return L[-1]
        
    mydict = {}
    
    for c in TF.input_cpdlist:
        mzlist = set([x[0] for x in TF.network.cpd_dict[c].hitlist 
                        if x[0] in TF.network.input_mzlist])
        if len(mzlist) > 1:
            logging.info("    - multiple match - " + c + ' - ' + ', '.join([str(x) for x in mzlist]))
            mydict[c] = max_fold_change([TF.network.input_mzfcdict[x] for x in mzlist])
        else:
            mydict[c] = TF.network.input_mzfcdict[mzlist.pop()]
    
    return mydict


def make_color_dict(TF, zcolors):
    '''
    make a colordict {cpd:color} from m/z fold changes.
    TF is an instance of TableFeatures
    Take highest fc if conflict, print warning
    '''
        
    mydict, fcdict = {}, makedict_cpd_foldchange(TF)
    max_fc = max( np.std(TF.network.input_mzfcdict.values()), 1.0 )
    for k,v in fcdict.items():
        mydict[k] = quote(zcolors[scale_color(v, max_fc)])
    
    return mydict, make_colorbar(zcolors, max_fc)


def scale_color(f, max_fc=10):
    '''
    Auto adjusted to user data
    '''
    if f < -max_fc: return 0
    elif abs(f) < 0.001: return 5
    elif f > max_fc: return 10
    else: return 5 + int(5*f/max_fc)

def make_colorbar(zcolors, max_fc):
    '''
    Make a DOT string for colorbar placement
    subgraph colorbar {rank=source; 
        html [shape=none, margin=0, label=<
        <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" CELLPADDING="4">
        <TR>
        <TD>-1</TD>
        <TD BGCOLOR="#FFEE44">  </TD>
        <TD BGCOLOR="#FFCC44">  </TD>
        <TD BGCOLOR="#FFBB44">  </TD>
        <TD>1</TD>
        </TR>
        </TABLE>>];
        }
    '''
    right_label = str(round(max_fc, 2))
    s = 'subgraph colorbar {rank=source;\n\
            html [shape=none, margin=0, label=<\n\
            <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="1" CELLPADDING="4"><TR>'
    s += '<TD>-' + right_label + '</TD>'
    for z in zcolors: s += '<TD BGCOLOR="%s">  </TD>;'%z
    s += '<TD>' + right_label + '</TD></TR></TABLE>>];\n}\n'
    return s


#
# main class for DOT/Graphviz based visualization
#

class FishEyeViz:
    '''
    Based on previous FishEye, PathViz classes.
    Visualize a module/pathway, bipartite network of metabolites and enzymes.
    Centered on manipulating dot strings.
    Edges are imported from input module, and the network is recovered from
    reference data, by filling enzymes in between metabolite edges.
    '''
    def __init__(self, edges, idstr, edgedict, dict_cpds_def, rxnlist, enzyme=0):
        self.dir_edges = []
        self.undir_edges = []
        self.eclist = []
        self.cmpds = []
        #dicts for node formatting 
        self.labeldict = {}
        self.shapedict = {}
        self.styledict = {}
        self.colordict = {}

        self.pathlabel = idstr
        self.dict_cpds_def = dict_cpds_def
        self.medges = [x for x in edges if x[0]!=x[1]]
        self.rxndict = self.build_rxndict(rxnlist)
        
        if enzyme: self.build_bi_network(edgedict)
        else: self.build_network(edgedict)
        
        self.tally_nodes()

    def read_sif(self, sif_input):
        edges = [x.rstrip().split() for x in open(sif_input).readlines()]
        return [(x[0], x[2]) for x in edges if x[0] != x[2]]
        
    def build_rxndict(self, rxnlist):
        rxndict = {}
        for r in rxnlist: rxndict[r['id']] = r
        return rxndict
        
    def build_bi_network(self, edgedict):
        '''
        Producing a new set of edges connecting metabolites with enzymes.
        Directional and undirectional edges are distinguished.
        '''
        dir_edges, undir_edges = [], []
        for e in self.medges:
            # in case edge direction is wrong for edgedict key
            if not edgedict.has_key(e): e = (e[1], e[0])
            for r in edgedict[e].split(';'):
                if self.rxndict[r]['ecs'] and self.rxndict[r]['ecs'] != [""]:  #has enzyme
                    if self.rxndict[r]['products']:
                        # directional
                        dir_edges += self.fill_bi_edge(e, self.rxndict[r])
                    else:
                        ecstr = ';'.join(self.rxndict[r]['ecs'])
                        undir_edges += [(e[0], ecstr), (ecstr, e[1])]
                        
                else:
                    undir_edges += [e]

        self.dir_edges = dir_edges
        self.undir_edges = undir_edges
        

    def build_network(self, edgedict):
        '''
        Metabolite network without enzymes.
        Directional and undirectional edges are distinguished.
        '''
        dir_edges, undir_edges = [], []
        for e in self.medges:
            if not edgedict.has_key(e): e = (e[1], e[0])
            for r in edgedict[e].split(';'):
                rxn = self.rxndict[r]
                if e[0] in rxn['reactants'] and e[1] in rxn['products']:
                    dir_edges.append( e )
                elif e[1] in rxn['reactants'] and e[0] in rxn['products']:
                    dir_edges.append( (e[1], e[0]) )
                else:
                    undir_edges.append( e )
                
        self.dir_edges = list(set(dir_edges))
        self.undir_edges = list(set(undir_edges))


    def fill_bi_edge(self, e, rxn):
        ecstr = ';'.join(rxn['ecs'])
        if e[0] in rxn['reactants']:
            return [(e[0], ecstr), (ecstr, e[1])]
        else:
            return [(e[1], ecstr), (ecstr, e[0])]

    def tally_nodes(self):
        self.cmpds = list(set(flatten(self.medges)))
        total = set(flatten(self.dir_edges + self.undir_edges))
        self.eclist = [x for x in total if x not in self.cmpds]


    def make_header(self, mark="FishEye: "):
        size = max(8*int(np.sqrt(len(self.cmpds))), 24)
        h = 'digraph G {\n     ratio = auto; size="%d,60"; concentrate=true' %size
        h += ';\n fontsize=28; label="' + mark + ": " + self.pathlabel
        h += '";\n node [fontname=Helvetica, fontsize=18]; edge[style=bold];\n'
        return h

    # [label='', shape=ellipse, style=filled, color=".7 .3 1.0"]
    def write_nodes(self, nodelist):
        buf_str = ""
        for ec in nodelist:
            buf_str += '    "' + ec +'" ' + self.dotformat(self.labeldict[ec],
                                             self.shapedict[ec],
                                             self.styledict[ec],
                                             self.colordict[ec]
                                             ) + ';\n'
        return buf_str


    def make_a_path_dot(self, user_color_dict={}, colorbar="", mark=""):
        self.default_style()
        self.make_node_face()
        self.colordict.update(user_color_dict)
        
        outstr = self.make_header(mark)
        outstr += colorbar
        outstr += self.write_nodes(self.eclist + self.cmpds)

        for edge in self.dir_edges:
            outstr += self.makeline(edge)
        for edge in self.undir_edges:
            outstr += self.makeline(edge, "dir=none")

        outstr += "    }\n"
        return outstr
    
    def write_dot(self, user_color_dict, colorbar, mark,  f='fisheye_viz_result.dot'):
        out = open(f, 'w')
        out.write( self.make_a_path_dot(user_color_dict, colorbar, mark))
        out.close()
    
    def dotformat(self, label, shape, style, color):
        buf_str = '[label='+label + ', shape='+shape
        if style:
            buf_str += ', style='+style
        if color:
            buf_str += ', color='+color
        return buf_str + ']'

    def default_style(self):
        for ec in self.eclist:
            self.labeldict[ec] = quote(ec)
            self.shapedict[ec] = 'plaintext'
            self.styledict[ec] = ''
            self.colordict[ec] = '' #".7 .3 1.0"
        for cp in self.cmpds:
            self.labeldict[cp] = quote(cp)
            self.shapedict[cp] = 'ellipse'
            self.styledict[cp] = 'filled'
            self.colordict[cp] = '"#FFFF99"'


    def make_node_face(self):
        '''
        labeldict is the face value or node label in the dot png
        '''
        for c in self.cmpds:
            self.labeldict[c] = self.beautifycpd(c)

    def makeline(self, (node1, node2), attr=""):
        '''
        write a line in DOT file, node1 -> node2
        '''
        return '    "' + node1 + '" -> "' + node2 + '" [' + attr + '];\n'
    
    def beautifycpd(self, cpd):
        '''
        adjust long compound names to multiple lines
        '''
        try:
            s = self.dict_cpds_def[cpd].split(";")[0]
        except KeyError:
            s = cpd
        slen = len(s)
        spare = [ ]
        chunknum = slen/20 + 1
        if chunknum > 2:
            return quote(cpd)
        elif chunknum > 1:
            blocks = s.split('-')
            newstr = blocks[0]
            for b in blocks[1:]:
                if len(newstr) < 23:
                    newstr += '-' + b
                else:
                    spare.append(newstr)
                    newstr = b
            spare.append(newstr)
            return quote("-\\n".join(spare))
        else:
            return quote(s)


