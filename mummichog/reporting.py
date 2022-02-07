# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
Generate reports from mummichog analysis.
a) Summary report in HTML
b) Result tables
c) Intermediate data
d) Web based presentation

@author: Shuzhao Li, Francisco Castellanos
With js code contributed by Andrei Todor

'''

import os
import shutil
import csv
import json
import xlsxwriter
import logging
import numpy as np

from . import resources
from .websnippets import *
from .config import VERSION, SIGNIFICANCE_CUTOFF

class WebReporting:
    '''
    To generate HTML summary report.
    
    This will change to template (Jinja2) based methods.
    '''
    def __init__(self, Local, PA, MA, AN):

        self.mixedNetwork = Local.mixedNetwork
        self.model = Local.mixedNetwork.model
        self.data = Local.mixedNetwork.data
        self.Local = Local
        
        self.PA = PA
        self.MA = MA
        self.AN = AN

        logging.info('Web reporting')
        
    def run(self):
        self.get_dict_cpd_statistic()
        self.collect_web_export_graphs()
        html = self.web_export()
        return html
        
    def run_nodisk(self):
        '''
        For web use, no local disk write
        '''
        self.get_dict_cpd_statistic()
        self.collect_web_export_graphs()
        return self.web_export_str()
        
    def get_dict_cpd_statistic(self):
        self.dict_cpd_statistic = {}
        for T in self.mixedNetwork.hit_Trios:
            self.dict_cpd_statistic[T[2]] = self.mixedNetwork.rowDict[T[0]].statistic
        
    def collect_web_export_graphs(self):
        '''
        Return activity network and up to top 5 modules.
        Moving activity network after modules
        
        AN yet to be revised
        
        '''
        self.web_export_networks = [M.graph for M in self.MA.top_modules[:5]] + [self.AN.activity_network]
        
        
        vis_nodes = []
        for g in self.web_export_networks: vis_nodes += g.nodes()
        self.vis_nodes = set(vis_nodes)


    def web_export(self):
        '''
        Return HTML without head/foot.
        '''
        outfile = os.path.join(self.Local.rootdir, 'result.html')
        webtextList = self.web_export_str()
        with open(outfile, 'w') as O:
            O.write(''.join(webtextList))
        return webtextList[1]

    def web_export_str(self):
        '''
        Write HTML and javascript based report.
        The visualization section includes activity network and up to top 5 modules.
        return webtextList as [head, middle HTML, end]
        '''

        HTML = HtmlExport(self.data, self.PA)
        title = 'Mummichog Report: ' + self.data.paradict['output']
        HTML.add_element(title, 'h1', '')
        
        HTML.add_element('', 'div', 'colorbar', 'colorbar')
        stats = "mummichog version: %s, using metabolic model %s.<br>" %(VERSION, self.model.version)
        HTML.add_element(stats, 'div', 'stats')
        
        # MWAS figure
        HTML.add_element('User input data', 'h2', '')
        details_userData = "User supplied %d features as reference list." %len(self.data.ListOfMassFeatures)
        details_userData += "We are using %d features (p < %f) as significant list. The feature level data are shown in the Manhattan plots below." %(
                                    len(self.data.input_featurelist), self.data.paradict['cutoff'] )
        HTML.add_element(details_userData, 'p', '')
        # HTML.add_element(self.Local.inline_plot_userData_MWAS, 'div', 'inline_plot_userData_MWAS') --  Using js plots instead
        HTML.add_element('', 'div', '', 'mz_user_input')
        HTML.add_element('', 'div', '', 'retention_time_input')


        # Pathway table and figure
        HTML.add_element('Top pathways', 'h2', '')
        details_pathwayData = "The reference feature list is mapped to %d EmpiricalCompounds, significant features to %d EmpiricalCompounds. " %(
                                    self.PA.total_number_EmpiricalCompounds, len(self.PA.significant_EmpiricalCompounds))
        details_pathwayData += "In the table below, pathway_size is number of detected EmpiricalCompounds for each pathway; overlap_size is number of significant EmpiricalCompounds.\
                                Empirical p-values are estimated by permutation test. Details on EmpiricalCompounds are in 'tables/ListOfEmpiricalCompounds.tsv'."
        HTML.add_element(details_pathwayData, 'p', '')
        pathwaystablly = self.write_pathway_table()
        HTML.add_element(pathwaystablly, 'div', 'pathwaystablly')
        # HTML.add_element(self.Local.inline_plot_pathwayBars, 'div', 'inline_plot_pathwayBars') --  Using js plots instead
        HTML.add_element('', 'div', '', 'pathwayBarPlot')

        # place to insert network visusalization
        HTML.add_element('Top modules', 'h2', '')
        Module_pharase = "Network modules may capture concerted metabolite activities that are missed by predefined pathways. \
        Use pull-down menus to visualize modules and the activity network, which tries to combine the results of \
        pathway/module analyses. The visualization requires interenet connection."
        HTML.add_element(Module_pharase, 'p', '')
        select_menu = HTML.make_select_menu(len(self.web_export_networks))
        HTML.add_element(select_menu, 'div', 'network_selection', 'network_selection')
        # another color bar?
        #HTML.add_element('', 'div', 'colorbar', 'colorbar')
        HTML.add_element('', 'div', 'networkvisual', 'networkvisual')
        

        modulestablly = self.write_module_table()
        HTML.add_element(modulestablly, 'div', 'modulestablly')
        
        HTML.add_element('EmpiricalCompounds favored by above analyses', 'h2', '')
        EmpCpd_phrase = "An EmpiricalCompound is a computational unit for a tentative metabolite. It can group multiple ions, \
                    and be any of the isobaric/isomeric species. This table contains EmpiricalCompounds prioritized by pathway/module analysis."
        HTML.add_element(EmpCpd_phrase, 'p', '')
        metabolitestablly = self.write_metabolite_table()
        HTML.add_element(metabolitestablly, 'div', 'metabolitestablly')
        
        HTML.add_element('More data', 'h2', '')
        moreinformation = '''Full data tables from pathway analysis, module analysis and activity network are stored in the 'tables/' directory. 
            The tsv files can be imported into a spreadsheet program, e.g. MS Excel.
            Under the 'figures/network_modules/' directory are files intended for Cytoscape (cytoscape.org) visualization. Please refer to Cytoscape's guides for details.
            Details of this run are recorded in file mummichog.log.txt.
            '''
        HTML.add_element(moreinformation, 'p', 'moreinformation')
        
        footer = '''Mummichog algorithms are described in Li et al. Predicting Network Activity from High Throughput Metabolomics. PLoS Computational Biology (2013); doi:10.1371/journal.pcbi.1003123.\
        This software is provided as is, without warranty of any kind. Visit <a href="http://mummichog.org">http://mummichog.org</a> for updates.
        '''
        HTML.add_element(footer, 'footer', '')
        
        # moved to self.export_sif_related()
        # self.export_cpd_attributes()
        
        # push network data into javascript
        HTML.make_js_data(self.web_export_networks, 
                          self.model.dict_cpds_def, 
                          self.dict_cpd_statistic)
        
        return HTML.export_text()
        
    def filter_vis_nodesdict(self, d1):
        d2 = {}
        for k in d1.keys():
            if k in self.vis_nodes:
                d2[k] = d1[k]
                
        return d2

    def write_pathway_table(self):
        s = '<table><tr><th>Pathways</th>\
                    <th>overlap_size</th><th>pathway_size</th><th>p-value</th>\
                    <th>overlap_EmpiricalCompounds</th></tr>'
        ii = 0
        for P in self.PA.resultListOfPathways:
            ii += 1
            if ii < 6 or P.adjusted_p < SIGNIFICANCE_CUTOFF:
                empCpds = [E.EID for E in P.overlap_EmpiricalCompounds]
                s += '<tr> <td>' + P.name + '</td> <td>'\
                        + str(P.overlap_size) + '</td><td>' + str(P.EmpSize
                        ) + '</td><td>' + str(round(P.adjusted_p, 5)) + '</td><td>' \
                        + ','.join(empCpds) + '</td></tr>'        

        return s + '</table>'

    
    def write_module_table(self):
        '''
        for M in self.MA.top_modules:
            counter += 1
            M.export_network_txt(self.model, 
                                 os.path.join(self.moduledir, 'module_' + str(counter) + '.txt'))
        
        '''
        s, counter = '', 0
        for M in self.MA.top_modules:
            counter += 1
            nodes = M.graph.nodes()
            names = [self.model.dict_cpds_def.get(x, '') for x in nodes]
            s += '<div class="moduleline">' + 'module_' + str(counter) + ", p=" + str(round(M.p_value, 5)) + ", " + str(len(nodes)) + " metabolites" + '</div>'
            s += '<div class="metabolites">' + ', '.join(names) + '</div>'
            
        return s + '\n'



    def write_metabolite_table(self):
        '''
        To export significant EmpiricalCompounds that has a face_compound, which was added during pathway/module analysis
        
        
        '''
        
        s = '<table><tr><th>EmpiricalCompound</th><th>Input m/z</th>\
                    <th>Retention time</th><th>ion</th><th>mz_diff</th><th>Statistic</th>\
                    <th>Significant</th></tr>'
        # will add theoretical MW?
        for E in self.PA.significant_EmpiricalCompounds:
            if E.face_compound:
                col2 = 'Best guess: ' + E.face_compound + ', ' + self.model.dict_cpds_def.get(E.face_compound, '').split(";")[0]
                s += '<tr > <td>' + E.EID + '</td> <td colspan="2">' + col2 + '</td><td colspan="4">' + ', '.join(E.compounds) + '</td></tr>'
                for f in E.massfeature_rows:
                    F = self.mixedNetwork.rowDict[f]
                    ion = E.row_to_ion[f]
                    
                    # should store earlier?
                    mz_diff = F.mz - self.model.Compounds[ E.face_compound ]['adducts'][ ion ]
                    
                    s += '<tr> <td> </td> <td>' + str(F.mz) + '</td> <td>' + str(F.retention_time) + '</td><td>' + ion + '</td><td>' + str(
                                round(mz_diff,4)) + '</td><td>' + str(round(F.statistic,2)) + '</td><td>' + write_yes_no_MassFeature(F) + '</td></tr>'
           
        return s + '</table>'
                


class LocalExporting:
    '''
    To write out local data

    '''
    def __init__(self, mixedNetwork, PA, MA, AN):
        self.mixedNetwork = mixedNetwork
        self.model = mixedNetwork.model
        self.data = mixedNetwork.data
        
        self.PA = PA
        self.MA = MA
        self.AN = AN
        
        self.create_dirs()
        
    def create_dirs(self):
        '''
        Directory tree, renamed for version 2
        outdir/ result.html
            mummichog.log.txt
            tables/
            figures/
                network_modules/ 
        
        logging.basicConfig(filename=os.path.join(optdict['workdir'], optdict['outdir'], 'mummichog.log'), 
                            format='%(message)s', 
                            level=logging.INFO)
        logging.info('\n'.join(["mummichog version: %s" %VERSION,
                                 "pwd: %s" %os.getcwd(),
                                 "user command: %s" %' '.join(sys.argv),
                                 "\n",
                                 ]))
        
        
        #print_and_loginfo("Pathway analysis report was written to \n%s (.tsv and .xlsx)" %outfile)
     


        Web export may avoid static file linking, as it's tricky to be correct on local computers
        '''
        self.rootdir = os.path.join(self.data.paradict['workdir'], self.data.paradict['outdir'])
        os.mkdir(self.rootdir)
        self.tabledir, self.figuredir, self.moduledir = os.path.join(
                                                self.rootdir, 'tables'), os.path.join(
                                                self.rootdir, 'figures'), os.path.join(
                                                self.rootdir, 'figures', 'network_modules')
        os.mkdir(self.tabledir)
        os.mkdir(self.figuredir)
        os.mkdir(self.moduledir)

        self.jsDir = os.path.join(self.rootdir, 'js')
        os.mkdir(self.jsDir)

    def run(self):
        '''
        '''
        # export tables
        self.export_userData()
        self.export_EmpiricalCompounds()
        self.export_pathway_enrichtest()
        self.writeTable_top_modules()
        
        # plot figures
        # self.plot_userData_MWAS() commenting out as plots are available via plotly js and rendered in the client browser
        # self.plot_pathwayBars() commenting out as plots are available via plotly js and rendered in the client browser
        self.plot_pathway_model()
        self.plot_module_model()
        
        # export mudules and activity network
        self.export_top_modules()
        self.export_activity_network()
        self.export_cpd_attributes()
        self.exportJsLibs()

        # testing JSON output
        self.export_json_all()
        
        
    def export_userData(self):
        '''
        to write out user data with row_numbers,
        to use as unique ID and help users to track data.
        '''
        s = "massfeature_rows\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user\n"
        for F in self.data.ListOfMassFeatures:
            s += F.make_str_output() + '\n'
                
        with open(os.path.join(self.tabledir, "userInputData.txt"), 'w') as O:
            O.write(s)
            

    def export_EmpiricalCompounds(self):
        '''
        This exports all tri relationships.
        In ActivityNetwork, the top predicted metaolite is determined.

        Output two files: empCpd to row#, row# to empCpd
        '''
        s = "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names\n"
        for E in self.mixedNetwork.ListOfEmpiricalCompounds:
            names = [self.model.dict_cpds_def.get(x, '') for x in E.compounds]
            s += '\t'.join([E.EID, ';'.join(E.massfeature_rows), E.str_row_ion, ';'.join(E.compounds), '$'.join(names)]
                ) + '\n'
        with open(os.path.join(self.tabledir, "ListOfEmpiricalCompounds.tsv"), 'w') as O:
            O.write(s)

        s = "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user\n"
        for row in self.mixedNetwork.mzrows:
            # not all input rows match to an empCpd
            try:
                for E in self.mixedNetwork.rowindex_to_EmpiricalCompounds[row]:
                    names = [self.model.dict_cpds_def.get(x, '') for x in E.compounds]
                    s += '\t'.join([row, E.EID, E.str_row_ion, ';'.join(E.compounds), '$'.join(names)]
                        ) + '\t' + self.mixedNetwork.rowDict[row].make_str_output() + '\n'
            except KeyError:
                pass

        with open(os.path.join(self.tabledir, "userInput_to_EmpiricalCompounds.tsv"), 'w') as O:
            O.write(s)

        with open(os.path.join(self.tabledir, 'exported_Compounds.json'), 'w', encoding='utf-8') as f:
            json.dump(self.model.Compounds, f, ensure_ascii=False, indent=2)
    
    # export functions for pathway analysis
    def export_pathway_enrichtest(self):
        '''
        result as sorted list of pathway instances
        Cutting out 'p-value (raw)' from version 1.
        '''
        resultstr = [['pathway', 'overlap_size', 'pathway_size', 'p-value', 
                      'overlap_EmpiricalCompounds (id)', 'overlap_features (id)', 'overlap_features (name)',] ]
        
        for P in self.PA.resultListOfPathways:
            empCpds = [E.EID for E in P.overlap_EmpiricalCompounds]
            cpds = [E.chosen_compounds for E in P.overlap_EmpiricalCompounds]
            names = [ [self.model.dict_cpds_def.get(x, '') for x in y] for y in cpds ]
            resultstr.append([str(x) for x in [P.name, P.overlap_size, P.EmpSize, P.adjusted_p]]
                             + [','.join(empCpds), ','.join(['/'.join(x) for x in cpds]), '$'.join(['/'.join(x) for x in names]) ])

        outfile = os.path.join(self.tabledir, "mcg_pathwayanalysis_") + self.data.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', "w", newline="")
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_pathwayanalysis_")
        

    # export functions for module analysis
    def writeTable_top_modules(self):
        '''
        Write .tsv report for top modules.
        '''
        resultstr = [
                     ['MODULE', 'p-value', 'size', 'members (name)', 'members (id)', 
                      'This module overlaps with'], ]
        
        for M in self.MA.top_modules:
            nodes = M.graph.nodes()
            names = [self.model.dict_cpds_def.get(x, '') for x in nodes]
            resultstr.append([M.id, str(M.p_value),
                              str(len(nodes)), '$'.join(names),
                              ','.join(nodes), self.find_top_pathways(nodes),
                               ])

        outfile = os.path.join(self.tabledir, "mcg_modularanalysis_") + self.data.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', "w", newline="")
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_modularanalysis_")
        
    def find_top_pathways(self, nodes):
        '''
        return top 3 pathways for each feature, used for annotation
        '''
        pathways = []
        for n in nodes: pathways += self.model.cpd2pathways.get(n, [])
        counts = [(pathways.count(p), p) for p in set(pathways)]
        counts.sort(reverse=True)
        return ';'.join([str(x) for x in counts[:3]])


    def plot_userData_MWAS(self):
        self.inline_plot_userData_MWAS = self.data.make_manhattan_plots( 
            os.path.join(self.figuredir, "mcg_MWAS_") + self.data.paradict['output']
            )

    def plot_pathwayBars(self):
        self.inline_plot_pathwayBars = self.PA.plot_bars_top_pathways(
            os.path.join(self.figuredir, "mcg_pathwayBars_") + self.data.paradict['output']
            )

    def plot_pathway_model(self):
        self.PA.plot_model_pvalues(
            os.path.join(self.figuredir, "plot_pathwayModel_") + self.data.paradict['output']
            )

    def plot_module_model(self):
        self.MA.plot_model_pvalues(
            os.path.join(self.figuredir, "plot_moduleModel_") + self.data.paradict['output']
            )

    def draw_top_modules(self):
        '''
        May add later
        '''
        pass

    def export_top_modules(self):
        '''
        export to figures/network_modules/
        for visualizatin using CytoScape
        '''
        for M in self.MA.top_modules:
            M.export_network_txt(self.model, 
                                 os.path.join(self.moduledir, M.id + '.txt'))

    def export_activity_network(self):
            self.AN.export_network_txt(self.model, 
                                os.path.join(self.moduledir, 'activity_network.txt'))

    
    def export_cpd_attributes(self):
        '''
        Only for significant Trios.
        This exports attributes for network visualization in Cytoscape 3.
        
        mixedNetwork.hit_Trios is updated after initiation of ActivityNetwork
        [(mzFeature, EmpiricalCompound, cpd),...]
        
        '''
        s = 'CPD_ID\tCPD_NAME\tmz_row_Statistic\tFrom_mz_row\tFrom_EmpiricalCompound\n'
        for T in self.mixedNetwork.hit_Trios:
            name = self.model.dict_cpds_def.get(T[2], '')
            mz_row_Statistic = self.mixedNetwork.rowDict[T[0]].statistic
            s += '\t'.join([T[2], name, str(mz_row_Statistic), T[0], T[1].EID]) + '\n'
            
        with open(os.path.join(self.moduledir, "Node_attributes.txt"), 'w') as O:
            O.write(s)
    
    
    def export_network(self, network):
        '''
        Write ActivityNetwork for Cytoscape 3 formats.
        Changed to .txt in mummichog v 1.0.
        '''
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in network.edges():
            s += e[0] + '\t' + e[1] + '\t' + self.network.MetabolicModel.edge2enzyme.get(e, '') + '\n'
        
        out = open(os.path.join(self.sifdir, 'ActivityNetwork.txt'), 'w')
        out.write(s)
        out.close()

    def exportJsLibs(self):
        '''
        To-do in v3: Should move all js visualization to separate package.
        Difficult to copy files. os.getcwd() is a wrong way 
        because mummichog may be installed to Python site-packages as library.
        '''
        package_dir = resources.__path__[0]
        # Use os.path.join to be OS independent.
        shutil.copyfile(os.path.join(package_dir, "plotly-graphs.js"), 
                        os.path.join(self.jsDir, "plotly-graphs.js"))
        shutil.copyfile(os.path.join(package_dir, "plotly-latest.min.js"), 
                        os.path.join(self.jsDir, "plotly-latest.min.js"))

    def export_json_all(self):
        '''
        This will be the export of all mummichog output in JSON,
        thus enabling visualization in a separate package.

        Manhattan plots are based on userData_original.
        userData_empirialCompound is link btw user input data and theoretical metabolites.

        Prototyping, will improve details [2020-04-21, SL]
        '''
        json_data = {

            'userData_original' : [
                {
                    'id': M.row_number, 'mz': str(M.mz), 'retention_time': str(M.retention_time), 'p_value': str(M.p_value),
                    'statistic': str(M.statistic), 'CompoundID_from_user': M.CompoundID_from_user,
                } for M in self.data.ListOfMassFeatures
            ],

            'userData_empirialCompound' : [
                {
                    'id': E.EID,
                    'listOfFeatures': E.listOfFeatures,     # This links to userData_original['id'], which is row_number
                    'chosen_compounds': E.chosen_compounds,
                    'face_compound': E.face_compound,
                    # will add compound names, def, etc.
                } for E in self.mixedNetwork.ListOfEmpiricalCompounds
            ],

            'meta_data' : {
                'input_parameters': self.data.paradict,
                'metabolic_model': self.model.version,
                'mummichog_version': VERSION,
                'significance_cutoff' : SIGNIFICANCE_CUTOFF
            },

            # a pathway is defined in models.metabolicPathway
            'result_pathwayAnalysis' : [
                {
                    'id': P.id,
                    'name': P.name,
                    'pathway_p_value': P.adjusted_p,
                    'overlap_size': P.overlap_size,
                    'pathway_size': P.EmpSize,
                    'overlap_EmpiricalCompounds': [E.EID for E in P.overlap_EmpiricalCompounds],
                    'significant_compounds': [E.chosen_compounds for E in P.overlap_EmpiricalCompounds],
                    # 'all_compounds': P.cpds,
                } for P in self.PA.resultListOfPathways
            ],

            'result_networkModules' : [
                {
                    'id': M.id,
                    'p_value': M.p_value,
                    'members': M.graph.nodes(), # This is on theoretical metabolites/compounds
                    #'size' is number of members
                    # Common names of members can be looked up, from metabolic model/metabolite definition
                } for M in self.MA.top_modules
            ],

        }

        # testing dumps
        # with open(os.path.join(self.rootdir, "result.json"), "w") as O:
        #    O.write(json.dumps(json_data))
        
        with open(os.path.join(self.jsDir, "result.js"), "w") as O:
            O.write("var alldata = '" + json.dumps(json_data) + "';")


class HtmlExport:
    '''
    This class is not tied to any data. 
    Serving as a skeleton to be used for export functions in AnalysisCentral.
    In future versions, this can use existing tools to convert JSON to HTML.
    '''
    def __init__(self, data, PA):
        self.elements = []
        self.jsdata = ''
        
        self.HTML_HEAD = HTML_HEAD
        self.HTML_END = HTML_END
        self.javascript_HEAD = javascript_HEAD
        self.javascript_END = javascript_END
        self.data = data
        self.PA = PA
        
        
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


    def make_js_data(self, networks, cpdnamedict, dict_cpd_statistic):
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

        cpdcolordict = self.rescale_color(dict_cpd_statistic)
        
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
    
    def rescale_color_testing(self, dict_cpd_statistic):
        '''
        rescale color to 0~20
        7 steps, centered at 0, +- 3 * stdev
        '''
        stdev = np.std(dict_cpd_statistic.values())
        step = stdev
        newdict = {}
        for k,v in dict_cpd_statistic.items(): 
            if v > 3*step:
                newdict[k] = 9+9
            elif v < -3*step:
                newdict[k] = 1
            else:
                newdict[k] = int(v/step) + 9
        
        return newdict

    def export_text(self):
        s = ''
        for element in self.elements:
            s += element

        return [self.HTML_HEAD, 
                s + self.javascript_HEAD + self.jsdata + self.javascript_END,
                self.HTML_END]



# --------------------------------------------------------
#
# a few function for utility or visualization scheme
#

def write_xlsx(data, filename, sheetname=''):
    '''
    data in regular lists, [ [], [], ... ]
    '''
    workbook = xlsxwriter.Workbook(filename, {'strings_to_numbers': True})
    ws = workbook.add_worksheet(sheetname)
    row = 0
    for line in data:
        col = 0
        for item in line:
            ws.write(row, col, item)
            col += 1
        row += 1
    
    workbook.close()

def write_yes_no_MassFeature(f):
    if f.is_significant:
        return "yes"
    else:
        return "no"


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
