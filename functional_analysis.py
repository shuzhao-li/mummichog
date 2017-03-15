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

@author: Shuzhao Li, Andrei Todor
'''

import logging, random
from scipy import stats
import ng_modularity as NGM


#from base import *
from get_user_data import *

from reporting import *

logging.basicConfig(format='%(message)s', level=logging.INFO)



class AnalysisCentral:
    '''
    Instance to organize analysis data;
    separation of data and presentation from functional units.
    
    MetabolicNetwork and input_mzlist are processed in TableFeatures.
    Visualization is transitioning to web based platform.
    
    '''

    def __init__(self, mixedNetwork):
        '''
        working on v2
        mixedNetwork = DataMeetModel(theoreticalModel, userData)
        self.
        rowDict
        cpd2mzFeatures
        ListOfEmpiricalCompounds
        '''
        
        self.mixedNetwork = mixedNetwork
        self.model = mixedNetwork.model
        self.data = mixedNetwork.data
        self.create_dirs()


    def create_dirs(self):
        '''
        Directory tree
        outdir_result.html
        tsv/ result data in tsv, xls files
        sif/ all sif visualization related files
        web/
        '''
        self.rootdir = os.path.join(self.data.paradict['workdir'], self.data.paradict['outdir'])
        self.csvdir, self.sifdir, self.webdir = os.path.join(
                                                self.rootdir, 'tsv'), os.path.join(
                                                self.rootdir, 'sif'), os.path.join(
                                                self.rootdir, 'web')
        os.mkdir(self.csvdir)
        os.mkdir(self.sifdir)
        os.mkdir(self.webdir)

    def run_all_analysis(self):
        '''
        Separation of analysis and presentation
        
        Good place to debug,
        by commenting out some steps 
        
        
        NI = NetworkInspector(self)
        self.Inspected = NI.inspect_network()
        self.mwstr2mnodes = NI.mwstr2mnodes
        
        
        '''
        metabolic_pathways = self.model.metabolic_pathways
        
        PA = PathwayAnalysis(metabolic_pathways, self.mixedNetwork)
        self.pathway_result = PA.cpd_enrich_test()
        
        MA = ModularAnalysis(self.mixedNetwork)
        MA.dispatch()
        self.top_modules = MA.top_modules
        
        
        

    def export_csv_data(self):
        '''
        Write all result data to tsv/
        '''
        self.write_tentative_annotation()
        self.write_pathway_enrichtest_result()
        self.write_top_modules()
        self.export_Mnodes_evidence()
        self.export_metabolite_worksheet()
        
    def export_sif_related(self):
        '''
        graphviz drawing is deprecated
        '''
        self.export_top_modules()
        self.export_network(self.Inspected)
        self.export_dict_cpd_statistic()


    def write_tentative_annotation(self):
        '''
        Write tentative annotation file for reference features,
        as iterated through CpdList        
        C.hitlist = [(input mz, match form, diff), ...]
        
        
        ver 2
        Changing data structure
        Initial map btw input features and Cpd
        
        '''
        resultstr = []
        for C in self.network.CpdList:
            for x in C.hitlist:
                resultstr.append([x[0].mz, x[0].retention_time, C.id, x[1], x[2], C.name] + [
                                        '$'.join(self.network.MetabolicModel.cpd2pathways.get(C.id, []))])
        
        resultstr.sort()
        data2write = [ ['m/z', 'rt', 'id', 'match_form', 'mz_difference', 'name', 'pathway']
                     ] + resultstr
        outfile = os.path.join(self.csvdir, "_tentative_featurematch_"
                               ) + self.network.paradict['output']
        
        # write .tsv
        f =  open(outfile + '.tsv', 'wb')
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( data2write )
        f.close()
        # export .xlsx too
        write_xlsx(data2write, outfile + '.xlsx', "_tentative_featurematch_")
        
        print_and_loginfo("\nAnnotation was written to \n%s (.tsv and .xlsx)" %outfile)


    #
    # export functions for pathway analysis
    #
    def write_pathway_enrichtest_result(self):
        '''
        result as sorted list of pathway instances
        '''
        resultstr = [['pathway', 'overlap_size', 'pathway_size', 'p-value (raw)',
                      'p-value', 'overlap_features (name)', 'overlap_features (id)', 'used_input_mzs'], ]
        
        usenodes = []   # to bag nodes from significant pathways
        for P in self.pathway_result:
            nodes = P.overlap_features
            if P.adjusted_p < 0.05:
                usenodes += list(nodes)
            used_input_mzstr = []
            for n in nodes:
                used_input_mzstr.append([str(x[0]) for x in self.network.cpd_dict[n].hitlist 
                                                    if x[0] in self.network.input_featurelist])
                
            names = [ self.network.MetabolicModel.dict_cpds_def.get(x, '') for x in nodes ]
            resultstr.append([str(x) for x in [P.name, P.overlap_size, P.cpd_num, 
                                               P.p_FET, P.adjusted_p]]
                             + ['$'.join(names), ';'.join(nodes), ';'.join([','.join(x) for x in used_input_mzstr])])

        resultstr.append([])
        resultstr.append(['Annotation for metabolites in significant pathways'])
        resultstr.append(['m/z', 'rt', 'id', 'match_form', 'mz_difference', 'name', 'pathway'])
        usenodes = set(usenodes)
        for C in self.network.CpdList:
            if C.id in usenodes:
                for x in C.hitlist:
                    resultstr.append( [x[0].mz, x[0].retention_time, C.id, x[1], x[2], C.name] + [
                                            '$'.join( self.network.MetabolicModel.cpd2pathways.get(C.id, []) )] )
        
        
        outfile = os.path.join(self.csvdir, "mcg_pathwayanalysis_") + self.network.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', 'wb')
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_pathwayanalysis_")
        
        print_and_loginfo("Pathway analysis report was written to \n%s (.tsv and .xlsx)" %outfile)
        
    
    # not used now
    def draw_color_pathways(self, result):
        '''
        draw pathways with p < SIGNIFICANCE_CUTOFF or top 10 pathways
        '''
        use_result = [r for r in result if r.adjusted_p < SIGNIFICANCE_CUTOFF]
        if len(use_result) < 10: use_result = result[:10]
        
        print_and_loginfo("\nDrawing pathway maps in %s/..." %self.csvdir)
        for P in use_result:
            self.draw_colorpath(self.network.MetabolicModel.pathdotdict[P.id], P.overlap_features, 
                                os.path.join(self.csvdir, 
                                        self.make_png_name(P.name, P.adjusted_p)))
        
        
    def make_png_name(self, name, pvalue):
        return str(-np.log10(pvalue))[:4] + '_' + name.replace('/', '') + '.png'
        
    def draw_png(self, dotstr, pngfile):
        G=pgv.AGraph(dotstr)
        G.draw(pngfile, prog='dot')
        
    def make_green(self, dotline):
        return dotline.replace('];', 
                        ', shape="egg", style="filled", color="lawngreen"];\n')
        
    # not used now
    def draw_colorpath(self, dotstr, featurelist, outfile):
        new = []
        for line in dotstr.splitlines():
            #avoid edge lines marked by '[]'
            if '[]' not in line and set(line.split('[')[0].replace('"', '').strip(
                                        ).split(';')).intersection(set(featurelist)):
                newline = self.make_green(line)
                new.append( newline )   
            else:
                new.append(line)
                    
        self.draw_png(''.join(new), outfile)
        


    #
    # export functions for module analysis
    #

    def write_top_modules(self):
        '''
        Write .tsv report for top modules.
        '''
        resultstr = [ ['Selected modules (p < %s)'%SIGNIFICANCE_CUTOFF],
                     ['MODULE', 'p-value', 'size', 'members (name)', 'members (id)', 
                      'top pathways', 'members'], ]
        allnodes, counter = [], 0
        for M in self.top_modules:
            counter += 1
            nodes = M.graph.nodes()
            names = [self.network.MetabolicModel.dict_cpds_def.get(x, '') for x in nodes]
            resultstr.append(['module_' + str(counter), str(M.p_value),
                              str(len(nodes)), '$'.join(names),
                              ','.join(nodes), self.find_top_pathways(nodes),
                              '$'.join([str(x) for x in zip(nodes, names)]) ])
                              
            allnodes += nodes
        
        resultstr.append([])
        resultstr.append(['Annotation for metabolites in above modules'])
        resultstr.append(['m/z', 'rt', 'id', 'match_form', 'mz_difference', 
                          'name', 'pathway'])
        allnodes = set(allnodes)
        for C in self.network.CpdList:
            if C.id in allnodes:
                for x in C.hitlist:
                    resultstr.append( [x[0].mz, x[0].retention_time, C.id, x[1], x[2], C.name] + [
                                            '$'.join( self.network.MetabolicModel.cpd2pathways.get(C.id, []) )] )
        
        outfile = os.path.join(self.csvdir, "mcg_modularanalysis_") + self.network.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', 'wb')
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_modularanalysis_")
        
        print_and_loginfo("\nModular analysis report was written to\n%s (.tsv and .xlsx)" %outfile)



    def find_top_pathways(self, nodes):
        '''
        return top 3 pathways for each feature, used for annotation
        '''
        pathways = []
        for n in nodes: pathways += self.network.MetabolicModel.cpd2pathways.get(n, [])
        counts = [(pathways.count(p), p) for p in set(pathways)]
        counts.sort(reverse=True)
        return ';'.join([str(x) for x in counts[:3]])

    def find_match_pathway(self, nodes):
        return [n+'('+str(self.network.MetabolicModel.cpd2pathways.get(n, []))+')' for n in nodes]


    def draw_top_modules(self):
        '''
        Will add m/z values per cpd in visualization?
        To add Scale bar
        '''
        print_and_loginfo("\nDrawing top module graphs in %s/..." %self.modules_dir)
        mycdict, colorbar = make_color_dict(self.tf, zcolors)
        ii = 0
        for M in self.top_modules:
            ii += 1
            self.draw_module(M, 'module_' + str(ii), mycdict, colorbar)


    def draw_module(self, M, idstr, mycdict, colorbar):
        '''
        Draw module pictures with and without enzymes;
        output to dir tsv/
        '''
        apath = FishEyeViz(M.graph.edges(), idstr, self.network.MetabolicModel.edge2rxn, 
                           self.network.MetabolicModel.dict_cpds_def,
                           self.network.MetabolicModel.metabolic_rxns, 1)
        dotfile = os.path.join(self.csvdir, idstr+'.dot')
        pngfile = os.path.join(self.csvdir, idstr+'.png')
        apath.write_dot(mycdict, colorbar, 'mummichog', dotfile)
        draw_path(dotfile, pngfile)
        
        # no enzyme
        apath = FishEyeViz(M.graph.edges(), idstr, self.network.MetabolicModel.edge2rxn, 
                           self.network.MetabolicModel.dict_cpds_def,
                           self.network.MetabolicModel.metabolic_rxns, 0)
        dotfile = os.path.join(self.csvdir, idstr+'_noenzyme.dot')
        pngfile = os.path.join(self.csvdir, idstr+'_noenzyme.png')
        apath.write_dot(mycdict, colorbar, 'mummichog', dotfile)
        draw_path(dotfile, pngfile)
        

    def export_top_modules(self):
        '''
        export to sif/
        for visualizatin using CytoScape
        '''
        print_and_loginfo("\nExporting top modules to %s/..." %self.sifdir)
        ii = 0
        for M in self.top_modules:
            ii += 1
            M.export_network_txt(self.network.MetabolicModel, 
                          os.path.join(self.sifdir, 'module_' + str(ii) + '.txt'))
            

    #
    # export functions for activity network analysis
    #

    def draw_network(self, network):
        '''
        Draw a single figure of activity network.
        Remaking the color scheme will yield new warnings of multiple matches.
        '''
        print_and_loginfo("\nDrawing ActivityNetwork graph in %s/..." 
                          %self.sifdir)
        mycdict, colorbar = make_color_dict(self, zcolors)

        apath = FishEyeViz(network.edges(), 'ActivityNetwork', 
                           self.network.MetabolicModel.edge2rxn, 
                           self.network.MetabolicModel.dict_cpds_def,
                           self.network.MetabolicModel.metabolic_rxns, 1)
        dotfile = os.path.join(self.sifdir, 'ActivityNetwork.dot')
        pngfile = os.path.join(self.sifdir, 'ActivityNetwork.png')
        apath.write_dot(mycdict, colorbar, 'mummichog', dotfile)
        draw_path(dotfile, pngfile)
        
    
    def export_dict_cpd_statistic(self):
        '''
        makedict_cpd_foldchange is legacy name, 
        where fold change can be any statistic users supply.
        This function populates self.dict_cpd_foldchange,
        and exports attributes for network visualization in Cytoscape 3.
        '''
        #self.dict_cpd_foldchange = self.filter_vis_nodesdict( makedict_cpd_foldchange(self) )
        self.dict_cpd_foldchange = makedict_cpd_foldchange(self)
        s = 'CPD_ID\tCPD_NAME\tCPD_STATISTIC\n'
        for k,v in self.dict_cpd_foldchange.items():
            s += '\t'.join([k, self.network.MetabolicModel.dict_cpds_def.get(k, k).split(';')[0],
                           str(v)]) + '\n'
            
        out = open(os.path.join(self.sifdir, 
                                'Node_attributes.txt'), 'w')
        out.write(s)
        out.close()
    
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


    def export_Mnodes_evidence(self):
        '''
        chosen Mnodes only, as the rest can be looked up in tentative annotation file.
        '''
        resultstr = [['Compound_ID', 'name', 'evidence_score', 'MW', 'all_possibles'],
                      ['', 'input_mz', 'match_form', 'mz_difference', ],
                     ]
        for M in self.mwstr2mnodes.values():
            if M.chosen:
                resultstr.append([M.chosen, 
                                  self.network.MetabolicModel.dict_cpds_def.get(M.chosen, '').split(";")[0], 
                                  M.evidence_score, M.mwstr,
                                  "among possibles %s" %M.id])
            else:
                resultstr.append([M.id, M.cpd.name, M.evidence_score, M.mwstr,])
            for x in M.cpd.hitlist:
                resultstr.append([''] + list(x))
                
        outfile = os.path.join(self.csvdir, 
                               'InspectedNodes_ActivityNetwork.tsv')
        writer = csv.writer( open(outfile, 'wb'), delimiter='\t' )
        writer.writerows( resultstr )
        print_and_loginfo("\nInspected network report was written to\n%s" %outfile)


    def export_metabolite_worksheet(self):
        '''
        collection of all significant metabolites from pathway & module analysis;
        output a table to work on testable metabolites.
        '''
        
        resultstr = [['Compound_ID', 'name', 'input_mz', 'rt', 'match_form', 'mz_difference', 'in_pathways', 
                      'evidence_score', 'MW', 'all_possibles',],
                    ]
        for M in self.mwstr2mnodes.values():
            for x in M.cpd.hitlist:
                if M.chosen:
                    resultstr.append([M.chosen, self.network.MetabolicModel.dict_cpds_def.get(M.chosen, '').split(";")[0], 
                                      #] + list(x) + 
                                      ] + [x[0].mz,x[0].retention_time,x[1],x[2]] + 
                                      ['$'.join(self.network.MetabolicModel.cpd2pathways.get(M.chosen, [])), 
                                       M.evidence_score, M.mwstr, "among possibles %s" %M.id])
                else:
                    resultstr.append([M.id, M.cpd.name] + [x[0].mz,x[0].retention_time,x[1],x[2]] +
                                      ['$'.join(self.network.MetabolicModel.cpd2pathways.get(M.id, [])),
                                      M.evidence_score, M.mwstr,])
            
        
        outfile = os.path.join(self.csvdir, "mcg_metabolite_worksheet_") + self.network.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', 'wb')
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_metabolite_worksheet_")
        
        print_and_loginfo("\nWorksheet of top metabolites was written to\n%s (.tsv and .xlsx)" %outfile)


    def collect_web_export_graphs_test(self):
        '''
        Add edges to activity network from top modules.
        Return expanded activity network and up to top 5 modules.
        '''
        edges_expand = []
        for n in self.top_modules[:5]: edges_expand += n.graph.edges()
        
        expanded_activity_network = self.Inspected
        for e in edges_expand:
            if e not in expanded_activity_network.edges():
                if e[0] in expanded_activity_network.nodes() or e[1] in expanded_activity_network.nodes():
                    expanded_activity_network.add_edge(*e)
        
        self.web_export_networks = [expanded_activity_network] + [M.graph for M in self.top_modules[:5]]


    def collect_web_export_graphs(self):
        '''
        Return activity network and up to top 5 modules.
        Moving activity network after modules
        '''
        self.web_export_networks = [M.graph for M in self.top_modules[:5]] + [self.Inspected]
        vis_nodes = []
        for g in self.web_export_networks: vis_nodes += g.nodes()
        self.vis_nodes = set(vis_nodes)


    def web_export(self):
        '''
        Write HTML and javascript based report.
        The visualization section includes activity network and up to top 5 modules.
        
        
        '''
        self.collect_web_export_graphs()
        
        HTML = HtmlExport()
        title = 'Mummichog Report: ' + self.network.paradict['output']
        HTML.add_element(title, 'h1', '')
        
        HTML.add_element('', 'div', 'colorbar', 'colorbar')

        # place to insert network visusalization
        select_menu = HTML.make_select_menu(len(self.web_export_networks))
        HTML.add_element(select_menu, 'div', 'network_selection', 'network_selection')
        HTML.add_element('', 'div', 'networkvisual', 'networkvisual')
        
        stats = "mummichog version: %s, using metabolic model %s.<br>" %(VERSION, self.network.MetabolicModel.version)
        # stats +=  ' Run parameters - ' + str(self.network.paradict) 
        HTML.add_element(stats, 'div', 'stats')
        
        HTML.add_element('Top pathways', 'h2', '')
        pathwaystablly = self.write_pathway_table()
        HTML.add_element(pathwaystablly, 'div', 'pathwaystablly')
        
        HTML.add_element('Top modules', 'h2', '')
        modulestablly = self.write_module_table()
        HTML.add_element(modulestablly, 'div', 'modulestablly')
        
        HTML.add_element('Top metabolite predictions', 'h2', '')
        metabolitestablly = self.write_metaboite_table()
        HTML.add_element(metabolitestablly, 'div', 'metabolitestablly')
        
        HTML.add_element('More data', 'h2', '')
        moreinformation = '''Full data tables from pathway analysis, module analysis and activity network are stored in the 'tsv/' directory. 
            The tsv files can be imported into a spreadsheet program, e.g. MS Excel.
            Under the 'sif/' directory are files intended for Cytoscape (cytoscape.org) visualization. Please refer to Cytoscape's guides for details.
            Details of this run are recorded in file mummichog.log.
            '''
        HTML.add_element(moreinformation, 'p', 'moreinformation')
        
        footer = '''Mummichog algorithms are described in Li et al. Predicting Network Activity from High Throughput Metabolomics. PLoS Computational Biology (2013); doi:10.1371/journal.pcbi.1003123.\
        This software is provided as is, without warranty of any kind.
        '''
        HTML.add_element(footer, 'footer', '')
        
        # moved to self.export_sif_related()
        # self.export_dict_cpd_statistic()
        
        # push network data into javascript
        HTML.make_js_data(self.web_export_networks, 
                          self.network.MetabolicModel.dict_cpds_def, 
                          self.dict_cpd_foldchange)
        
        outfile = os.path.join(self.rootdir, 'result.html')
        out = open(outfile, 'w')
        out.write(HTML.export_text())
        out.close()
        print_and_loginfo("\nHTML report was written to\n%s" %outfile)
        
        
    def filter_vis_nodesdict(self, d1):
        d2 = {}
        for k in d1.keys():
            if k in self.vis_nodes:
                d2[k] = d1[k]
                
        return d2

    def write_pathway_table(self):
        s = '<table><tr><th>Pathways</th>\
                    <th>overlap_size</th><th>pathway_size</th><th>p-value (raw)</th>\
                    <th>p-value</th></tr>'
        ii = 0
        for P in self.pathway_result:
            ii += 1
            if ii < 6 or P.adjusted_p < SIGNIFICANCE_CUTOFF:
                s += '<tr> <td>' + P.name + '</td> <td>'\
                        + str(P.overlap_size) + '</td><td>' + str(P.cpd_num
                        ) + '</td><td>' + str(round(P.p_FET, 5)) + '</td><td>' \
                        + str(round(P.adjusted_p, 5)) + '</td></tr>'        

        return s + '</table>'

    
    def write_module_table(self):
        s, counter = '', 0
        for M in self.top_modules:
            counter += 1
            nodes = M.graph.nodes()
            names = [self.network.MetabolicModel.dict_cpds_def.get(x, '') for x in nodes]
            s += '<div class="moduleline">' + 'module_' + str(counter) + ", p=" + str(round(M.p_value, 5)) + ", " + str(len(nodes)) + " metabolites" + '</div>'
            s += '<div class="metabolites">' + ', '.join(names) + '</div>'
            
        return s + '\n'



    def write_metaboite_table(self):
        '''
        if M.chosen, metabolite ID is from multiple possibility and M.id can be too long.
        M.id.replace(',', ', ') to give formatting flexibility.
        '''
        s = '<table><tr><th>Compound_ID</th>\
                    <th>name(input_mz)</th><th>evidence_score(match_form)</th><th>MW(mz_diff)</th>\
                    <th>all_possibles</th></tr>'
        
        for M in self.mwstr2mnodes.values():
            if M.chosen:
                s += '<tr> <td>' + M.chosen + '</td> <td>'\
                            + self.network.MetabolicModel.dict_cpds_def.get(M.chosen, ''
                            ).split(";")[0] + '</td><td>' + str(M.evidence_score
                            ) + '</td><td>' + M.mwstr + '</td><td>' \
                            + M.id.replace(',', ', ') + '</td></tr>'
            else:
                s += '<tr> <td>' + M.id + '</td> <td>' + M.cpd.name + '</td><td>' + str(M.evidence_score
                            ) + '</td><td>' + M.mwstr + '</td><td> </td></tr>'
                
            for x in M.cpd.hitlist:
                s += '<tr> <td> </td> <td>%f</td><td>%s</td><td>%f</td><td> </td></tr>' %(x[0].mz, x[1], x[2])
            
        return s + '</table>'
                



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
    
    version 2
    moving everything into EmpiricalCompound space
    
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
        
        self.ListOfEmpiricalCompounds = mixedNetwork.ListOfEmpiricalCompounds
        self.total_number_EmpiricalCompounds = len(self.ListOfEmpiricalCompounds)
        
        
        print_and_loginfo("\nPathway Analysis...")
        
        
    def get_pathways(self, pathways):
        '''
        convert pathways in JSON formats (import from .py) to list of Pathway class.
        currency metabolites removed.
        
        Add list of EmpiricalCompounds
        
        '''
        new = []
        for j in pathways:
            P = metabolicPathway()
            P.json_import(j)
            P.adjusted_p = ''
            
            P.EmpiricalCompounds = self.cpd2empiricalCompounds(P.cpds)
            
            new.append(P)
        return new
        

    def do_permutations(self, pathways, num_perm):
        '''
        Modified from Berriz et al 2003 method.
        After collecting p-values from resampling, do a Gamma fit.
        
        Permutation is simplified in version 2; no more new TableFeatures instances.
        
        '''
        self.permutation_record = []
        print_and_loginfo("Resampling, %d permutations to estimate background ..." 
                          %num_perm)
        
        # this is feature number not cpd number
        N = len(self.mixedNetwork.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii + 1))
            sys.stdout.flush()
            query_EmpiricalCompounds = self.mixedNetwork.compile_significant_EmpiricalCompounds(
                                        random.sample(self.mixedNetwork.mzrows, N) )
            
            self.permutation_record += (self.calculate_permutation_value(set(query_EmpiricalCompounds), pathways))

        # now do Gamma fit
        print_and_loginfo("\nPathway background is estimated on %d random pathway values" 
                          %len(self.permutation_record))
        self.gamma = stats.gamma.fit(1-np.array(self.permutation_record))


    def calculate_permutation_value(self, query_EmpiricalCompounds, pathways):
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


    def adjust_p_by_permutations(self, pathways):
        '''
        EASE score is used as a basis for adjusted p-values,
        as mummichog encourages bias towards more hits/pathway.
        pathways were already updated by first round of Fisher exact test,
        to avoid redundant calculations
        '''
        self.do_permutations(pathways, self.paradict['permutation'])
        
        if self.paradict['modeling'] == 'gamma':
            for P in pathways: 
                P.adjusted_p = self.calculate_gamma_p(P.p_EASE)
        else:
            for P in pathways: P.adjusted_p = self.calculate_p(P.p_EASE, self.permutation_record)
        return pathways
        

    def calculate_p(self, x, record):
        '''
        calculate p-value based on the rank in record of permutation p-values
        '''
        total_scores = [x] + record
        total_scores.sort()
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
    
    def calculate_min_p(self, x, record):
        '''
        calculate p-value based on the rank in record of p-values
        '''
        total_scores = [x] + record
        total_scores.sort()
        D = float(len(total_scores))
        
        return (total_scores.index(x) + 1)/D
    
    def calculate_gamma_p(self, x):
        a, loc, scale = self.gamma
        return 1 - stats.gamma.cdf(1-x, a, loc, scale)
    
    
    def cpd2empiricalCompounds(self, cpds):
        '''
        Mapping cpds to empirical_cpds.
        '''
        cpds_empirical = []
        for c in cpds: cpds_empirical += self.mixedNetwork.Compounds_to_EmpiricalCompounds.get(c, [])
        return set(cpds_empirical)
        
    
    def cpd_enrich_test(self):
        '''
        Fisher Exact Test in cpd space, after correction of detected cpds.
        Fisher exact test is using scipy.stats.fisher_exact
        for right-tail p-value:
        >>> stats.fisher_exact([[12, 5], [29, 2]], 'greater')[1]
        0.99452520602188932
        
        query size is now counted by EmpiricalCompounds.
        
        adjusted_p should be model p-value, not fdr.
        '''
        FET_tested_pathways = []
        qset = self.mixedNetwork.significant_EmpiricalCompounds
        query_set_size = len(qset)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        print_and_loginfo("Query number of significant compounds = %d compounds" %query_set_size)
        
        for P in self.pathways:
            # use the measured pathway size
            P.overlap_EmpiricalCompounds = P.overlap_features = qset.intersection(P.EmpiricalCompounds)
            P.overlap_size = overlap_size = len(P.overlap_EmpiricalCompounds)
            ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                # Fisher's exact test
                P.p_FET = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                   [ecpd_num - overlap_size, negneg]], 'greater')[1]
                # EASE score as in Hosack et al 2003
                P.p_EASE = stats.fisher_exact([[max(0, overlap_size - 1), query_set_size - overlap_size],
                                   [ecpd_num - overlap_size + 1, negneg]], 'greater')[1]
            else:
                P.p_FET = P.p_EASE = 1
                
            FET_tested_pathways.append(P)
            #  (enrich_pvalue, overlap_size, overlap_features, P) 
            
        result = [(P.adjusted_p, P) for P in 
                                        self.adjust_p_by_permutations(FET_tested_pathways)]
        result.sort()
        result = [x[1] for x in result]
        self.record_sigpath_cpds(result)
        
        return result
    
    
    def record_sigpath_cpds(self, result):
        use_result = [r for r in result if r.adjusted_p < SIGNIFICANCE_CUTOFF]
        for P in use_result: self.mixedNetwork.hit_EmpiricalCompounds += P.overlap_EmpiricalCompounds
        print self.mixedNetwork.hit_EmpiricalCompounds[:5]
    


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
    Just need tracking mapping btw compound and EmpiricalCompounds
    
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
        
        significant_EmpiricalCompounds = self.mixedNetwork.significant_EmpiricalCompounds
        cpdlist=[]
        for E in significant_EmpiricalCompounds: cpdlist += E.compounds
        self.seed_cpds = set(cpdlist)
        
        
        
        self.modules_from_input = []
        self.permuation_mscores = []
        self.top_modules = []


    def dispatch(self):
        '''
        Only real modules are saved in total. 
        Permutated modules are not saved but their scores are recorded. 
        '''
        s = "\nModular Analysis, using %d permutations ..." %self.paradict['permutation']
        print s
        logging.info(s)
        self.run_analysis_real()
        self.do_permutations(self.paradict['permutation'])
        self.rank_significance()

    def run_analysis_real(self):
        self.modules_from_input = self.find_modules(self.seed_cpds)

    def run_analysis(self, input_mzlist):
        '''
        This function is used by permutation only as the real input is analyzed by
        self.run_analysis_real(), which stores more information.
        '''
        mscores = [x.A for x in self.mzFeatures2modules(input_mzlist)] or [0]
        
        self.permuation_mscores += mscores


    def mzFeatures2modules(self, mzFlist):
        '''
        This can use mzFeature to Compound as a shortcut for fast permuations tests.
        '''

        cpdlist=[]
        for m in mzFlist:
            cpdlist += self.mixedNetwork.rowindex_to_Compounds.get(m, [])

        return self.find_modules(set(cpdlist))
        

    def mz2modules(self, input_mzlist):
        '''
        wrapper function, mz2cpds then find_modules.
        The permutated input list is fuzzy matched via self.network.mz2cpds().
        
        
        not used now.
        '''
        TFX = TableFeatures(self.network, input_mzlist)
        cpdlist=[]
        for c in TFX.input_cpdlist:
            cpdlist.append(TFX.network.cpd_dict[c].cpd.id)

        return self.find_modules(set(cpdlist))#TFX.input_cpdlist)





    def find_modules(self, cpds):
        '''
        get connected nodes in up to 4 steps.
        modules are set of connected subgraphs plus split moduels within.
        A shaving step is applied to remove excessive nodes that do not 
        connect seeds (thus Mmodule initiation may reduce graph size). 
        A module is only counted if it contains more than one seeds.
        '''
        global SEARCH_STEPS
        seeds, modules, modules2, module_nodes_list = cpds, [], [], []

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
                if 3 < sub.number_of_nodes() < 500:
                    M = Mmodule(self.network, sub, cpds, self.mixedNetwork
                    
                    
                    
                    )
                    modules.append(M)
                
        # add modules split from modules
        if USE_DEBUG:
            logging.info( '# initialized network size = %d' %len(seeds) )
            # need export modules for comparison to heinz
            self.export_debug_modules( modules )
            
        for sub in modules:
            if sub.graph.number_of_nodes() > 5:
                modules2 += [Mmodule(self.network, x, cpds, self.mixedNetwork
                
                
                                )
                             for x in self.split_modules(sub.graph)]
        
        new = []
        for M in modules + modules2:
            if M.graph.number_of_nodes() > 3 and M.nodestr not in module_nodes_list:
                new.append(M)
                module_nodes_list.append(M.nodestr)
                if USE_DEBUG: logging.info( str(M.graph.number_of_nodes()) + ', ' + str(M.A) )

        return new




    def find_modules_v2testing(self, cpds):
        '''
        get connected nodes in up to 4 steps.
        modules are set of connected subgraphs plus split moduels within.
        A shaving step is applied to remove excessive nodes that do not 
        connect seeds (thus Mmodule initiation may reduce graph size). 
        A module is only counted if it contains more than one seeds.
        '''
        global SEARCH_STEPS
        seeds, modules, modules2, module_nodes_list = cpds, [], [], []
        ecpds=[]
        for c in cpds:
            #for ec in self.tf.network.theoretical2empirical[c]:
            #ecpds.append(ec)
            ecpds.append(self.tf.network.theoretical2empirical[c][0])
        for ii in range(SEARCH_STEPS):
            edges = nx.edges(self.network.network, seeds)
            if ii == 0:
                # step 0, counting edges connecting seeds
                edges = [x for x in edges if x[0] in seeds and x[1] in seeds]
                new_network = nx.from_edgelist(edges)
                #print(new_network.number_of_nodes())
            else:
                # step 1, 2, 3, ... growing to include extra steps/connections
                new_network = nx.from_edgelist(edges)
                seeds = new_network.nodes()
                #print(new_network.number_of_nodes())
            for sub in nx.connected_component_subgraphs(new_network):
                if 3 < sub.number_of_nodes() < 500:
                    M = Mmodule(self.network, sub, ecpds, self.tf)
                    if M.graph.number_of_nodes() > 0:
                        modules.append(M)
                
        # add modules split from modules
        if USE_DEBUG:
            logging.info( '# initialized network size = %d' %len(seeds) )
            # need export modules for comparison to heinz
            self.export_debug_modules( modules )
        for sub in modules:
            if sub.graph.number_of_nodes() > 5:
                modules2 += [Mmodule(self.network, x, ecpds, self.tf)
                             for x in self.split_modules(sub.graph)]
        
        new = []
        for M in modules + modules2:
            if M.graph.number_of_nodes() > 3 and M.nodestr not in module_nodes_list:
                new.append(M)
                module_nodes_list.append(M.nodestr)
                if USE_DEBUG: logging.info( str(M.graph.number_of_nodes()) + ', ' + str(M.A) )
        return new

    def export_debug_modules(self, modules):
        '''
        write out initial modules, to be split by alternative algorithm
        '''
        s = ''
        for M in modules: s += M.make_sif_str()
        out = open( os.path.join(self.modules_dir, 'debug_modules.sif'), 'a' )
        out.write(s + '#\n')
        out.close()
        

    def split_modules(self, g):
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
    def split_modules_nemo(self, g):
        '''
        Alternative function using NeMo algorithm for module finding.
        Not used for now.
        '''
        net = NGM.nemo_network(g)
        return [nx.subgraph(g, x) for x in net.find_modules() if len(x) > 3]


    def do_permutations(self, num_perm):
        '''
        Run num_perm permutations on mzlist;
        populate activity scores in self.permuation_mscores
        
        version 2:
        self.ref_featurelist is now using row_numbers
        
        '''
        N = len(self.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii+1))
            sys.stdout.flush()
            self.run_analysis(random.sample(self.ref_featurelist, N))

    def rank_significance(self):
        '''
        compute p-values of modules:
        scores of random modules are fitted to a Gamma distribution,
        p-value is calculated from CDF.
        '''
        print_and_loginfo("\nNull distribution is estimated on %d random modules" 
                          %len(self.permuation_mscores))
        print_and_loginfo("User data yield %d network modules" 
                          %len(self.modules_from_input))
        
        if self.paradict['modeling'] == 'gamma':
            a, loc, scale = stats.gamma.fit(self.permuation_mscores)
            if USE_DEBUG:
                logging.info( 'Gamma fit parameters a, loc, scale = ' + ', '.join([str(x) for x in [a, loc, scale]]) )
            
            for M in self.modules_from_input:
                M.p_value = 1 - stats.gamma.cdf(M.A, a, loc, scale)
            
        else:
            for M in self.modules_from_input:
                M.p_value = self.calculate_p(M.A, self.permuation_mscores)
            
        chosen = [(M.p_value, M) for M in self.modules_from_input]
        chosen.sort()
        self.top_modules = [x[1] for x in chosen if x[0] < SIGNIFICANCE_CUTOFF]
        
        self.record_sigmodule_cpds(self.top_modules)
        
        #if len(self.top_modules) < 10: self.top_modules = [x[1] for x in chosen[:10]]

    def calculate_p(self, x, record):
        '''
        calculate p-value based on the rank in record of scores
        '''
        total_scores = [x] + record
        total_scores.sort(reverse=True)
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D

    def record_sigmodule_cpds(self, top_modules):
        for M in top_modules: self.network.significant_cpdlist += M.graph.nodes()
        


# --------------------------------------------------------
#
# activity network analysis
#

class NetworkInspector:
    '''
    Tally the m/z-cpds responsible for significant pathways and modules;
    inspect their quality in spectral data;
    build a network of high confidence cpds. 
    
    1) weight on evidence from m/z patterns.
    2) resolve conflicting cpds to same m/zs 
    
    There are two layers of resolving conflicting cpds to m/z:
    a) same m/z mapped to cpds with different theoretical MW. 
       Higher M.evidence_score wins.
    b) Cpds with identical theoretical MW will have identical matches.
       Cpd with higher network degree wins.
    
    Draw inspected network.
    Export m/z - cpd mapping data, plus evidence deduction
    
    
    Activity network needs to be redone in version 2.
    Evidence scores should be moved upfront, to initial mapping btw features and metabolites
    
    
    
    
    
    '''
    def __init__(self, TF):
        '''
        class needs a reference metabolic network: hsanet.
        paradict passes user supplied arguments/parameters.
        TF is an instance of TableFeatures, processed input data.
        '''
        self.tf = TF
        self.network = TF.network
        self.paradict = TF.network.paradict


    def inspect_network(self):
        '''
        Inspect cpds responsible for significant pathways and modules.
        self.collect_Mnodes() resolves same m/z matched to multiple cpds; 
        Here, cpds with identical MW are resolved by network degrees. 
        Draw and export inspected network.
        '''
        self.cpds = set(self.tf.input_cpdlist).intersection(
                                        set(self.network.significant_cpdlist))

        self.collect_Mnodes()
        start_nodes = [M.id for M in self.mwstr2mnodes.values() if ',' not in M.id]
        for M in self.mwstr2mnodes.values():
            if ',' in M.id:
                tmp = []
                for c in M.id.split(','):
                    tmp.append((nx.subgraph(self.network.network, [c] + start_nodes
                                           ).number_of_edges(), c))
                tmp.sort()
                start_nodes += [tmp[-1][1]]
                logging.info('    - chosen %s over %s' %(tmp[-1][1], M.id))
                M.chosen = tmp[-1][1]
        
        print_and_loginfo("\nGot ActivityNetwork of %d metabolites." %len(start_nodes))
        Inspected = nx.subgraph(self.network.network, [self.network.cpd_dict[x].cpd.id for x in start_nodes]) #x.split('_')[0]
        return Inspected



    def collect_Mnodes(self):
        '''
        run cpds through Mnode,
        1) filter Mnode by M.evidence_score > threshold (and primary ion if enforced)
        2) merge Mnodes with identical theoretical MW
        
        
        Mnode replaced by EmpiricalCompound
        
        
        
        '''
        mwstr2mnodes = {}
        #
        # to fix
        # EmpiricalCompound is initiated right at reading input data, not be initiatd here
        #
        all_mnodes = [] # [EmpiricalCompound(self.network.cpd_dict[c].id,self.network.cpd_dict[c].mw) for c in self.cpds]

        for ec in self.network.ecs:
            if ec.id in self.cpds:
                all_mnodes.append(ec)

        for M in all_mnodes:
            if M.evidence_score >= self.paradict['evidence']:
                if not self.paradict['force_primary_ion'] or (
                            self.paradict['force_primary_ion'] and M.primary_ion_present):
                    if not mwstr2mnodes.has_key(M.mwstr): mwstr2mnodes[M.mwstr] = M
                    else: mwstr2mnodes[M.mwstr].id += ',' + M.id
        
        self.mwstr2mnodes = mwstr2mnodes



