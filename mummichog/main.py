#!/usr/bin/env python
# Copyright (c) 2010-2020 Shuzhao Li.
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
mummichog 2,
pathway and network analysis for metabolomics

@author: Shuzhao Li

Online documentation: http://mummichog.org
Major change in data structure from version 1. 
Only the default metabolic model (human_model_mfn) and a worm model are ported to version 2 so far.
This is intended to be the branch of command line version. 
A server version is upcoming, where metabolic models are stored a separate online database.
'''


from .functional_analysis import *
from .reporting import *


def main():
    
    print (fishlogo)
    print ( "mummichog version %s \n" %VERSION )
    optdict = dispatcher()

    print_and_loginfo("Started @ %s\n" %time.asctime())
    userData = InputUserData(optdict)
    
    #specify which metabolic model 
    if userData.paradict['network'] in ['human', 'hsa', 'Human', 'human_mfn', 'hsa_mfn', '']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'human_model_mfn' ])
    elif userData.paradict['network'] in ['worm', 'C. elegans', 'icel1273', 'Caenorhabditis elegans']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'worm_model_icel1273' ])
        
    else:
        raise KeyError( "Unsupported species/model. Pls contact author." )
    
    mixedNetwork = DataMeetModel(theoreticalModel, userData)

    # getting a list of Pathway instances, with p-values, in PA.resultListOfPathways
    PA = PathwayAnalysis(mixedNetwork.model.metabolic_pathways, mixedNetwork)
    PA.cpd_enrich_test()
    
    # Module analysis, getting a list of Mmodule instances
    MA = ModularAnalysis(mixedNetwork)
    MA.dispatch()
    
    # do activity network
    AN = ActivityNetwork( mixedNetwork, set(PA.collect_hit_Trios() + MA.collect_hit_Trios()) )
    
    Local = LocalExporting(mixedNetwork, PA, MA, AN)
    Local.run()
    
    Web = WebReporting(Local, PA, MA, AN)
    Web.run()
    
    print_and_loginfo("\nFinished @ %s\n" %time.asctime())




#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    main()


