#!/usr/bin/env python
# Copyright (c) 2010-2018 Shuzhao Li.
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
pathway and network analysis for untargeted and targeted metabolomics

@author: Shuzhao Li

online documentation:
http://mummichog.org



'''


from functional_analysis import *
from reporting import *


def main():
    
    print fishlogo
    print ( "mummichog version %s \n" %VERSION )
    optdict = dispatcher()

    print_and_loginfo("Started @ %s\n" %time.asctime())
    userData = InputUserData(optdict)
    
    # can specify which model in metabolicModels[]
    
    theoreticalModel = metabolicNetwork(metabolicModels[ 'human_model_mfn' ])
    
    #theoreticalModel = metabolicNetwork(metabolicModels[ 'recon2' ])
    
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


