Mummichog
=========

Mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics.
It leverages the organization of metabolic networks to predict functional activity directly from feature tables,
bypassing metabolite identification. The features include

* computing significantly enriched metabolic pathways
* identifying significant modules in the metabolic network
* visualization of top networks in web browser
* visualization that also plugs into Cytoscape
* tentative annotations
* metabolic models for different species through plugins

This is mummichog package version 2. Version 3 is under development.

Please note that mummichog-server is a different package/repository.


Installation
------------

*Python 3 is required for Mummichog version 2.3 and beyond.*

*Mummichog version 2.2 was the last version using Python 2; new branch as mummichog-python2*

Mummichog can be installed using pip (pip Installs Packages), the Python package manager. The command below will install the default (version 2):

    pip install mummichog

This is OS independent. To read more on pip `here <https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py>`.

One can also run mummichog without installing it. Direct python call on a downloaded copy can work.

The initial paper on mummichog is described in Li et al. Predicting Network Activity from High Throughput Metabolomics. PLoS Computational Biology (2013); doi:10.1371/journal.pcbi.1003123.. More on `project website <http://mummichog.org>`.
