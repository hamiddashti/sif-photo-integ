author: C. van der Tol
date: 10 July 2018

With this code, one can reproduce the analysis and the figures in the paper of Van der Tol et al., submitted to RSE on 10 July 2018, entitled: The scattering and re-absorption of red and near-infrared chlorophyll fluorescence in the models Fluspect and SCOPE

The following codes are masters:
do_analysis		this code carries out the retrieval of leaf properties and the fluorescence emission spectrum from measurements	
makeplots_revision1	this code carried out forward simulations, RMSE's etc, and produces the graph

in do_analysis, set some options to make it run the necessary parts. All parts have been run in preparation of the paper.

SCOPE simulations for the TOC sensitivity analysis have been carried out with SCOPE version 1.70, see github.com/christiaanvandertol
The output is located in ../data/output/, along with the parameter settings.

Leaf level retrievals are stored in ../data/output/fluspect_output/2015/. Look at do_analysis to see how they have been produced from the raw measurements.

The field measurements are located in ../data/measured.
The file ../data/optical_coefficients provides an input file needed for Fluspect.

