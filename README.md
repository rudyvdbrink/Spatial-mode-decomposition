# Spatial-mode-decomposition
Example code to compute spatial modes, and run statistical analyses of mode variance. The code is by no means optimized, for the sake of interpretability, but should run relatively quickly nonetheless. I wrote this code in MATLAB 2012a, and it should be forward compatible. 

The code is written so as to produce modes that are more strongly expressed in the atomoxetine condition than in the placebo condition, and compare the variance of the first mode between conditions using cross-validation (both the average and with ROC analysis). 

There are four relevant files:
1) Spatial_mode_decomposition.m:
	This contains the relevant code to compute spatial modes. 
	I've tried to comment most lines to explain what is 	computed and how it's done. When running this, make sure 	that MATLAB has access to all other files.
2) permtest.m:
	Function for permutation testing.
3) sroc.m:
	Function for ROC analysis.
4) M.mat:
	Data file that contains a matrix of z-scored BOLD time-	series to which the mode decomposition is applied. In 	principle you could replace this with your own data, but 	you'd need to make sure that it is of the right size: 	participants by conditions (two max) by time points by 	brain regions. So four dimensions in total. You'd	also 	need to make sure that the binning for ROC analysis is 	done correctly if the number of time points differs from 	my own.   


Ruud van den Brink, 2017

If you use any of this this code for a publication, please cite: 
van den Brink, Nieuwenhuis, & Donner (2018) Amplification and Suppression of Distinct Brain-wide Activity Patterns by Catecholamines 
Preprint: https://www.biorxiv.org/content/early/2018/06/04/270645.full.pdf
