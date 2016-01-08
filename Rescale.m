% Rescale.m is a routine that rescales variables from one range to another; 
% see "Smolyak method for solving dynamic economic models: Lagrange interpo-
% lation, anisotropic grid and adaptive domain" by Kenneth L. Judd, Lilia  
% Maliar, Serguei Maliar and Rafael Valero, (2014), Journal of Economic  
% Dynamics and Control 44, 92–123 (henceforth, JMMV (2014)), Appendix C 
%
% This version: Novenber 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Input:   "orig_matrix" is a matrix whose values we want to rescale to be 
%          in a different inteval
%          "min_resc" and "max_resc" are, respectively, the vectors of the   
%          lower and upper limits of the values of the rescaled matrix (output)
%          "min_orig" and "max_orig" are, respectively, the vectors of the   
%          lower and upper limits of the values of the original matrix 
%          "orig_matrix"
%
% Output:  "resc_matrix" is the resulting rescaled matrix
% 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function resc_matrix = Rescale(orig_matrix,min_resc,max_resc,min_orig,max_orig)

[fil,~] = size(orig_matrix); % Compute the number of columns in "orig_matrix" 

resc_matrix = (ones(fil,1)*((max_resc-min_resc)./(max_orig-min_orig))).*(orig_matrix-ones(fil,1)*min_orig)+ones(fil,1)*min_resc;
% Compute the rescaled matrix of data points; see equation (49) in JMMV (2014)
