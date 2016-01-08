% RunMeExample.m a simple example illustrating the Smolyak anisotropic 
% method for interpolation of a 2-dimensional function; see "Smolyak method 
% for solving dynamic economic models: Lagrange interpolation, anisotropic 
% grid and adaptive domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar 
% and Rafael Valero, (2014), Journal of Economic Dynamics and Control 44, 
% 92–123 (henceforth, JMMV (2014)), Section 2.2.3 
%
% This version: November 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Input:   Choose anisotropy by changing "vector_mus_dimensions"
%
% Output:  Graphs that compare the true function and Smolyak interpolation 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Smolyak anisotropic method for 2 dimensions;
% -----------------------------------------------
d     = 2;                  % Number of dimensions 
                            
vector_mus_dimensions = [5 2];
                            % Introduce the level of approximation in every 
                            % dimension from 1 to 10; see Section 4 of 
                            % JMMV (2014)

mu_max  = max(vector_mus_dimensions);
                            % Compute the maximum level of approximation
                            % across all dimensions

Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu_max);
    % Construct the matrix of indices of multidimesional Smolyak elements 
    % (grid points and polynomial basis functions) that satisfy the usual 
    % isotropic Smolyak rule for the approximation level equal to "mu_max"

Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
    % Select from the matrix of isotropic indices "Smol elem_iso" a subset of 
    % indices that correspond to the given anisotropic "vector_mus_dimensions"
   
Smol_grid_ani = Smolyak_Grid(d,mu_max,Smol_elem_ani); 
    % Construct the Smolyak grid for the given subindices of anisotropic  
    % Smolyak elements "Smol_elem_ani"

Smol_polynom = Smolyak_Polynomial(Smol_grid_ani,d,mu_max,Smol_elem_ani); 
    % Matrix of the polynomial basis functions evaluated in the grid
    % points

%scatter(Smol_grid_ani(:,1),Smol_grid_ani(:,2),'filled'),title('Smolyak grid')
    % Enable it to draw the Smolyak grid
    
    
% 2. Interpolation of a function y=2*x1 .*exp(-4*x1.^2-16*x2.^2);
% ---------------------------------------------------------------
y_Smolyak = 2*Smol_grid_ani(:,1) .* exp(-4*Smol_grid_ani(:,1).^2 - 16*Smol_grid_ani(:,2).^2);
    % Evaluate the function on the Smolyak grid
    
b = Smol_polynom\y_Smolyak; 
    % Compute the coefficients of Smolyak interpolating polynomial

    
% 3. Compare the true and interpolated functions on a dense grid
% ---------------------------------------------------------------
[x1,x2] = meshgrid(-1:.05:1, -1:.05:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

y_true = 2*x1.* exp(-4*x1.^2 - 16*x2.^2);
    % Evaluate the true function on the grid    

subplot(1,2,1), surf(x1,x2,y_true),title('True function')
    % Plot the true function

for j=1:size(x1(:,1),1)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

subplot(1,2,2), surf(x1,x2,y_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation




