% Smolyak_Polynomial.m is a routine that constructs the multidimensional 
% basis functions of Smolyak polynomial of the approximation level that  
% corresponds to the previously constructed Smolyak (multidimensional) grid 
% points; see "Smolyak method for solving dynamic economic models: Lagrange 
% interpolation, anisotropic grid and adaptive domain" by Kenneth L. Judd, 
% Lilia Maliar, Serguei Maliar and Rafael, (2014). Journal of Economic 
% Dynamics and Control 44, 92–123 (henceforth, JMMV (2014)). 
%
% This version: Novenber 5, 2014. First version: May 30, 2011.
% -------------------------------------------------------------------------
% Inputs:  "points"    is the matrix of points in which the polynomial basis 
%                      functions must be evaluated; numb_pts-by-d
%          "d"         is the number of dimensions (state variables)
%          "Smol_elem" is the matrix of the subindices of the Smolyak
%                      unidimensional elements; these elements can be either 
%                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
%                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
%                      in the former case, the indices i1,...,id that jointly 
%                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
%                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
%                      in the later case, they are a subset of the above 
%                      indices 
%
% Output:  "Smol_bases" is the matrix of multidimensional basis functions of 
%                      Smolyak polynomial of the given level of approximation, 
%                      evaluated in data matrix "points"
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


 function Smol_bases = Smolyak_Polynomial(points,d,mu,Smol_elem)
 
% Smolyak polynomial is given by the sum of multidimensional basis functions, 
% multiplied by the coefficients; see formula (15) in JMMV (2014). By 
% convention, the first basis function is given by 1 (unity). 

% Unidimensional basis functions are given by Chebyshev polynomial bases; 
% in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
% degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
% this notation here

% 1. Construct the unidimensional basis functions and evaluate them in
% all the points of matrix "points"
% -------------------------------------------------------------------------
i_max = mu+1;                   % The maximum subindex of unidimensional
                                % set S_i whose points are used to
                                % construct Smolyak grid; e.g., for mu=1, 
                                % we consider up to S_i_max={0,-1,1} where
                                % i_max=mu+1=1+1=2
                                
% Compute the number of elements in the i_max-th unidimensional set of 
% elements S_i_max; this coincides with a maximum subindex of elements 
% (unidimensional grid point or unidimensional basis function); 
% see Section 2.2.1 in JMMV (2014)

m_i_max(i_max==1) = 1;          % If i_max=1, then m(i_max)=1, i.e., set  
                                % S_1={0} and the maximum subindex is 1
m_i_max(i_max>1) =  2.^(i_max(i_max>1)-1) + 1;
                                % If i_max>1, then m(i_max)= 2^(i_max-1)+1;
                                % e.g., for S_2={0,-1,1}, the maximum
                                % subindex is 3

                                 
numb_pts = size(points,1);      % Compute the number of points (rows),   
                                % "numb_pts" in the matrix of points 
                                % "points", in which the polynomial bases  
                                % must be evaluated              
phi = ones(numb_pts,d,m_i_max); 
                                % Allocate memory to a matrix of the
                                % unidimensional polynomial bases "phi_n",  
                                % evaluated in all the points 
                     
% For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
% our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
% construction
                               
phi(:,:,2) = points;            % For a polynomial bases "phi_n" with n=2, 
                                % we have phi_n(x) is x; evaluating it in 
                                % all the points gives us matrix "points"; 
                                % numb_pts-by-d                    
 for j = 3:m_i_max              % For polynomial bases "phi_n", from n=3 to
                                % n=m_i_max, ...
    phi(:,:,j) = 2.* phi(:,:,2).*phi(:,:,j-1) - phi(:,:,j-2); 
                                % Use the recurrence formula to compute the 
                                % Chebyshev polynomial basis functions of 
                                % the degrees from 2 to m_i_max-1
 end 

 
% 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
% required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
% and 3.4.2 for examples
% ----------------------------------------------------------------------
 Smol_bases = [];              % Initially, the matrix of multidimensional 
                               % polynomial bases is empty
                               
 numb_terms = size(Smol_elem,1);
                               % Compute the number of terms (i.e., multi-
                               % dimensional polynomial bases) in Smolyak 
                               % polynomial                                
 
 for jt = 1:numb_terms         % For each term of Smolyak polynomial, ...
                                
     index_row = Smol_elem(jt,:);
                               % Identify the subindices of the unidimensional
                               % basis function that constitute a jt-th multi-
                               % dimensional basis function; this is a jt-th
                               % row of matrix "Smol_elem"; 1-by-d
     product = ones(numb_pts,1); 
                               % Initialize a vector of products of unidimen-
                               % sional basis functions; numb_pts-by-1
     for jd = 1:d              % For each dimension (state variable), ...
         n = index_row(jd);    % A subindex of a unidimensional basis function 
                               % phi_n in a dimension jd is denoted n
         if n ~= 1;            % If the subindex of unidimensional basis 
                               % function is not equal to unity, ...
         product = product.*phi(:,jd,n);
                               % Compute the product of basis functions
                               
         % Otherwise, i.e., if n = 1, there is no need to compute the 
         % product of unidimensional basis functions, as it's equal to unity 
        end
     end
     Smol_bases = cat(2,Smol_bases,product);
                               % Attach to the previously obtained matrix of 
                               % multidimensional basis functions a new
                               % product of unidimensional basis functions;
                               % e.g., for mu=1 and d=2, basis_bs is of
                               % size numb_pts-by-5
 end
