% Accuracy_Test_Smolyak.m is a routine for evaluating the accuracy of     
% Smolyak solutions of the multi-country model: it computes approximation  
% errors in the optimality conditions on a given set of points in the state 
% space; see "Smolyak method for solving dynamic economic models: Lagrange 
% interpolation, anisotropic grid and adaptive domain" by Kenneth L. Judd, 
% Lilia Maliar, Serguei Maliar and Rafael Valero (2014). Journal of Economic 
% Dynamics and Control. Volume 44, July 2014, Pages 92–123 (henceforth, 
% JMMV, 2014).  
%
% This version: November 5, 2014. First version: December 17, 2012.
% -------------------------------------------------------------------------
% Inputs:    "k" and "a" are, respectively, current-period capital and 
%                        productivity levels, in the given set of points  
%                        on which the accuracy is tested; 
%            "bk"        are the coefficients of the capital policy functions 
%                        of N countries;
%            "IM"        is the integration method for evaluating accuracy, 
%                        IM=1,2,..,10=Gauss-Hermite quadrature rules with  
%                        1,2,...,10 nodes in each dimension, respectively, 
%                        IM=11=Monomial rule with 2N nodes,
%                        IM=12=Monomial rule with 2N^2+1 nodes;
%            "alpha", "gam", "phi", "beta", "A", "tau", "rho" and "vcv"
%                        are the parameters of the model;
%            "discard"   is the number of data points to discard 
%            "min_Cheb" and "max_Cheb" 
%                        are, respectively, the lower and upper 
%                        bounds of Chebyshev polynomials (equal to -1 and 1 
%                        for each dimension
%            "min_Data" and "max_Data" 
%                        are, respectively, the minimum and maximum values 
%                        of the simulated series
%            "mu"        is the level of approximation
%            "Smol_elem" is the subindices of multidimensional Smolyak elements 
%                        (points or polynomial basis functions
%            
%
% Outputs:   "Errors_mean" and "Errors_max" 
%                        are, respectively, the mean and maximum residuals 
%                        across all optimality conditions;
%            "Errors_max_EE", "Errors_max_MUC", "Errors_max_MUL", and 
%            "Errors_max_RC" 
%                        are the maximum approximation errors disaggregated 
%                        by optimality conditions 
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Errors_mean Errors_max time_test] = Accuracy_Test_Smolyak(k,a,bk,IM,alpha,gam,delta,beta,A,tau,rho,vcv,discard,min_Cheb,max_Cheb,min_Data,max_Data,d,mu,Smol_elem)

tic                   % Start counting time needed to run the test        

[P,N] = size(a);      % Infer the number of points P on which the accuracy 
                      % is evaluated and the number of countries N
                      

% 2. Integration method for evaluating accuracy 
% ---------------------------------------------
if (IM>=1)&&(IM<=10)
    [n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(IM,N,vcv);
                             % Compute the number of integration nodes, 
                             % n_nodes, integration nodes, epsi_nodes, and 
                             % integration weights, weight_nodes, for Gauss-
                             % Hermite quadrature integration rule with IM 
                             % nodes in each dimension
elseif IM == 11
    [n_nodes,epsi_nodes,weight_nodes] = Monomials_1(N,vcv);
                             % Monomial integration rule with 2N nodes
elseif IM == 12
    [n_nodes,epsi_nodes,weight_nodes] = Monomials_2(N,vcv);
                             % Monomial integration rule with 2N^2+1 nodes
end    

% 3. Polynomial bases for the test
%--------------------------------- 
X = Smolyak_Polynomial( Rescale([k a],min_Cheb,max_Cheb,min_Data,max_Data),d,mu,Smol_elem);
% Form a matix of Smolyak polynomial basis functions

% 4. Given the solution for capital, compute consumption on the given set 
% of points  
%------------------------------------------------------------------------
for p = 1:P;                 % For each given point, ...     
    
       p                     % Display the point (with the purpose of  
                             % monitoring the progress) 
                                     
        % 4.1 Variables in point p 
        % ------------------------        
        k0 = k(p,1:N);       % N capital stocks of period t
        a0 = a(p,1:N);       % N productivity levels of period t
        X0 = X(p,:);         % Complete (second-degree) polynomial 
                             % bases at t
                                    
        % 4.2 Capital and consumption choices at t
        % ----------------------------------------
        k1 = X0*bk; 
        % Compute a row-vector of capital of period t+1 (chosen at t) using
        % the corresponding capital policy functions; 1-by-N
        
        C0 = (A*k0.^alpha.*a0 - k1+k0*(1-delta))*ones(N,1);
        % C is computed by summing up individual consumption, which in turn, is 
        % found from the individual budget constraints; 1-by-1

        c0 = C0*ones(1,N)/N;  % Individual consumption is the same for all
                              % countries; 1-by-N 

        % 4.3 Capital and consumption choices at t+1
        %-------------------------------------------
        a1 = (ones(n_nodes,1)*a0).^rho.*exp(epsi_nodes);    
        % Compute the next-period productivity levels in each integration node
        % using condition (?) in the online appendix; n_nodes-by-N

        k1_dupl = ones(n_nodes,1)*k1; 
        % Duplicate k1 n_nodes times to create a matrix with n_nodes identical
        % rows; n_nodes-by-N 

        X1 = Smolyak_Polynomial( Rescale([k1_dupl a1],min_Cheb,max_Cheb,min_Data,max_Data),d,mu,Smol_elem);
        % Form a matrix of the polynomial basis functions  
        
        k2 = X1*bk; 
        % Compute capital of period t+2 (chosen at t+1) using the second-
        % degree capital policy functions; n_nodes-by-N 

        C1 = (A*k1_dupl.^alpha.*a1 - k2+k1_dupl*(1-delta))*ones(N,1);
        % Aggregate consumption is computed by summing up individual consumption, 
        % which in turn, is found from the individual budget constraints; 
        % n_nodes-by-1

        c1 = C1*ones(1,N)/N;    % Individual consumption is the same for 
                                % all countries; n_nodes-by-N

% 5. Approximation errors in point p
%-----------------------------------
        
        % 5.1 Lagrange multiplier associated with the aggregate resource
        % constraint
        %---------------------------------------------------------------
        for j = 1:N
            MUC0j(1,j) = tau*c0(1,j).^(-gam); 
            % Compute a country's marginal utility of consumption multiplied 
            % by its welfare weight
        end
        lambda0 = mean(MUC0j,2);
        % An optimality condition w.r.t. consumption of period t equates 
        % the Lagrange multiplier of the aggregate resource constraint of 
        % period t and each country's marginal utility of consumption 
        % multiplied by its welfare weight; to infer the Lagrange multiplier,  
        % we average across N countries; 1-by-1
        
        for j = 1:N
            MUC1j(1:n_nodes,j) = tau*c1(1:n_nodes,j).^(-gam);
            % Compute a country's marginal utility of consumption multiplied 
            % by its welfare weight
        end
        lambda1 = mean(MUC1j,2);
        % Similarly, the Lagrange multiplier of the aggregate resource 
        % constraint of period t+1 is equal to a country's marginal utility 
        % of consumption multiplied by its welfare weight; to infer the 
        % Lagrange multiplier, we average across N countries; 1-by-n_nodes
        
        % 5.2 Unit-free Euler-equation errors
        %------------------------------------
        for j = 1:N
            Errors(p,j) = 1-weight_nodes'*(beta*lambda1/lambda0.*(1-delta+alpha*A*k1(1,j)^(alpha-1)*a1(1:n_nodes,j)));
        % A unit-free Euler-equation approximation error of country j
        end
        
        % 5.2 Unit-free errors in the optimality conditions w.r.t. consumption
        %---------------------------------------------------------------------
        for j = 1:N;
            Errors(p,N+j) = 1-lambda0./(tau*c0(1,j)^(-gam)); 
        % A unit-free approximation error in the optimality condition w.r.t. 
        % consumption of country j (this condition equates marginal utility 
        % of consumption, multiplied by the welfare weight, and the 
        % Lagrange multiplier of the aggregate resource constraint)
        end
        
        % 5.3 Unit-free errors in the optimality conditions w.r.t. labor 
        %---------------------------------------------------------------
        Errors(p,2*N+1:3*N) = zeros(N,1);
        % These errors  are zero by construction 
        
        % 5.4 Unit-free approximation error in the aggregate resource constraint
        %-----------------------------------------------------------------------
        Errors(p,3*N+1) = 1-(c0(1,1:N) + k1(1,1:N)-k0(1,1:N)*(1-delta))*ones(N,1)/((A*k0(1,1:N).^alpha.*a0(1,1:N))*ones(N,1));
        % This error is a unit-free expression of the resource constraint  
        % in the online appendix of JMM (2011)
        
        % 5.5 Approximation errors in the capital-accumulation equation
        %--------------------------------------------------------------
        Errors(p,3*N+2:4*N+1) = zeros(N,1);
        % These errors are always zero by construction
        
        % For this model, GSSA produces zero errors (with machine accuracy)
        % in all the optimality conditions except of the Euler equation. 
        % Here, the errors in all the optimality conditions are introduced 
        % to make our accuracy measures comparable to those in the February 
        % 2011 special issue of the Journal of Economic Dynamics and Control 
end

% 6. Mean and maximum approximation errors computed after discarding the  
% first "discard" observations
%-----------------------------------------------------------------------

        % 6.1 Approximation errors across all the optimality conditions
        %--------------------------------------------------------------
        Errors_mean = log10(mean(mean(abs(Errors(1+discard:end,:))))); 
        % Average absolute approximation error 
        
        Errors_max = log10(max(max(abs(Errors(1+discard:end,:)))));    
        % Maximum absolute approximation error
 
        % 6.2 Maximum approximation errors disaggregated by the optimality 
        % conditions
        %-----------------------------------------------------------------
%         Errors_max_EE = log10(max(max(abs(Errors(1+discard:end,1:N)))));    
%         % Across N Euler equations
% 
%         Errors_max_MUC = log10(max(max(abs(Errors(1+discard:end,N+1:2*N)))));    
%         % Across N optimality conditions w.r.t. consumption (conditions on   
%         % marginal utility of consumption, MUC)
% 
%         Errors_max_MUL = log10(max(max(abs(Errors(1+discard:end,2*N+1:3*N)))));    
%         % Across N optimality conditions w.r.t. labor (conditions on marginal 
%         % utility of labor, MUL)
% 
%         Errors_max_RC = log10(max(max(abs(Errors(1+discard:end,3*N+1)))));    
%         % In the aggregate resource constraint 
 
% 7. Time needed to run the test
%-------------------------------
time_test = toc;         