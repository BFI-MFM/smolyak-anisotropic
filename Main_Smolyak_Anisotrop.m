% MATLAB software that solves a multi-country model using the anisotropic
% version of the Smolyak method, as described in the article  "Smolyak 
% method for solving dynamic economic models: Lagrange interpolation, 
% anisotropic grid and adaptive domain" by Kenneth L. Judd, Lilia Maliar, 
% Serguei Maliar and Rafael Valero (2014), Journal of Economic Dynamics 
% and Control, 44, 92–123 (henceforth, JMMV, 2014). 
% 
% This software is based on that of Lilia  Maliar and Serguei Maliar for 
% solving the multi-country model using the generalized stochastic simulation 
% algorithm (GSSA) method, as described in the paper "Numerically Stable and 
% Accurate Stochastic Simulation Approaches for Solving Dynamic Economic 
% Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011), 
% Quantitative Economics 2/2, 173–210 (henceforth, JMM, 2011). The modifi-
% cations made are concerned with the construction of the grid and polynomials. 
% 
% This version: November 24, 2014. First version: December 17, 2012.
% 
% ------------------------------------------------------------------------
% The software uses the following files: 
% ------------------------------------------------------------------------
% 0. "RunMeExample.m"            is a simple example illustrating the Smolyak 
%                                anisotropic method for interpolation of a
%                                2-dimensional function
% 1. "Main_Smolyak_Anisotrop.m"  is a main file for computing a solution to 
%                                the multi-country model using the anisotropic 
%                                version of the Smolyak method; the standard 
%                                isotropic Smolyak method is a special case
%                                of the described anisotropic method
% 2. "Accuracy_Test_Smolyak.m"   computes residuals of the optimality 
%                                conditions of the multi-country model on a   
%                                given set of points in the state space  
% 3. "Smolyak_Elem_Isotrop.m"    constructs the subindices of Smolyak elements 
%                                (grid points and basis functions) for the 
%                                isotropic case; 
% 4. "Smolyak_Elem_Anisotrop.m"  selects a subset of the subindices of Smolyak 
%                                elements that correspond to the given 
%                                anisotropic case from a set of subindices 
%                                of the Smolyak isotropic elements                                 
% 5. "Smolyak_Grid"              constructs a multidimensional Smolyak grid, 
%                                given a set of the indices (either isotropic
%                                or anisotropic) of the Smolyak elements
% 6. "Smolyak_Polynomial.m"      constructs the Smolyak basis functions that
%                                constitute a Smolyak polynomial function
%                                given a set of the indices (either isotropic
%                                or anisotropic) of the Smolyak elements
% 7. "Rescale"                   rescales variables from one range to another 
% 8. "Productivity.m"            generates random draws of the productivity  
%                                shocks and simulates the corresponding series  
%                                of the productivity levels; borrowed from JMM 
%                                (2011)
% 9. "Monomials_1.m"             constructs integration nodes and weights for 
%                                an N-dimensional monomial (non-product) 
%                                integration rule with 2N nodes; borrowed from 
%                                JMM (2011) 
% 10. "Monomials_2.m"            constructs integration nodes and weights for 
%                                an N-dimensional monomial (non-product) 
%                                integration rule with 2N^2+1 nodes; borrowed 
%                                from JMM (2011)
% 10. "GH_Quadrature.m"          constructs integration nodes and weights for  
%                                the Gauss-Hermite rules with the number of  
%                                nodes in each dimension ranging from one to  
%                                ten; borrowed from JMM (2011)                     
% 11. "aT20200N10.mat"           contains the series of the productivity  
%                                levels of length 20,200 for 10 countries that 
%                                are used for computing solutions and for 
%                                evaluating accuracy; borrowed from JMM (2011)
% -------------------------------------------------------------------------
% Copyright © 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
% rights reserved. The code may be used, modified and redistributed under  
% the terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;
 
% 1. Choose the number of countries
% ---------------------------------
N     = 3;       % Choose the number of countries, 1<=N<=10 (note that the 
                 % code also works for the one-country case, N=1)
d     = 2*N;     % Number of dimensions (state variables); per each country, 
                 % the state variables are capital and productivity                 

% 2. Model's parameters
% ---------------------
gam     = 1;      % Utility-function parameter
alpha   = 0.36;   % Capital share in output
beta    = 0.99;   % Discount factor
delta   = 0.025;  % Depreciation rate 
rho     = 0.95;   % Persistence of the log of the productivity level
sigma   = 0.01;   % Standard deviation of shocks to the log of the 
                  % productivity level
vcv = sigma^2*(eye(N)+ones(N,N)); 
                  % Variance-covariance matrix of the countries' productivity 
                  % shocks in which diagonal terms are equal to 2*sigma^2   
                  % and in which off-diagonal terms are equal to sigma^2; 
                  % this vcv follows from the assumption that a country's 
                  % shock has both common-for-all-countries and country-
                  % specific components; N-by-N

                                     
% 3. The normalizing constant, A, and welfare weight, tau
% -------------------------------------------------------
A       = (1-beta+beta*delta)/alpha/beta;  % The normalizing constant in output  
tau     = 1;                               % The welfare weight of country j 

% The above normalization ensures that steady state of capital of all 
% countries is equal to one 



% _________________________________________________________________________
%                               
% Compute an initial guess for capital policy functions using GSSA. Namely, 
% compute a first-degree polynomial solution using a one-node Monte Carlo  
% integration method (this solution will be used as an initial guess for
% the Smolyak method); borrowed from JMM (2011) 
% _________________________________________________________________________
%
tic;                  % Start counting time needed to compute the initial 
                      % guess

% 4. Initial condition for simulation
% -----------------------------------

k(1,1:N) = 1;  % Initial condition for capital (is equal to steady state)
a(1,1:N) = 1;  % Initial condition for the productivity level (is equal to 
               % steady state)

% 5. Construct the productivity levels, a 
% ---------------------------------------
T     = 10000;   % Choose the simulation length for the solution procedure,
                 % T<=10,000   
                 
% To solve models with N>10 or T>10,000, one needs to simulate new series
% of the productivity levels by enabling the code in paragraph 5.  

% a20200 = Productivity(T,N,a(1,1:N),sigma,rho);
                               % Generate a random draw of the productivity 
                               % shocks and simulate the corresponding  
                               % series of the productivity levels of length   
                               % T periods for N countries 
% save aT20200N10 a20200;       % Save the series of the productivity levels  
                               % into a file "aT20200N10.mat" 
load aT20200N10;               % Load the previously saved series of the 
                               % productivity levels of length 20,200 for 
                               % 10 countries (the first 10,000 observations
                               % are used for finding a solution, and the 
                               % remaining 10,200 observations are used for
                               % evaluating accuracy)
a = a20200(1:T,1:N);           % Restrict the series of the productivity 
                               % levels for the solution procedure to the 
                               % given T<=10,000 and N<=10
                            
% 6. The GSSA parameters  
% ----------------------
kdamp     = 0.1;     % Damping parameter for (fixed-point) iteration on 
                     % the coefficients of the capital policy functions
dif_GSSA_1d = 1e+10; % Set the initial difference between the series from
                     % two iterations in the convergence criterion (condition
                     % (10) in JMM, 2011) to a very large number
% To achieve convergence under N>10, one may need to modify the values of 
% the damping parameter kdamp or refine the initial guess 

% 7. Initialize the first-degree capital policy functions of N countries 
%-----------------------------------------------------------------------                          
bk_1d  = [zeros(1,N); diag(ones(1,N)*0.9);diag(ones(1,N)*0.1)]; 
% Matrix of polynomial coefficients of size (1+2N)-by-N: for each country  
% (each column), 1+2N rows correspond to a constant, N coefficients on the
% capital stocks, k(t,1),...,k(t,N), and N coefficients on the productivity 
% levels, a(t,1),...,a(t,N)

% As an initial guess, assume that a country's j capital depends only on 
% its own capital and productivity level as k(t+1,j)=0.9*k(t,j)+0.1*a(t,j); 
% (note that in the steady state, we have k(t+1,j)=0.9*k(t,j)+0.1*a(t,j)=1)

% Note that diag(ones(1,N)*q) delivers an N-by-N matrix with  diagonal 
% entries equal to q. 

% 8. Initialize the capital series
% --------------------------------
k_old = ones(T+1,N);   % Initialize the series of next-period capital of N
                       % countries; these series are used to check the
                       % convergence on the subsequent iteration (initially, 
                       % capital can take any value); (T+1)-by-N
                       
% 9. The main iterative cycle of GSSA
% -----------------------------------              
while dif_GSSA_1d > 1e-4*kdamp;   % 10^4*kdamp is a convergence parameter,
                                  % adjusted to the damping parameter; see 
                                  % JMM (2011) for a discussion

% 9.1 Generate time series of capital
% -----------------------------------
for t = 1:T 
    x(t,:) = [1 k(t,:) a(t,:)];   % The basis functions of the first-degree 
                                  % polynomial of at time t
    k(t+1,:) = x(t,:)*bk_1d;      % Compute next-period capital using bk_1d
end

% 9.2 Compute consumption series 
%-------------------------------
C = (A*k(1:T,:).^alpha.*a(1:T,:) - k(2:T+1,:)+k(1:T,:)*(1-delta))*ones(N,1);
% Aggregate consumption is computed by summing up individual consumption, 
% which in turn, is found from the individual budget constraints; T-by-1

c = C*ones(1,N)/N;                % Individual consumption is the same for
                                  % all countries; T-by-N 

% 9.3 Evaluate the percentage (unit-free) difference between the series  
% from the previous and current iterations
% ---------------------------------------------------------------------
dif_GSSA_1d = mean(mean(abs(1-k./k_old)))
                % Compute a unit-free difference between the capital series 
                % from two iterations; see condition (10) in JMM (2011)
                   
% 9.4 Monte Carlo realizations of the right side of the Euler equation, Y, 
% in condition (C4) in the online Appendix C
%-------------------------------------------------------------------------
for j = 1:N
   Y(1:T-1,j) = beta*c(2:T,j).^(-gam)./c(1:T-1,j).^(-gam).*(1-delta+alpha*A*k(2:T,j).^(alpha-1).*a(2:T,j)).*k(2:T,j);
end  % (T-1)-by-N

% 9.5 Compute and update the coefficients of the capital policy functions 
% -----------------------------------------------------------------------
bk_hat_1d  = (x(1:T-1,:)'*x(1:T-1,:))\x(1:T-1,:)'*Y(1:T-1,:);
                              % Compute new coefficients of the capital 
                              % policy functions using the OLS
bk_1d = kdamp*bk_hat_1d + (1-kdamp)*bk_1d; 
                              % Update the coefficients of the capital  
                              % policy functions using damping 
                                     
% 9.6 Store the capital series 
%-----------------------------
k_old = k;         % The stored capital series will be used for checking 
                   % the convergence on the subsequent iteration
end;

% 10. Time needed to compute the initial guess 
% --------------------------------------------
time_GSSA_1d  = toc; 

% 11.  Obtain the simulated data from the GSSA solution
% -----------------------------------------------------
k_prime_GSSA = x*bk_1d;     % Simulate the GSSA solution for next-period  
                            % capital; this will be needed in Section 16.1  
                            % (in order to find an initial vector of polynomial 
                            % coefficients for different levels of Smolyak 
                            % approximation) 

simul_GSSA = [k(2:end,:) a];% The simulated series for the state variables 
                            % (capital is given by a GSSA solution)
min_simul = min(simul_GSSA);% Compute the minimum values of simulated series;
                            % 1-by-d
max_simul = max(simul_GSSA);% Compute the maximum values of the simulated series
                            % 1-by-d

min_Cheb = ones(1,d).*-1;   % Lower bound of unidimensional Chebyshev poly-
                            % nomial is -1; create a vector of lower bounds 
                            % on Chebyshev polynomials for all dimensions;
                            % 1-by-d
max_Cheb = ones(1,d).*1;    % Upper bound of unidimensional Chebyshev poly-
                            % nomial is 1; create a vector of upper bounds 
                            % on Chebyshev polynomials for all dimensions;
                            % 1-by-d
simul_norm = Rescale(simul_GSSA,min_Cheb,max_Cheb,min_simul,max_simul);
                            % Rescale the simulated series to be in the 
                            % interval [-1, 1]; this will be used in Section 
                            % 14.1 and is needed because Chebyshev polynomials 
                            % are defined in this interval                           

% _________________________________________________________________________                              
%
% Compute a solution to the multi-country model using the anisotropic version 
% of the Smolyak method or the isotropic version as a special case
% 
% The code computes solutions up to approximation level mu_max = 3. To 
% compute solutions of approximation level mu_max = 4, it is necessary
% to use a more accurate initial guess such as a solution under mu_max = 3. 
% Finally, to compute solutions of approximation level mu_max = 5, it
% is necessary to use an accurate initial guess such as a solution under
% mu_max = 4 and to reduce the damping parameter "kdamp" as needed. 
% _________________________________________________________________________


tic;                        % Start counting the time needed to compute the 
                            % Smolyak solution
                            
% 12. The Smolyak algorithm's parameter and grids 
% -----------------------------------------------

  % 12.1 Anisotropic structure 
  %---------------------------
  vector_mus_dimensions = [1,2,3,1,2,3]; 
                            % Introduce a level of approximation in every 
                            % dimension; see Section 4 of JMMV (2014);  
                            % d-by-1, where d=2*N; the first N terms refer 
                            % to the approximation levels of capital in  
                            % countries 1,...,N, and the last N terms refer  
                            % to those of the productivity levels; e.g.,for
                            % N=2, we have approximation levels for [capital 
                            % country 1, capital country 2, productivity 
                            % country 1, productivity country 2];

  mu_max  = max(vector_mus_dimensions);
                            % Compute the maximum level of approximation
                            % across all dimensions

  % 12.2 Construct Smolyak grids and compute the number of polynomial
  % coefficients
  %------------------------------------------------------------------
  Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu_max);
  % Construct the matrix of indices of multidimesional Smolyak elements 
  % (grid points and polynomial basis functions) that satisfy the usual 
  % isotropic Smolyak rule for the approximation level equal to "mu_max"

  Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
  % Select from the matrix of isotropic indices "Smol elem_iso" a subset of 
  % indices that correspond to the given anisotropic "vector_mus_dimensions"
 
  Smol_grid_ani = Smolyak_Grid(d,mu_max,Smol_elem_ani); 
  % Construct the Smolyak grid for the given subindices of anisotropic Smolyak 
  % elements "Smol_elem_ani"
  
  Smol_grid_transf_ani = Rescale(Smol_grid_ani,min_simul,max_simul,min_Cheb,max_Cheb); 
  % Using a linear transformation, transform  the anisotropic Smolyak grid  
  % into the measurement units of the model's state variables (we use the 
  % bounds from the previously simulated GSSA solution)
    
  npol = size(Smol_grid_ani,1);  
  % Compute the number of polynomial coefficients; it is equal to the number 
  % of grid points

  % 12.3 Damping
  %--------------
  kdamp     = 0.1;    % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy functions
  % To achieve convergence under some parameterizations, one may need to 
  % modify the values of the damping parameter kdamp or to refine the 
  % initial guess 
  
% 13. Choose an integration method 
% --------------------------------                             
IM    = 11;      % 1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
                 % nodes in each dimension, respectively;
                 % 11=Monomial rule with 2N nodes;
                 % 12=Monomial rule with 2N^2+1 nodes
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
% Note that under the one-node Gauss-Hermite quadrature rule, the conditional 
% expectation (integral) is approximated by the value of the integrand 
% evaluated in one integration node in which the next-period productivity 
% shock is zero, i.e., the next-period productivity level is 
% a(t+1,:)=a(t,:).^rho*exp(0)=a(t,:).^rho
  
% 14. The remaining choices for the Smolyak algorithm (such as an initial guss 
% on coefficients, data matrix X, grid points for capital and productivity,
% an initial difference in the solutions on two iterations)
% ------------------------------------------------------------------------- 
    % 14.1 Compute the initial guess on the Smolyak polynomial coefficients
    % by using the previously obtained GSSA solution
    % ---------------------------------------------------------------------
    bk_mu = Smolyak_Polynomial(simul_norm,d,mu_max,Smol_elem_ani)\k_prime_GSSA;
    % Recall that k_prime_GSSA=x*bk_1d is next-period capital computed by 
    % GSSA; see Section 11
                   
    % 14.2 Pre-compute the (inverse of) matrix of the polynomial basis 
    % functions evaluated in the grid points 
    %-----------------------------------------------------------------
    X = Smolyak_Polynomial(Smol_grid_ani,d,mu_max,Smol_elem_ani); 
    % Matrix of the polynomial basis functions evaluated in the grid
    % points; npol-by-npol
    
    X_inv = X\eye(npol);
    % Pre-compute the inverse of matrix, where "npol" is the number of the 
    % grid points (polynomial terms) for the given anisotropic approximation
    
    
    % 14.3 Grid points for capital and for productivity 
    %--------------------------------------------------
    k0  =  Smol_grid_transf_ani(:,1:N);   
                        % Grid points for current capital stocks of N 
                        % countries; npol-by-N  
    a0  =  Smol_grid_transf_ani(:,N+1:end);   
                        % Grid points for current productivity levels of N 
                        % countries; npol-by-N  
                        
    % 14.4 Initialize the values of next-period capital of N countries                
    %------------------------------------------------------------------                    
    k_old = ones(npol,N);   
                        % This variable will be used to check the convergence  
                        % on the subsequent iteration (initially, capital  
                        % can take any value); npol-by-N
                        
    % 14.5 Choose an initial difference between the solutions for capital  
    % from two iterations in the convergence criterion                    
    %----------------------------------------------------------------------                   
    dif_Smol_mu  = 1e+10;
                        % Set it to a very large number
                        
     
  % 15 The main iterative cycle of the Smolyak algorithm
  % ------------------------------------------------------  
    while dif_Smol_mu > 1e-10; % 10^(-10) is a convergence criterion                     
    
    % 15.1 Generate time series of capital
    % ---------------------------------------
    k1 = X*bk_mu;
    % Compute next-period capital using polynomial coefficients bk_mu

    % 15.2 Compute consumption series of all countries 
    %--------------------------------------------------- 
    C = (A*k0.^alpha.*a0 - k1+k0*(1-delta))*ones(N,1);
    % Aggregate consumption is computed by summing up individual consumption, 
    % which in turn, is found from the individual budget constraints; npol-by-1
    
    c = C*ones(1,N)/N;            
    % Individual consumption is the same for all countries; npol-by-N 
 
    % 15.3 Approximate the conditional expectations using the integration 
    % method chosen
    %--------------------------------------------------------------------
    % Deterministic integration methods approximate the values of 
    % conditional expectations, Y, in the Euler equation as a weighted average 
    % of the values of the integrand in the given nodes with the given weights 
        Y = zeros(npol,N); % Allocate memory for the variable Y
        for i = 1:n_nodes         
            a1  =  a0.^rho.*exp(ones(npol,1)*epsi_nodes(i,:));   
            % Compute the next-period productivity levels for each integration       
            % node using the process for productivity; n_nodes-by-N
            X2 = Smolyak_Polynomial( Rescale([k1 a1],min_Cheb,max_Cheb,min_simul,max_simul),d,mu_max,Smol_elem_ani); 
            % Matrix of the polynomial basis functions evaluated in the grid 
            % points; npol-by-npol
            k2  = X2*bk_mu;  
            % Compute capital of period t+2 (chosen at t+1) using the
            % capital policy functions; n_nodes-by-N 

            C1 = (A*k1.^alpha.*a1 - k2+k1*(1-delta))*ones(N,1);
            % C is computed by summing up individual consumption, which in
            % turn, is found from the individual budget constraints; npol-by-1

            c1 = C1*ones(1,N)/N;                 
            % Compute next-period consumption for N countries; npol-by-N
      
            for j = 1:N
                Y(:,j) = Y(:,j)+weight_nodes(i,1)*beta*c1(:,j).^(-gam)./c(:,j).^(-gam).*(1-delta+alpha*A*k1(:,j).^(alpha-1).*a1(:,j)).*k1(:,j);
            end % npol-by-N

         end
 
    % 15.4 Evaluate the percentage (unit-free) difference between the 
    % capital series from the previous and current iterations
    % -----------------------------------------------------------------

    dif_Smol_mu = mean(mean(abs(1-k1./k_old)))
                % Compute a unit-free difference between the solutions for
                % capital from two iterations; see condition (10) in JMM (2011)

    % 15.5 Compute and update the coefficients of the capital policy 
    % functions 
    % ----------------------------------------------------------------
    bk_hat_mu = X_inv*Y;
    % Compute new coefficients of the capital policy functions by solving 
    % the inverse problem using the chosen approximation method
    bk_mu = kdamp*bk_hat_mu + (1-kdamp)*bk_mu; 
    % Update the coefficients of the capital policy functions using damping 
                                     
    % 15.6 Store the capital series 
    %--------------------------------
    k_old = k1;      % The stored capital series will be used for checking 
                     % the convergence on the subsequent iteration

    end;

     % 15.7 The Smolyak algorithm's output
     % --------------------------------------
     time_Smol = toc;            % Time needed to compute the solution

     
% 16. Accuracy test of the GSSA solutions: errors on a stochastic simulation 
% --------------------------------------------------------------------------

     % 16.1 Specify a set of points on which the accuracy is evaluated
     %----------------------------------------------------------------
     T_test = 10200;             % Choose the simulation length for the test 
                                 % on a stochastic simulation, T_test<=10,200 

     a_test = a20200(T+1:T+T_test,1:N); 
                                 % Restrict the series of the productivity 
                                 % levels for the test on a stochastic 
                                 % simulation to the given T_test<=10,200  
                                 % and N<=10                             
          
     k_test(1,1:N) = 1;          % Initial condition for capital (equal to 
                                 % steady state)

     % 16.2 Choose an integration method for evaluating accuracy of solutions
     %-----------------------------------------------------------------------
     IM_test = 11;               % See Section 13 for the integration 
                                 % options 

     % To implement the test on a stochastic simulation with T_test>10,200,
     % one needs to simulate new series of the productivity levels with 
     % larger T_test  by enabling the code in paragraph 5.

     % 16.3 Compute errors on a stochastic simulation for the Smolyak 
     % polynomial solution 
     % --------------------------------------------------------------
     % 16.3.1 Simulate the time series solution under the given capital-
     % policy-function coefficients, bk_mu  
     %------------------------------------------------------------------
     bk = bk_mu;                 % The vector of coefficients of the capital
                                 % decision functions of N countries;
                                 % npol-by-N
                          
    for t = 1:T_test-1   
    k_test(t+1,:)  =  Smolyak_Polynomial( Rescale([k_test(t,:) a_test(t,:)],min_Cheb,max_Cheb,min_simul,max_simul),d,mu_max,Smol_elem_ani)*bk; 
    end    
    
    % 16.3.2 Errors across 10,000 points on a stochastic simulation
    % -------------------------------------------------------------
    discard = 200; % Discard the first 200 observations to remove the effect
                   % of the initial conditions 
    [Errors_mean(1),Errors_max(1), time_test(1)] = Accuracy_Test_Smolyak(k_test,a_test,bk,IM_test,alpha,gam,delta,beta,A,tau,rho,vcv,discard,min_Cheb,max_Cheb,min_simul,max_simul,d,mu_max,Smol_elem_ani);
    % Errors_mean    is the unit-free average absolute approximation error  
    %                across 4N+1 optimality conditions (in log10) 
    % Errors_max     is the unit-free maximum absolute approximation error   
    %                across 4N+1 optimality conditions (in log10) 


% 17. Display the results for the polynomial solutions of the degrees from 
% one to mu_max  
% ------------------------------------------------------------------------
disp(' '); disp('           SMOLYAK OUTPUT:'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); 
disp(time_Smol);
disp('b) for implementing the accuracy test'); 
disp(time_test);
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error across 4N+1 optimality conditions'); 
disp(Errors_mean)
disp('b) max error across 4N+1 optimality conditions'); 
disp(Errors_max)

% save Results_N time_Smol time_test Errors_mean Errors_max kdamp RM IM N T bk_mu k_test a_test IM_test alpha gam delta beta A tau rho vcv discard npol mu_max T_test ;
% save Results_Smolyak