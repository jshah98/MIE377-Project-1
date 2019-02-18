%% MIE377 (Winter 2019) - Tutorial 2
% The purpose of this program is to implement a single-factor model in
% Matlab. We will be implementing the CAPM and using the estimated
% parameters for portfolio optimization. 
%
% TA: Ricardo Pillaca
% Instructor: Giorgio Costa

clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1: Data pre-processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the sample historical data
load('lab2data.mat')

% Calculate the asset and factor returns (factor models use returns, not
% prices)
rets    = prices( 2:end, : ) ./ prices( 1:end - 1, : ) - 1;
facRets = sp500price( 2:end , 1 ) ./ sp500price( 1:end - 1, 1 ) - 1;

% Number of assets
n = size(rets,2);

% Number of observations;
N = size(rets, 1);

% Define the risk-free rate
rf= 0.001;

% CAPM requires the use of excess returns
exRets    = rets - rf;
exFacRets = facRets - rf;

% Calculate the factor expected excess return from historical data using
% the geometric mean
expExFacRets = geomean(exFacRets + 1) - 1;

% Calculate the factor variance
sigmaF = var(exFacRets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Beta estimate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate beta using the three alterntive methods (use the excess returns)

% 2.1 beta = sigma_im / sigma_m^2
betaCAPM = zeros(n,1);

for i = 1 : n
     
     covMat = cov(exRets(:,i),exFacRets);
     %tempTranspose = transpose(temp);
     betaCAPM(i) = covMat(2,1)/sigmaF;
    
end
%display(betaCAPM);
% 2.2 Use the regress(y,X) function from matlab
betaReg = zeros(n,1);

% The function regress(y,X) requires a response vector y and a predictor
% matrix X. The matrix X must include a column of ones to account for the
% intercept alpha.

 X = [ones(210,1) exFacRets];
display(X);
for i = 1 : n
      temp = exRets(:,i);
     b = regress(temp,X);
    betaReg(i) = b(2);
end

% 2.3 Use the closed-form (CF) solution to find the collection of alphas 
% and betas for all assets
temp = inv(transpose(X)* X)*transpose(X)*exRets;
 alpha  = temp(1,:);
 temp_betaCF = temp(2,:);
 betaCF = transpose(temp_betaCF);
% Display your three estimates of beta for all assets
display([betaCAPM betaReg betaCF]);

%Conclude that using the three methods yields the same Betas for all assets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 3: Portfolio optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For portfolio optimization we need the asset expected returns, mu, and
% the covariance matrix, Q. Therefore, we need the following parameters
% from our regression model
% - alpha (intercept)
% - beta (factor loading)
% - factor expected return
% - factor variance
% - residual variance (i.e., idiosyncratic risk)

% Calculate the residuals
epsilon = exRets - X*temp;
for i = 1:n
    for t = 1:210
 epsilon(t,i) = (exRets(t,i) - rf) - (alpha(i) + betaCF(i)*(expExFacRets)); 
    end
   end
% Calculate the residual variance with "N - p - 1" degrees of freedom
% sigmaEp = 

% Calculate the asset expected excess returns, mu
% mu = 

% Calculate the diagonal matrix of residuals and the asset covariance 
% matrix
% D = 
% Q = 

% Setup the MVO problem as we did in Lab 1. Matlab requires the input in 
% the following format
% 
%   min     (1/2) x' H x + c' x
%   s.t.    A x <= b
%           Aeq x = beq
%           lb <= x <= ub
% 
% And, based on our problem, we have that
% 
%   min     x' Q x
%   s.t.    mu' x >= 0.003
%           sum( x ) = 1

% Add the expected return constraint
A = -mu';
b = -0.0025;

%constrain weights to sum to 1
Aeq = ones(1,50);
beq = 1;

% It might be useful to increase the tolerance of 'quadprog'
options = optimoptions('quadprog','TolFun',1e-9);

% Solve this problem using 'quadprog'
x = quadprog( 2 * Q, [], A, b, Aeq, beq, [], [], [], options );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End














