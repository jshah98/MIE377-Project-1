function  [mu, Q] = CAPM(returns, factRet)
    
    % Use this function to perform a linear regression model for CAPM.

% Number of observations;
N = size(returns, 1);
%Number of columns in returns
n = size(returns,2);
%Number of columns in factor table
nFactCols = size(factRet,2);
% Calculate the factor expected excess return from historical data using
% the geometric mean
expExFactRet = geomean(factRet + 1) - 1;
% Calculate the factor variance
sigmaF = var(factRet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Beta estimate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function regress(y,X) requires a response vector y and a predictor
% matrix X. The matrix X must include a column of ones to account for the
% intercept alpha.
 X = [ones(N,1) factRet(:,1)];
Use the closed-form (CF) solution to find the collection of alphas 
% and betas for all assets
temp = inv(transpose(X)* X)*transpose(X)*returns;
 alpha  = temp(1,:);
 temp_betaCF = temp(2,:);
 betaCF = transpose(temp_betaCF);
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
epsilon = returns - X*temp;
for i = 1:n
    for t = 1:N
    epsilon(t,i) = (returns(t,i)) - (alpha(i) + betaCF(i)*(expExFactRet(1))); 
    end
   end
% Calculate the residual variance with "N - p - 1" degrees of freedom
denom = N -2;
sigmaEp = zeros(n,1);
for i = 1:n
    sum = 0;
    for t = 1:N
       sum = sum + (epsilon(t,i))^2; 
    end
    variance = sum/denom;
    sigmaEp(i,1) = variance;
end


% Calculate the asset expected excess returns, mu
 mu = transpose(alpha) + betaCF*(expExFactRet(1));

% Calculate the diagonal matrix of residuals and the asset covariance 
% matrix
D = diag(sigmaEp);
temp =factRet(:,1);
 Q = var(temp)*betaCF *transpose(betaCF) + D;

   
end
