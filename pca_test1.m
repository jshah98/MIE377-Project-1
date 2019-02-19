clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices
adjClose = readtable('Project1_Data_adjClose.csv');
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Date = [];

% Load the factors weekly returns
factRet = readtable('Project1_Data_FF_factors.csv');
factRet.Properties.RowNames = cellstr(datetime(factRet.Date));
factRet.Date = [];

riskFree = factRet(:,4);
factRet = factRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factRet.Properties.RowNames));
% Initiate counter for the number of observations per investment period
toDay = 0;


    periodReturns = table2array( returns);
    periodFactRet = table2array( factRet);

  returns = periodReturns;
  factRet = periodFactRet;
%%%End of Comment out section
X =pca(returns);