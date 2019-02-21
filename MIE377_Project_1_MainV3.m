%% MIE377 (Winter 2019) - Project 1
% The purpose of this program is to implement the following factor models
% a) CAPM
% b) Fama-French 3-factor model
% c) Principal component analysis model
% 
% and to use these models to estimate the asset expected returns and 
% covariance matrix. These parameters will then be used to test the 
% out-of-sample performance of two portfolio optimization models: 
% 1. MVO
% 2. MVO with cardinality constraint
%
% Use can use this template to write your program.
%
% Student Name:
% Student ID:

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
factorRet = readtable('Project1_Data_FF_factors.csv');
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
factorRet.Date = [];

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Factor models
% Note: You must populate the functios CAPM.m, FF.m and PCA.m with your
% own code.
FMList = {'CAPM' 'FF' 'PCA'};
FMList = cellfun(@str2func, FMList, 'UniformOutput', false);
NoModels = length(FMList);

% Investment strategies
% Note: You must populate the functios MVO.m and MVO_card.m with your own
% code to construct the optimal portfolios. 

invList = {'MVO' 'MVO_card'};
invList = cellfun(@str2func, invList, 'UniformOutput', false);
NoStrats = length(invList);

% Cardinality constraint: state the maximum number of assets
card = 12;

% Tags for factor models under different investment strategies
tags = {'MVO (CAPM)' 'Card MVO (CAPM)' 'MVO (FF)' 'Card MVO (FF)' ...
        'MVO (PCA)' 'Card MVO (PCA)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns, cov. matrix,
% etc) from the Fama-French factor models. You will have to re-estimate 
% your parameters at the start of each rebalance period, and then 
% re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        currentVal= ones(1, 6);
        currentVal(t,:) = initialVal;
        
        
    else
        for k = 1 : (NoStrats * NoModels)
            
            currentVal(t,k) = currentPrices .* NoShares{k};
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % Calculate 'mu' and 'Q' using the 3 factor models.
    % Note: You need to write the code for the 3 factor model functions. 
    for i = 1 : NoModels
        
        [mu{i}, Q{i}] = FMList{i}(periodReturns, periodFactRet);
        
    end
            
    % Optimize your portfolios to get the weights 'x'
    % Note: You need to write the code for the 2 portfolio optimization 
    % functions.
    for i = 1 : NoModels
        
        % Define the target return
        targetRet = mean(mu{i});
        
        for j = 1: NoStrats
            
            k = j + (i - 1) * 2;
            
            x{k}(:,t) = invList{j}(mu{i}, Q{i}, targetRet, card); 
            
        end
    end
    
    % Calculate the optimal number of shares of each stock you should hold
    for k = 1 : NoStrats * NoModels
        % Number of shares your portfolio holds per stock
        NoShares{k} = x{k}(:,t) .* currentVal(t,k) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,k) = periodPrices * NoShares{k};
        
        %------------------------------------------------------------------
        % Calculate your transaction costs for the current rebalance
        % period. The first period does not have any cost since you are
        % constructing the portfolios for the 1st time. 
        
        if t ~= 1
           
           tCost(t-1, k) =  (x{k}(:,t)- x{k}(:,t-1) )*0.005 * periodPrices(t,:)
            
        end
        
        NoSharesOld{k} = NoShares{k};
        %------------------------------------------------------------------
        
    end

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 4.1 Calculate the portfolio average return, variance (or standard 
% deviation), and any other performance and/or risk metric you wish to 
% include in your report.
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% 4.2 Plot the portfolio values 
% 
% Note: The code below plots all portfolios onto a single plot. However,
% you may want to split this into multiple plots for clarity, or to
% compare a subset of the portfolios. 
%--------------------------------------------------------------------------
plotDates = dates(dates >= testStart);

fig1 = figure(1);

for k = 1 : NoModels * NoStrats
    
    plot( plotDates, portfValue(:,k) )
    hold on
    
end

legend(tags, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'fileName','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.3 Plot the portfolio weights 
%--------------------------------------------------------------------------

% MVO (CAPM) Plot
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('MVO (CAPM) portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'fileName2','-dpng','-r0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End