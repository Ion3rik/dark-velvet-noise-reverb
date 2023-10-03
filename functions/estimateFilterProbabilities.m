% Function for mapping the analysed magnitude responses to DVN Filter probabilities
% Jon Fagerström
% 13.5.2022

% Arguments:
%   <targetResponses>:  Matrix of target magnitude responses with shape (taps, numResponses)
%   <filterResponses>:  Matrix of filter magnitude responses with shape (taps, numFilters)

% Returns:
%   <p>:                Filter probabilities with shape (numFilters, numResponses)
%   <pSparse>:          Sparse version of p
function [p,pSparse, fittedResponses, gNorm] = estimateFilterProbabilities(targetResponses, filterResponses)

%     % LASSO
%     lambda = 0.3;
%     yMean = mean(targetResponse);
%     xMean = mean(filterResponses);
%     y0 = targetResponse - yMean;
%     X0 = filterResponses - xMean;
%     
%     X0 = X0 ./ std(X0); 
%     
%     pSparse = lasso(X0,y0,'Lambda',lambda,'Standardize',false); 
    
    % Least squares with lower bound 0
    %lb = zeros(size(filterResponses,2),1);
    %p = lsqlin(filterResponses,targetResponse,[],[],[],[],lb,[]);
    numResponses = size(targetResponses,2);
    numFilters = size(filterResponses,2);
    fittedResponses = zeros(size(targetResponses));
    
    p = zeros(numFilters, numResponses);
    for n = 1:numResponses
        p(:,n) = lsqnonneg(filterResponses, targetResponses(:,n));  % solve the optimization problem
        fittedResponses(:,n) = filterResponses*p(:,n);               % save the fitted responses
        gNorm(n) = sum(p(:,n));
       
    end
    pSparse = p;
    pSparse(p<0.0001) = 0;
end

