function f = ODESysFunS(~, X, Mact, Minh, NumOfNodes)
% ODESysFun: Mendoza-type ODEs with internal alpha/beta weights.
% X is the vector of node activities.
%
% Parallel-safe version:
%   All required model data are passed explicitly instead of using evalin.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decay constant (assume 1)
gamma = ones(1, NumOfNodes);

% Steepness of activation
h = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIFORM α AND β VALUES (same for all edges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_val = 1;    % activation weight
beta_val  = 1;    % inhibition weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate
w = zeros(1, NumOfNodes);
f = zeros(1, NumOfNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xt = X(:)';   % ensure row vector for element-wise products

for i = 1:NumOfNodes
    
    % Binary adjacency
    ActMask = Mact(i,:);
    InhMask = Minh(i,:);
    
    hasAct = any(ActMask);
    hasInh = any(InhMask);
    
    % Weighted versions
    Ract = alpha_val * ActMask;
    Rinh = beta_val  * InhMask;
    
    if (~hasInh && hasAct)
        % ACTIVATION ONLY
        sum_alpha   = sum(Ract);
        sum_alpha_X = sum(Ract .* Xt);
        
        w(i) = ((1 + sum_alpha)/(sum_alpha)) * (sum_alpha_X/(1 + sum_alpha_X));
        
    elseif (~hasAct && hasInh)
        % INHIBITION ONLY
        sum_beta   = sum(Rinh);
        sum_beta_X = sum(Rinh .* Xt);
        
        w(i) = 1 - ((1 + sum_beta)/sum_beta) * (sum_beta_X/(1 + sum_beta_X));
    
    elseif (hasAct && hasInh)
        % BOTH ACTIVATORS AND INHIBITORS
        sum_alpha   = sum(Ract);
        sum_beta    = sum(Rinh);
        sum_alpha_X = sum(Ract .* Xt);
        sum_beta_X  = sum(Rinh .* Xt);
        
        w(i) = (((1+sum_alpha)/sum_alpha) * (sum_alpha_X/(1+sum_alpha_X))) * ...
               (1 - ((1+sum_beta)/sum_beta) * (sum_beta_X/(1 + sum_beta_X)));
    else
        % NO REGULATORS
        w(i) = 0;
    end
    
    % FINAL ODE
    f(i) = (-exp(0.5*h) + exp(-h*(w(i)-0.5))) / ...
           ((1 - exp(0.5*h)) * (1 + exp(-h*(w(i)-0.5)))) ...
           - gamma(i)*X(i);
end

f = f(:);  % ensure column vector
end