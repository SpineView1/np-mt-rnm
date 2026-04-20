function f = ODESysFun(~, X)
% ODESysFun: Mendoza-type ODEs with internal alpha/beta weights.
% X is the vector of node activities.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of Nodes
NumOfNodes = evalin('base','NumOfNodes');

% Decay constant (assume 1)
gamma = ones(1, NumOfNodes);

% Steepness of activation
h = 10;

% Activation / inhibition adjacency matrices
Mact = evalin('base','Mact');  
Minh = evalin('base','Minh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: UNIFORM α AND β VALUES (same for all edges)
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

for i = 1:NumOfNodes
    
    % Binary adjacency
    ActMask = Mact(i,:);
    InhMask = Minh(i,:);
    
    hasAct = any(ActMask);
    hasInh = any(InhMask);
    
    % Weighted versions
    Ract = alpha_val * ActMask;   % ← uniform α
    Rinh = beta_val  * InhMask;   % ← uniform β
    
    if (~hasInh && hasAct)
        % ACTIVATION ONLY
        sum_alpha   = sum(Ract);
        sum_alpha_X = sum(Ract .* X');
        
        % --- ORIGINAL FORMULA ---
        w(i) = ((1 + sum_alpha)/(sum_alpha)) * (sum_alpha_X/(1 + sum_alpha_X));
        
        % --- KEEPING YOUR ORIGINAL COMMENTED OPTIONS ---
        % w(i) = some alternate activation-only formula...
    
    elseif (~hasAct && hasInh)
        % INHIBITION ONLY
        sum_beta   = sum(Rinh);
        sum_beta_X = sum(Rinh .* X');
        
        % --- ORIGINAL FORMULA ---
        w(i) = 1 - ((1 + sum_beta)/sum_beta) * (sum_beta_X/(1 + sum_beta_X));

        % --- KEEPING YOUR ORIGINAL COMMENTED OPTIONS ---
        % w(i) = 1 - other inhibition-only variant...
    
    elseif (hasAct && hasInh)
        % BOTH ACTIVATORS AND INHIBITORS
        sum_alpha   = sum(Ract);
        sum_beta    = sum(Rinh);
        sum_alpha_X = sum(Ract .* X');
        sum_beta_X  = sum(Rinh .* X');
        
        % --- YOUR SELECTED WORKING FORMULA ---
        w(i) = (((1+sum_alpha)/sum_alpha) * (sum_alpha_X/(1+sum_alpha_X))) * ...
               (1 - ((1+sum_beta)/sum_beta) * (sum_beta_X/(1 + sum_beta_X)));
        %w(i) = (((1+sum_alpha)/sum_alpha) * (sum_alpha_X/(1+sum_alpha_X))) * ...
        %       (1 - ((1+sum_beta)/sum_beta) * (sum_beta_X/(1 + sum_beta_X + sum_alpha_X)));
        
        % --- YOUR ORIGINAL COMMENTED-OUT OPTIONS : KEPT EXACTLY ---
        % w(i)=(((1+sum_alpha)/(sum_alpha))*((sum_alpha_X)/(1+sum_alpha_X)))*(1-((1+sum_beta)/(sum_beta))*((sum_beta_X)/(1+sum_beta_X)));
         %w(i)=(((1+sum_alpha)/(sum_alpha))*((sum_alpha_X)/(1+sum_alpha_X)))*(1-((sum_beta_X)/(sum_beta_X+sum_alpha_X))*((1+sum_beta)/(sum_beta))*((sum_beta_X)/(1+sum_beta_X)));
        % w(i)=((sum_alpha_X)/(sum_beta+sum_alpha))*(((1+sum_alpha)/(sum_alpha))*((sum_alpha_X)/(1+sum_alpha_X)))*(1-((sum_beta_X)/(sum_beta+sum_alpha))*((1+sum_beta)/(sum_beta))*((sum_beta_X)/(1+sum_beta_X)));

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
