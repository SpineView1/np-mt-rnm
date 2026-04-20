function [Mact, Minh, NodeNames, NumOfNodes, stimuli_names] = CreateMatrices_new(NW2)
% CreateMatrices (fixed)
% ------------------------------------------------------------
% Reads network from Excel table with columns:
%   Nodes, Activators, Inhibitors, Stimuli
% Builds:
%   Mact(i,j)=1 if Stimuli{j} activates Nodes{i}
%   Minh(i,j)=1 if Stimuli{j} inhibits Nodes{i}
%
% Fixes vs old version:
%   - trims whitespace around tokens
%   - case-insensitive matching (strcmpi)
%   - robust handling of empty/missing activator/inhibitor fields
%   - uses cell indexing safely

    T = readtable(NW2);

    % Basic checks (optional but helpful)
    requiredColumns = {'Nodes','Activators','Inhibitors','Stimuli'};
    if ~all(ismember(requiredColumns, T.Properties.VariableNames))
        error('Excel must contain columns: Nodes, Activators, Inhibitors, Stimuli.');
    end

    NodeNames    = T.Nodes;
    stimuli_names = T.Stimuli;

    NumOfNodes   = length(NodeNames);
    NumOfStimuli = length(stimuli_names);

    Mact = zeros(NumOfNodes, NumOfStimuli);
    Minh = zeros(NumOfNodes, NumOfStimuli);

    for i = 1:NumOfNodes

        % ---------------- ACTIVATORS ----------------
        ActInLine = T.Activators(i);
        ActTokens = local_parseTokens(ActInLine);

        PresentInNode = zeros(NumOfStimuli, 1);
        for j = 1:numel(ActTokens)
            PresentInNode = PresentInNode + strcmpi(stimuli_names, ActTokens{j});
        end
        Mact(i,:) = PresentInNode.';

        % ---------------- INHIBITORS ----------------
        InhInLine = T.Inhibitors(i);
        InhTokens = local_parseTokens(InhInLine);

        PresentInNode = zeros(NumOfStimuli, 1);
        for j = 1:numel(InhTokens)
            PresentInNode = PresentInNode + strcmpi(stimuli_names, InhTokens{j});
        end
        Minh(i,:) = PresentInNode.';

    end
end

function tokens = local_parseTokens(cellValue)
% Turns one table cell (string/char/cell/NaN/missing) into a clean cellstr list.
% Splits by comma, trims whitespace, removes empty tokens and placeholder
% entries such as NOTHING / NONE / NA / N/A.

    tokens = {};

    if isempty(cellValue)
        return;
    end

    % Convert to string safely
    try
        s = string(cellValue);
    catch
        s = "";
    end

    % Handle missing/NaN-like entries
    if ismissing(s) || strlength(strtrim(s)) == 0
        return;
    end

    % Split + trim
    parts = strsplit(strtrim(s), ',');
    parts = strtrim(parts);

    % Remove empty tokens
    parts = parts(parts ~= "");

    % Remove placeholder terms meaning "no regulator"
    badTokens = ["NOTHING","NONE","NA","N/A","NULL","-"];
    keepMask = ~ismember(upper(parts), badTokens);
    parts = parts(keepMask);

    tokens = cellstr(parts);
end
