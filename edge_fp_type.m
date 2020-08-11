function fp_val = edge_fp_type(A,strat_num)
%fp_val = edge_fp_type(A,strat_num)
%   Takes as input:
%       A - square 2-by-2 game matrix
%       strat_num - strategy to analyze (1 or 2; default: 1)
%   Returns the stability value (fp_val) of monomorphic population of type
%   strat_num when invaded by the other strategy. Possible values for
%   fp_val:
%       2 - strict source
%       1 - source
%       0 - neutral
%      -1 - sink
%      -2 - strict sink

%check strat_num and number of inputs
if (nargin == 1)
    strat_num = 1;
elseif (nargin == 0)
    fprintf('Edge analysis failed: provide a 2-by-2 game matrix A \n')
    return
else
    if (strat_num ~= 2) && (strat_num ~= 1)
        fprintf('Edge analysis warning: unclear strat_num, defaulting to 1')
        strat_num = 1; %default to analyzing first strategy
    end
end


%check matrix input
sA = size(A);
if (sA(1,1) ~= 2) || (sA(1,2) ~= 2)
    fprintf('Edge analysis failed: Game matrix A has to be 2-by-2 \n')
    return
end

%transform matrix if analyzing strat 2
if strat_num == 2
    A = [0, 1; 1, 0]*A*[0, 1; 1, 0]; %rotate the matrix
end

%%ANALYSIS STARTS

inv_fit = A(2,1) - A(1,1); %invasion fitness of type 2 against 1

if inv_fit > 0
    fp_val = 2;
elseif inv_fit < 0
    fp_val = -2;
else %neutral first-order, so source/sink not strict
    %check if type 2 can invade type 1 by drift (i.e. not strict)
    drift_fit = A(2,2) - A(1,2);
    
    if drift_fit > 0
        fp_val = 1;
    elseif drift_fit < 0
        fp_val = -1;
    else
        fp_val = 0;
    end 
end

end

