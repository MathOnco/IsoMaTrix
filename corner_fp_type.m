function cfp_val = corner_fp_type(A,strat_num)
%cfp_val = corner_fp_type(A,strat_num)
%   Takes as input:
%       A - square 3-by-3 game matrix
%       strat_num - strategy to analyze (1, 2, or 3; default: 1)
%   Returns the stability value (cfp_val) of monomorphic population of type
%   strat_num when invaded by any combination of other two strategies. 
%   Possible values for cfp_val:
%       2 - source
%       1 - semi-source
%       0 - saddle
%      -1 - semi-sink
%      -2 - sink
%      NaN - neutral

%check strat_num and number of inputs
if (nargin == 1)
    strat_num = 1;
elseif (nargin == 0)
    fprintf('Corner analysis failed: provide a 3-by-3 game matrix A \n')
    return
else
    if ((strat_num ~= 3) && (strat_num ~= 2) && (strat_num ~= 1))
        fprintf('Corner analysis warning: unclear strat_num, defaulting to 1')
        strat_num = 1; %default to analyzing first strategy
    end
end

%check matrix input
sA = size(A);
if (sA(1,1) ~= 3) || (sA(1,2) ~= 3)
    fprintf('Corner analysis failed: Game matrix A has to be 3-by-3 \n')
    return
end

%check if matrix is degenerage (i.e. equivalent to the neutral zero-game)
if min(A,[],'all') == max(A,[],'all')
    fprintf('Corner analysis warning: A is a neutral game')
    cfp_val = NaN;
    return
end

%transform matrix if analyzing strat 2 or 3
if strat_num == 2
    A = [0, 1, 0; 1, 0, 0; 0, 0, 1]*A*[0, 1, 0; 1, 0, 0; 0, 0, 1]'; %rotate the matrix
elseif strat_num == 3
    A = [0, 0, 1; 0, 1, 0; 1, 0, 0]*A*[0, 0, 1; 0, 1, 0; 1, 0, 0]'; %rotate the matrix
end

%%ANALYSIS STARTS

fp_val_2 = edge_fp_type(A([1,2],[1,2]));
fp_val_3 = edge_fp_type(A([1,3],[1,3]));

if (abs(fp_val_2) == 2) || (abs(fp_val_3) == 2) 
    %if the corner is a strict sink or source along some edge
    cfp_val = round((fp_val_2 + fp_val_3)/2);
%otherwise both edges are not first-order stable or invadable
%this means the first row of A has equal entries
elseif fp_val_2*fp_val_3 < 0 
    %if one edge is second-order sink & other second-order source
    cfp_val = 0; %then corner is a saddle point
else
    %if both edges are second-order sinks or both second-order sources
    %  or one of them is neutral and other second-order
    %  or both are neutral
    inv_fit = ...%fitness fn for invader mix of p of type 2 & (1 - p) of type 3
        [A(2,2) - A(2,3),A(2,3),0] ... %p*f_2
        + ...
        ([0, A(3,2) - A(3,3),A(3,3)] - [A(3,2) - A(3,3),A(3,3),0]); %(1 - p)*f_3
    res_fit = [0,A(1,2) - A(1,3),A(1,3)]; %fitness fn for resident (f_1)
    
    %we need to check if inv_fit - res_fit changes signs
    gain_roots = roots(inv_fit - res_fit);
    
    if length(gain_roots) <= 1 %if no quadratic term in inv_fit - res_fit
        cfp_val = fp_val_2 + fp_val_3; %then same stability as edge
        if (fp_val_2 == 0) && (fp_val_3 == 0) %if both neutral
            cfp_val = NaN;
        end
    elseif isreal(gain_roots) %if roots are real
        roots_in_range = (0 < gain_roots) & (gain_roots < 1);
        if max(roots_in_range) %is at least one root for achieavable p in (0,1)
            if gain_roots(1) == gain_roots(2) %if double root
                %then semi-stable if at least one stable, or semi-sink if at least is sink
                cfp_val = sign(fp_val_2 + fp_val_3) ;
            else %distinct roots
                cfp_val = 0; %saddle point
            end
        elseif ((min(gain_roots) == 0) || (max(gain_roots) == 1))
            %if degenerate case touching at 0 or 1
            %then semi-source or semi-sink based on sign of invader's quadratic term
            cfp_val = sign(inv_fit(1))
        else %no internal root is achievable
            cfp_val = fp_val_2 + fp_val_3;
            if (fp_val_2 == 0) && (fp_val_3 == 0) %if both neutral
                cfp_val = NaN;
            end
        end
    else %roots are imaginary, no real crossing of 0
        cfp_val = fp_val_2 + fp_val_3;
    end
end

end
