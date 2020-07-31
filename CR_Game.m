function A = CR_Game(c)
%A = CR_Game(c)
%   Takes as input:
%       c - cost of resistance in [0,1]
%   Returns a 3-by-3 game matrix A of the normal-sensitive-resistant game
%   from West, Ma & Netwon (2018) with a cost of resistance of c as used in
%   the example in the IsoMaTrix manual.  

    %N           %S   %R
A = [1.2,        1,   1;           %N
    1.4 - 0.4*c, 1,   1.1 - 0.1*c; %S
    1.4,         0.7, 1.1];          %Ris
end

