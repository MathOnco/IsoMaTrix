function outA = Ohtsuki_Nowak_transform(A,k,nuType)
%Ohtsuki_Nowak_transform(A,k,r)
%   Takes as input:
%       A - square game matrix
%       k - degree of graph (> 3)
%       r - choice of dynamics:
%           "DB": death-birth
%           "IM": immitation
%           "BD": birth-death
%   Returns the Ohtsuki-Nowak transform of A for degree k with update rule r.
%
%   For overview of transform see: https://egtheory.wordpress.com/2012/10/25/ohtsuki-nowak-transform/

outA = A; %failure mode, return matrix unchanged

%check that the matrix has the right dimensions
sA = size(A);
if sA(1,1) ~= sA(1,2)
    fprintf('Transform Failed: Game matrix A has to be square \n')
    return
end

%check that k is large enough
if k < 3
    fprintf('Transform Failed: k < 3 is too small \n')
    return
end

%set nu
if nuType == "DB" %death-birth dynamics
    nu = 1/(k + 1);
elseif nuType == "IM" %immitation dynamics
    nu = 3/(k + 3);
elseif nuType == "BD" %birth-death dynamics
    nu = 1;
elseif isnumeric(nuType)
    nu = nuType;
else
    fprintf('Transform Failed: update rule misspecified \n')
    return
end

dA = diag(A); %get the diagonal vector
dO = ones(size(dA)); %make a ones vector of right length

outA = A + (1/(k - 2))*(dA*dO' - dO*dA') + (nu/(k - 2))*(A - A');

end

