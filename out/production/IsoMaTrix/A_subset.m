function Asub = A_subset(A,types)
    Asub = zeros(3,3);
    for i = 1:3
        for j = 1:3
            Asub(i,j) = A(types(i),types(j));
        end
    end
end



