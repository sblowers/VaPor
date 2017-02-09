function New = findPoints(Original,Vessel)

New = zeros(1,size(Original,1));
for M = 1:numel(New)
    InletDist = zeros(size(Vessel,1),1);
    for N = 1:size(Vessel,1)
        InletDist(N) = norm(Vessel(N,3:5)-Original(M,3:5));
    end
    index = find(InletDist == min(InletDist));
    New(M) = Vessel(index(1),1);
end

end