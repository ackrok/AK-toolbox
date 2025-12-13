function [combos] = getUnitCombos(units_m, units_n)
%getUnitCombos: Obtain all possible combinations of unit numbers between
%two sets of units. Default is set to avoid repeat combinations and to
%remove self-combinations. 

[a,b] = meshgrid(units_m, units_n);
c = cat(2,a',b'); 
d = reshape(c,[],2);
d = unique(sort(d,2),'rows'); %Remove repeat combinations 
d(d(:,1) == d(:,2),:) = []; %Remove self-combinations
combos = d; 

end

