function [LM] = location(ID, IEN)
%Create 'location' matrix
LM = zeros(2*length(IEN(1,:)),length((IEN(:,1)))); %Initialize location matrix
for a = 1:length(IEN(1,:)) % Loop through local basis functions
    for e = 1:length(IEN(:,1)) % Loop through elements
        p = 2*(a-1) + 1; % Local equation number for x-direction 
        LM(p,e) = ID(1,IEN(e,a)); % Assign global equation number
        p = 2*(a-1) + 2; % --||-- y-direction
        LM(p,e) = ID(2,IEN(e,a));
    end
end

