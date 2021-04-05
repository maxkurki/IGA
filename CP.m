function [C,CP] = CP(INC,IEN,m)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
C = [];
dofnr = 1;
for i = 1:length(INC)
    if INC(2,i) == m 
        for e = 1:length(IEN(:,1))
           for a = 1:length(IEN(1,:)) 
               if(IEN(e,a)==i)
                  C = [C; dofnr,e,a,2];
               end
           end
        end
        dofnr = dofnr + 1;
    end
end

CP = zeros(max(C(:,1)),2);

for i = 1:length(C(:,1))
    CP(C(i,1),:) = [IEN(C(i,2),C(i,3)), C(i,4)]; 
end
end

