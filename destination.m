function [ID] = destination(Ai, INC)
%Create 'destination' array
ID = zeros(2,length(INC(1,:))); %Initialize ID array
P = 0; %Initialize equation number
for A = 1:length(INC(1,:)) %Loop through global functions
    P = P + 1; %Increment equation number
    ID(1,A) = P; %Assign equation number to global function
    P = P + 1;
    ID(2,A) = P;
    for e = 1:length(Ai(:,1)) %Loop through boundary conditions
        if (Ai(e,1) == A) %If current global function has bc
            ID(Ai(e,2),A) = 0; %Set equation number for current global function to 0
            P = P - 1; %Reduce equation number to compensate
            if Ai(e,2) == 1 %If in 1-direction, correct equation number for 2-direction
                ID(2,A) = P;
            end
        end
    end
end

