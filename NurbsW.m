function [W] = NurbsW(Nu, Nv, w)
%Calculate NURBS weighting function from B-spline basis functions and weights
steps = length(Nu(1,:));
W = zeros(steps,steps);
for j = 1:steps
for k = 1:steps %Coordinates for B-spline curve
    W(j,k) = Nu(:,j)'*w*Nv(:,k);   
end
end
end

