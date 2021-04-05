function [U, Bxn, Byn, nu, wn] = knotIns(U, Un, p, n, Bx, By, dir, w)
%Inserts knot into knot vector and calculates new control points
Bx = Bx.*w;
By = By.*w;
if dir == 2
   Bx = Bx';
   By = By';
   w = w';
end
dim = size(Bx);

Bxn = zeros(dim(1)+1,dim(2));
Byn = zeros(dim(1)+1,dim(2));
wn = zeros(dim(1)+1,dim(2));



for i = 1:dim(2)
s = FindSpan(n,p,Un,U);
Bxn(1:s-p+1,i) = Bx(1:s-p+1,i);
Byn(1:s-p+1,i) = By(1:s-p+1,i);
wn(1:s-p+1,i) = w(1:s-p+1,i);
    for j = s-p+2:s+1
    a = (Un-U(j))/(U(j+p)-U(j));
    Bxn(j,i) = a*Bx(j,i) + (1-a)*Bx(j-1,i);
    Byn(j,i) = a*By(j,i) + (1-a)*By(j-1,i);
    wn(j,i) = a*w(j,i) + (1-a)*w(j-1,i);
    end
    Bxn(s+2:end,i) = Bx(s+1:end,i);
    Byn(s+2:end,i) = By(s+1:end,i);
    wn(s+2:end,i) = w(s+1:end,i);
end

if(dir == 2)
    Bxn = Bxn';
    Byn = Byn';
    wn = wn';
end

run = true;
if (Un==U(1))
    U = [Un, U];
    run = false;
end
j = 0;
while run
    j = j + 1;
    if (Un<=U(j))
        U = [U(1:j-1), Un, U(j:end)];
        run = false;
    end
end

Bxn = Bxn./wn;
Byn = Byn./wn;
nu = length(U) - p - 1;



