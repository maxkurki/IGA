function [s] = FindSpan(n,p,u,U) %From page 68 in NURBS book
%Finds the knot span relevant for the point u
%n - knot span highest index
%p - base function order
%u - evaluation point
%U - knot vector

if (u == U(n+2))
    s = n-1;
else
    low = p;
    high = n;
    mid = floor((low+high)/2);
    while ((u < U(mid) || u >= U(mid+1)) && low ~= mid)
        if (u < U(mid))
            high = mid;
        else
            low = mid;
        end
        mid = floor((low + high)/2);
    end
    if U(mid) == U(mid+1)
        mid = mid + 1;
    end
    if u >= U(mid+1)
       mid = mid + 1; 
    end
    s = mid-1;
end

