function [ders] = BasisFuns2(u,p,n,k,U) % From page 70 in NURBs book
% Finds value of basis functions at u and values of derivatives
i = FindSpan(n,p,u,U);
N = zeros(p+1,p+1);
ders = zeros(p+1,k+1);
a = zeros(2,p+1);
N(1,1) = 1;
left = zeros(1,p+1);
right = zeros(1,p+1);
for j = 2:p+1
   left(j) = u - U(i+3-j);
   right(j) = U(i+j) - u;
   saved = 0;
   
   for r = 1:j-1
      N(j,r) = (right(r+1)+left(j+1-r));
      temp = N(r,j-1)/N(j,r);
      N(r,j) = saved + right(r+1)*temp;
      saved = left(j+1-r)*temp;
   end
   N(j,j) = saved;
end

ders(:,1) = N(:,p+1);

for r = 1:p+1
    s1 = 1;
    s2 = 2;
    a(1,1) = 1;
    for x = 1:k
        d = 0;
        rk = r-x;
        pk = p-x;
        if (r >= x+1)
            a(s2,1) = a(s1,1)/(N(pk+2,rk));
            d = a(s2,1)*N(rk,pk+1);
        end
        if (rk >= 0)
            j1 = 2;
        else
            j1 = 2 - rk;
        end
        if (r-2 <= pk)
            j2 = x;
        else
            j2 = p-r+2;
        end
        
        for j = j1:j2
            a(s2,j) = (a(s1,j)-a(s1,j-1))/(N(pk+2, rk+j-1));
            d = d + a(s2,j)*N(rk+j,pk+1);
        end
        if (r-1 <= pk)
            a(s2,x+1) = -a(s1,x)/N(pk+2,r);
            d = d + a(s2,x+1)*N(r,pk+1);
        end
        ders(r,x+1) = d; %Error?
        j = s1;
        s1 = s2;
        s2 = j;
    end
end

r = p;
for x = 1:k
    for j = 1:p+1
        ders(j,x+1) = ders(j,x+1)*r;
    end
    r = r*(p-x);
end

end
