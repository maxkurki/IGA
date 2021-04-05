function [N] = BasisFuns(u,p,n,U) % From page 70 in NURBs book
% Finds value of basis functions at u
i = FindSpan(n,p,u,U);
N = zeros(1,p+1);
N(1) = 1;
left = zeros(1,p+1);
right = zeros(1,p+1);
for j = 2:p+1
   left(j) = u - U(i+3-j);
   right(j) = U(i+j) - u;
   saved = 0;
   
   for r = 1:j-1
      temp = N(r)/(right(r+1)+left(j+1-r));
      N(r) = saved + right(r+1)*temp;
      saved = left(j+1-r)*temp;
   end
   N(j) = saved;
end
end

