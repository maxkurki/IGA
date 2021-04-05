function [INC, IEN] = connect(n, m, p ,q)
% Create INC and IEN arrays
% From page 317 in Hughes
INC = zeros(2,n*m);
IEN = zeros((n-p)*(m-q), (p+1)*(q+1));
A = 0;
e = 0;
for j = 1:m
    for i = 1:n
        A = A + 1;
        INC(1,A) = i;
        INC(2,A) = j;
        
        if (i > p && j > q)
            pos = 1;
            e = e + 1;
            for jj = 0:q
                for ii = 0:p
                    B = A - jj*n - ii;
                    IEN(e,pos) = B;
                    pos = pos + 1;
                end
            end
        end
    end
end
end
