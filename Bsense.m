function [dg0,dg1,vol] = Bsense(d,U,V,INC,IEN,LM,gp,gw,p,q,Bx,By,w,D,C,mx,my)
%Calculate sensitivities of objective function and volume constraint for
%shape optimization, as well as total volume of the structure

na = length(IEN(1,:));
dg0 = zeros(max(C(:,1)),1);
dg1 = zeros(max(C(:,1)),1);
vol = 0;

for h = 1:length(C(:,1))
dof = C(h,1);
e = C(h,2);
a = C(h,3);
dir = C(h,4);
dx=mx;

if dir == 2
    dx = my;
end

DX = zeros(na,2);

DX(a,dir) = dx;

neq = sum(LM(:,e) ~= 0);
Ke = zeros(neq,neq);
indx = nonzeros(LM(:,e));
A = IEN(e,:);

for i = 1:length(gp)
    for j = 1:length(gp)
        [~,dR_dx,J] = Shape_function(gp(i),gp(j),e,p,q,Bx,By,w,U,V,INC,IEN);
        G=dR_dx';
        dG_da = -G*DX*G;
        Jmod = J*gw(i)*gw(j);
        dJ_da = Jmod*sum(diag(G*DX));
        dg1(dof) = dg1(dof) + dJ_da;
        Ke = Ke + Ksense(IEN,dR_dx,dG_da',Jmod,dJ_da,D,e,LM);     
    end
end

dg0(dof) = dg0(dof) -d(indx)'*Ke*d(indx);



end

E = unique(C(:,2));

for h = 1:length(E)
    e = E(h);
    
    for i = 1:length(gp)
        for j = 1:length(gp)
            [~,~,J] = Shape_function(gp(i),gp(j),e,p,q,Bx,By,w,U,V,INC,IEN);
            Jmod = J*gw(i)*gw(j);
            vol = vol + Jmod;
        end
    end
end
end

