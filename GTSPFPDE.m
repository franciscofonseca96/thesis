function [Price] = GTSPFPDE(S0,K,r,q,T,Nx,Nt,Cp,Cm,Lp,Lm,Ap,Am)
xmax=log(4*K);
xmin=-xmax;

dx=(xmax-xmin)/(Nx-1);
Xi=(xmin:dx:xmax)';

dt=T/(Nt-1);
tauj=0:dt:T;
nu = Cp*gamma(-Ap)*((Lp-1)^Ap-Lp^Ap+Ap*Lp^(Ap-1)) + Cm*gamma(-Am)*((Lm+1)^Am-Lm^Am-Am*Lm^(Am-1));

a = -dt/(2*dx) * (r - q - nu  + gamma(-Ap)*Cp*Ap*Lp^(Ap-1)-gamma(-Am)*Cm*Am*Lm^(Am-1));
b = dt * (r + gamma(-Ap)*Cp*Lp^(Ap)+gamma(-Am)*Cm*Lm^(Am));
c = dt/(2*dx) * (r - q - nu  + gamma(-Ap)*Cp*Ap*Lp^(Ap-1)-gamma(-Am)*Cm*Am*Lm^(Am-1));

qsi=zeros(Nx,1);
eta=zeros(Nx,1);


for i=1:Nx
    qsi(i,1)=Cm*gamma(-Am)*exp(-Lm*Xi(i))*dt/(dx)^Am;
end

for i=1:Nx
    eta(i,1)=Cp*gamma(-Ap)*exp(Lp*Xi(i))*dt/(dx)^Ap;  
end

gm=ggcoef(Nx,Am);
gp=ggcoef(Nx,Ap);

A1=Tridiag(Nx,b,a,c);

A1(1,1)=1;A1(end,end)=1;A1(1,2)=0;A1(end,end-1)=0;
A2 = zeros(Nx);
for i=2:Nx
    for j=2:Nx
        if j==i
            A2(i,j)=1-(exp(Lm*Xi(j))*qsi(i)*gm(2)+exp(-Lp*Xi(j))*eta(i)*gp(2));
        end
        if j==i-1
            A2(i,j)=-(exp(Lm*Xi(j))*qsi(i)*gm(3)+exp(-Lp*Xi(j))*eta(i)*gp(1));
        end
        if j==i+1
            A2(i,j)=-(exp(Lm*Xi(j))*qsi(i)*gm(1)+exp(-Lp*Xi(j))*eta(i)*gp(3));
        end
        if j<i-1
            A2(i,j)=-exp(Lm*Xi(j))*qsi(i)*gm(i-j+2);
        end
        if j>i+1
            A2(i,j)=-exp(-Lp*Xi(j))*eta(i)*gp(j-i+2);
        end
    end

end

A2(1,:)=0;
A2(end,:)=0;
A2(1,1)=0;
A2(end,end)=0;

A = A1 + A2;

Usol = zeros(Nx,Nt);

U1=max(K-exp(Xi),0);
LBC=(K-exp(xmin))*exp(-r*tauj);
UBC=tauj*0;

Usol(:,1)=U1;
Uold=U1;


for i=2:Nt
    
    Uold(1)=LBC(i);Uold(end)=UBC(i);
    Unew=A\Uold;
    
    Usol(:,i)=Unew;
    Uold=Unew;
    
    
    
end
plot(exp(Xi),Usol(:,end),'K')

Price=spline(exp(Xi),Usol(:,end),S0);
end