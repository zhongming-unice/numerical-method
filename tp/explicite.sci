
clear() 

//function u=Uexact(x,y,t)
//    u = 7*x^2 + 4*y^3 + 2*sin(t);
//endfunction
//
//function u=U0(x,y)
//    u = 7*x^2 + 4*y^3;
//endfunction
//
//// avec D = 1
//function f=f(x,y,t)
//    f = 2*cos(t) - 14 - 24*y;
//endfunction
//
//// valeur de U(0,0,t)
//function u=a(t)
//    u = 2*sin(t);
//endfunction
//
//// valeur de U(L,H,t)
//function u=b(t)
//    u = 7 + 32 + 2*sin(t);
//endfunction


function u=Uexact(x,y,t)
    u = sin(2*%pi*x) * sin(%pi*y) * exp(-t);
endfunction

function u=U0(x,y)
    u = sin(2*%pi*x) * sin(%pi*y);
endfunction

// avec D = 1
function f=f(x,y,t)
    f = -sin(2*%pi*x) * sin(%pi*y)*exp(-t) + 4*%pi^2*sin(2*%pi*x)* sin(%pi*y) * exp(-t) + %pi^2*sin(%pi*y)*sin(2*%pi*x)*exp(-t);
endfunction

// valeur de U(0,0,t)
function u=a(t)
    u = 0;
endfunction

// valeur de U(L,H,t)
function u=b(t)
    u = 0;
endfunction


D = 1
//la discrétisation du problème en espace : 
//Espace total : 
L = 1
H = 2
//Nombre de points d'espace : 
Nx = 10
Ny = 30
//Pas d ' espace : 
dx = L / Nx
dy = H / Ny

//la discrétisation du problème en temps :
//Temps initial : 
T0 = 0
//Nombre de points de temps : 
Nt = 100
//Pas de temps : 
dt = 0.002
//Temps final : 
Tfin = T0 + dt * Nt

// la stabilité du schéma, la condition CFL : lambda(1 et 2) < 0.5
lambda1 = D * dt / dx / dx;
lambda2 = D * dt / dy / dy;

//valeur V un tenseur de dimension 3 
V = zeros(Nx+1,Ny+1,Nt)
x = linspace(0,L,Nx+1)
y = linspace(0,H,Ny+1)
t = linspace(T0,Tfin,Nt)


//quand t = n, x = 0,1,2,...,Nx , y = 0,1,2,...,Ny
Un=zeros(Nx+1,Ny+1)
//quand t = n+1, x = 0,1,2,...,Nx , y = 0,1,2,...,Ny
Unp1=zeros(Nx+1,Ny+1)

for n = 1:Nt
    Unp1(1,1)=a(t(n))
    Unp1(Nx+1,Ny+1)=b(t(n))
    for i = 2: Nx
        for j = 2:Ny
            Unp1(i,j) = Un(i,j) + D*((Un(i+1,j)-2*Un(i,j)+Un(i-1,j))/dx^2 + (Un(i,j+1)-2*Un(i,j)+Un(i,j-1))/dy^2)+ dt * f(x(i),y(j),t(n))
        end
    end
    //    if norm(Un-Unp1)<10^(-6) break;
    Un = Unp1
    V(:,:,n)=Un
end

//valeur exacte
VE=zeros(Nx+1,Ny+1,Nt)

for n=1:Nt
    for i=1:Nx+1
        for j=1:Ny+1
        VE(i,j,n)=Uexact(x(i),y(j),t(n))
        end
    end
end

//scf(1)
//contour(x,y,Un)
