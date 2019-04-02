// author : ZHONG Ming
// equation de transport: schema Euler Decentre, Lax_Friedrichs, Lax-Wendroff
clear()

sigma = 0.1
x0 = 0.25
m = 7
c = 1

function y=Uinit(x)
//    y = exp(-(x-x0)**2/sigma**2)

//    y = sin((m*%pi*x)) + sin((m*%pi*x)/3)

            if x < x0 then
                y = 1;
            else
                y = 0;
            end
endfunction

function y=Uexact(x,t)
    y = Uinit(x-c*t)
endfunction

function y=a(t)
    y= Uexact(0,t)
endfunction

function y=b(t)
    y= Uexact(1,t)
endfunction


//la condition CFL : lambda = c * dt / dx
lambda = 0.3

//la discrétisation du problème en espace : 
//Espace total : 
L = 1
//Nombre de points d'espace : 
Nx = 50
//Pas d ' espace : 
dx = L / Nx


//la discrétisation du problème en temps :
//Temps initial : 
T0 = 0
//Nombre de points de temps : 
Nt = 20
//Pas de temps : 
dt = lambda * dx / c
//Temps final : 
Tfin = T0 + dt * Nt

//valeur 
V_dec=zeros(Nx+1,Nt)
V_laxf=zeros(Nx+1,Nt)
V_laxw=zeros(Nx+1,Nt)
x = linspace(0,L,Nx+1)
t = linspace(T0,Tfin,Nt)

//quand t = n, x = 0,1,2,...,Nx
U_dec_n=zeros(Nx+1,1)
U_laxf_n=zeros(Nx+1,1)
U_laxw_n=zeros(Nx+1,1)
//quand t = n+1, x = 0,1,2,...,Nx
U_dec_np1=zeros(Nx+1,1)
U_laxf_n=zeros(Nx+1,1)
U_laxw_n=zeros(Nx+1,1)

for n = 1:Nt

    U_dec_np1(1)=a(t(n))
    U_dec_np1(Nx+1)=b(t(n))
    U_laxf_np1(1)=a(t(n))
    U_laxf_np1(Nx+1)=b(t(n))
    U_laxw_np1(1)=a(t(n))
    U_laxw_np1(Nx+1)=b(t(n))

    for j=2:Nx
        U_dec_np1(j) = c*dt/dx*(U_dec_n(j-1)-U_dec_np1(j)) + U_dec_n(j)
        U_laxf_np1(j) = 0.5*(U_laxf_n(j-1)+U_laxf_n(j+1)) - c*dt/2/dx*(U_laxf_n(j+1)-U_laxf_n(j-1))
        U_laxw_np1(j) = U_laxw_n(j) - c*dt/2/dx*(U_laxw_n(j+1)-U_laxw_n(j-1)) + c^2*dt^2/2/dx^2*(U_laxw_n(j+1)-2*U_laxw_n(j)+U_laxw_n(j-1))
    end

//    if norm(U_dec_np1-U_dec_n) < 0.005
//        then break
//    end

    U_dec_n = U_dec_np1
    U_laxf_n = U_laxf_np1
    U_laxw_n = U_laxw_np1
    V_dec(:,n) = U_dec_n
    V_laxf(:,n) = U_laxf_n
    V_laxw(:,n) = U_laxw_n
end

//valeur exacte
VE=zeros(Nx+1,Nt)

for i=1:Nt
    for j=1:Nx+1
        VE(j,i)=Uexact(x(j),t(n))
    end
end

y = linspace(0,L,Nx+1)
for i=1:Nx+1
    y(i) = Uexact(x(i),0.15)
end
scf(1)
plot2d(x,U_dec_np1,style=1)
plot2d(x,U_laxf_np1,style=2)
plot2d(x,U_laxw_np1,style=3)
plot2d(x,y,style=6)
legend("Decentre","Lax_Friedrichs","Lax-Wendroff","Valeur exacte")

//l'erreur en fonction des pas de temps
for j=1:Nt
    err_dec_t(j)=norm(VE(:,j)-V_dec(:,j),2)
    err_laxf_t(j)=norm(VE(:,j)-V_laxf(:,j),2)
    err_laxw_t(j)=norm(VE(:,j)-V_laxw(:,j),2)
end
scf(3)
xtitle("erreur en fonction des pas de temps")
plot2d(t,err_dec_t,style=11)
plot2d(t,err_laxf_t,style=12)
plot2d(t,err_laxw_t,style=13)
legend("Decentre","Lax_Friedrichs","Lax-Wendroff")

