//MAM_OPT_TP7
//auteur : ZHONG Ming
//
//exec('optim_fonction.sci')

function [J,G]=cost(v)
    [J,G]=costR(v)
endfunction


N = 10
epsg = 10^-6
Kmax = 100
uk = rand(N,1)
J_cost_gf = zeros(Kmax,1)
J_cost_go = zeros(Kmax,1)


pas_fixe = 0.05
for k = 1:Kmax
    [J,G] = cost(uk)
    J_cost_gf(k) = J
//    if(norm(G)<epsg) then break; end
        u1 = uk - pas_fixe*G
        uk = u1
end

//parabolique : f(t) = a0 + a1*t + a2*t^2 = J(uk - t*G(uk))
// f'(t) = a1 + 2*a2*t  :  f'(tk) = 0  ====>   tk = -a1/(2*a2)
//f(0) = a0 = J(uk)
//f'(0) = a1 = -G(uk)*G(uk)
//f(tk-1) = a0 + a1*tk-1 + a2*tk-1^2 = J(uk - tk-1*G(uk))   ===> a2 = (J(uk - tk-1*G(uk)) - a0 - a1*tk-1) / tk-1^2 
//                                                                  = (J(uk - tk-1*G(uk)) - J(uk) - G(uk)*tk-1) / tk-1^2


u = rand(N,1)
kappa = 0.05
t = zeros(Kmax,1)
t(1) = 1
for k = 2:Kmax
    
    [J,G] = cost(u)
    if(norm(G)<epsg) then break
        else
        J_cost_go(k) = J
    
        a0 = J
        a1 = -norm(G)^2
    
        [J1,G1] = cost(u-t(k-1)*G)
        a2 = (J1-a0-a1*t(k-1)) / t(k-1)^2 
        t(k) = -a1/(2*a2)
        t(k) = t(k) * kappa
        u2 = u - t(k)*G
        u = u2
    end
    
end



scf(1)
plot2d(J_cost_gf,style=11)
plot2d(J_cost_go,style=22)
legend("pas fixe","pas optimal")
