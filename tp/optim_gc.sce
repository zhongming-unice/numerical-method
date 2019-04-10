//MAM_OPT_TP8
//auteur : ZHONG Ming
//
//exec('optim_fonction2.sci')

function [J,G]=cost(v)
    [J,G]=cost4(v)
endfunction


N = 40
epsg = 10^-6
Kmax = 100
uk = rand(N,1)
J_cost_gf = zeros(Kmax,1)
J_cost_gc = zeros(Kmax,1)
err_gf = zeros(Kmax,1)
err_gc = zeros(Kmax,1)

u_exact = zeros(N,1)
for i = 1:N
    u_exact(i) = 0.5*i*(N+1-i)
end
[J,d] = cost(uk)
for k = 1:Kmax
    [J,G] = cost(uk)
//    d = G
    J_cost_gc(k) = J
    if(norm(G)<epsg) then break
    else
        rho = (G'*d)/(Av(d)'*d)
        u1 = uk - rho*G
        [J1,G1] = cost(u1)
        betak = -(G1'*Av(d))/(d'*Av(d))
        d1 = G1 + betak*d
        d = d1
        uk = u1
    end
    err_gc(k) = norm(uk-u_exact)
end

pas_fixe = 0.1
u = rand(N,1)

for k = 1:Kmax
    [J,G] = cost(u)
    J_cost_gf(k) = J
    if(norm(G)<epsg) then break
    else
        u = u - pas_fixe*G
    end
    err_gf(k) = norm(u-u_exact)
end

scf(1)
plot2d(J_cost_gf,style=11)
plot2d(J_cost_gc,style=22)
legend("pas fixe","gradient conjugue")

scf(2)
plot2d(err_gf,style=33)
plot2d(err_gc,style=5)
legend("err pas fixe","err gradient conjugue")
