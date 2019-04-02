clear()
//author : ZHONG Ming
//mopen('C:\Users\kamibukuro233\numeri\tp\MAM4_MNE_08_TP_carre5.txt','r')
[Mat] = read('C:\Users\kamibukuro233\numeri\tp\MAM4_MNE_08_TP_mask10.msh',-1,3);


Nombre_total = iconvert(Mat(1,1:3),2)
NN = Nombre_total(1)
NE = Nombre_total(2)
NA = Nombre_total(3)

x(1:NN) = double(Mat(2:NN+1,1))
y(1:NN) = double(Mat(2:NN+1,2))
refn = iconvert(Mat(2:NN+1,3),1)

// refe(i,j) : ref du jème (j=1,2,3) noeud dans l'élément i
refe = iconvert(Mat(NN+2:NN+NE+1,1:3),2)

disp(NN)
disp(NE)
disp(NA)

//table dans l'élément n
function [Ke,Fe] = table(n)
    y23 = y(refe(n,2)) - y(refe(n,3))
    y31 = y(refe(n,3)) - y(refe(n,1))
    y12 = y(refe(n,1)) - y(refe(n,2))
    x32 = x(refe(n,3)) - x(refe(n,2))
    x13 = x(refe(n,1)) - x(refe(n,3))
    x21 = x(refe(n,2)) - x(refe(n,1))
    
    a_T = (x21*y31-y12*x13)/2.0
    B(1,:) = [y23, y31, y12]
    B(2,:) = [x32, x13, x21]
    Ke = (B'*B)/(4*a_T)
    
    f = ones(3,1)
    Fe = (1/3)*a_T*f
endfunction

function [Kg,Fg] = assemblage()
    Kg = zeros(NN,NN)
    Fg = zeros(NN,NN)
    
    for n = 1:NE
        [Ke, Fe] = table(n)
        for i = 1:3
            for j = 1:3
                Kg(refe(n,i), refe(n,j)) = Kg(refe(n,i), refe(n,j)) + Ke(i,j)
            end
        end
        
        for k = 1:3
            Fg(refe(n,k),1) = Fg(refe(n,k),1) + Fe(k,1)
        end
    end
    
    for i = 1:NN
        if refn(i,1)==1 then
            Kg(i,i) = Kg(i,i)*10^9
        end
    end
endfunction


[K, F] = assemblage()
u = inv(K)*F
z = zeros(NE,5)
for i = 1:NE
z(i,:) = [i,refe(i,1),refe(i,2),refe(i,3),i]
end
fec(x,y,z,u)

