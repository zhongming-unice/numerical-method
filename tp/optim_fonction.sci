//MAM_OPT_TP7
//auteur : ZHONG Ming
//
clear()

function [J,G]=cost1(v)
    J = 0
    n = length(v)
    for i=1:n
        J = J + (v(i)-1)^2
    end
    G = 2*v-2
endfunction



function [J,G]=cost2(v)
    J = 0
    n = length(v)
    G = zeros(n)
    for i=1:n
        J = J + (v(i)-i)^2
        G(i) = 2*(v(i)-i)
    end

endfunction




function [J,G]=costR(v)
    J = 0
    n = length(v)
    G = zeros(n-1)
    for i=1:n-1
        J = J + (v(i+1)-v(i)^2)^2 + (v(i)-i)^2
        G(i) = -4*v(i)*(v(i+1)-v(i)^2)+2*(v(i)-1)
    end
    J = J + (v(n)-1)^2
    G(n) = 2*(v(n)-i)
endfunction    





function [J,G]=costH(v)
    x = v(1)
    y = v(2)
    J = (x^2+y-2)^2+(y^2-2*x+1)^2
    G(1) = 4*x*(x^2+y-2)+4*(y^2-2*x+1)
    G(2) = 2*(x^2+y-2)+4*y*(y^2-2*x+1)  
endfunction
