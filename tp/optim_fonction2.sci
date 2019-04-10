//MAM_OPT_TP8
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



function Av1 = Av(v)
    n = length(v)
//    Av1 = 2*v // cost = J1
//    Av1 = Av3(v) // cost = J3
    Av1 = Av5(v) // cost = J4
endfunction

function Av1 = Av3(v)
    n = length(v)
    Av1(1) = 2*v(1) - v(2)
    for i = 2:n-1
        Av1(i) = -v(i-1) + 2*v(i) - v(i+1)
    end
    Av1(n) = -v(n-1) + 2*v(n)
endfunction


function [J,G]=cost3(v)
    n = length(v)
    f = ones(n,1)
    J = 0.5*Av3(v)'*v - f'*v
    G = Av3(v) - f

endfunction


function Av1 = Av5(v)
    n = length(v)
    Av1(1) = 4*v(1) - v(2) - v(3)
    Av1(2) = -v(1) + 4*v(2) - v(3) - v(4)
    for i = 3:n-2
        Av1(i) = -v(i-2) - v(i-1) + 4*v(i) - v(i+1) - v(i+2)
    end
    Av1(n-1) = -v(n-3) - v(n-2) + 4*v(n-1) - v(n)
    Av1(n) = -v(n-2) - v(n-1) + 4*v(n)
endfunction

function [J,G]=cost4(v)
    n = length(v)
    f = ones(n,1)
    J = 0.5*Av5(v)'*v - f'*v
    G = Av5(v) - f
endfunction
