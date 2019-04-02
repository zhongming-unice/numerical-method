function y = k(x)
    y = 1;
endfunction

function y = f(x)
    y = -(%pi) ^2 * sin(%pi*x) 
endfunction

function [Xi, Ki]=maillage(n,L)
    NN=n+1; //nombre des noeuds
    NE=n;  //nombre des elements
    h=L/n; //longeur d'intervalle
    Xi=zeros(1,NN);  //veucteur des noeuds
    Ki=zeros(2,NE);  //matrice des elements
    for i= 1:NN;
        Ki(1,i)=h*(i-1);
    end
    for j= 1:NE  //ki = (xi,xi+1)
        Ki(1,j)=Xi(1,j);
        Ki(2,j)=Xi(1,j+1);
    end
endfunction


function [Ke,Fe]=table(ie,n,L)
    Ke = zeros(2,2) ;
    Fe = zeros(2,1) ;
    [Xi,Ki] = maillage (n,L) ;
    
    h = L/n ;
    v = zeros(2,2) ;
    v (:,1) = [1,-1] ;
    v (:,2) = [-1,1] ;
    Ke = k(Xi(1,ie))*v/h ; // la matrice de e rigidit elementaire
    Fe (1,1) = 0.5*f(Ki(1,ie))*h ;
    Fe (2,1) = 0.5*f(Ki(2,ie))*h ; // la matrice de seconde membre elementaire
endfunction

function [KG, FG] = assemblage(n,L)
    KG = zeros(n+1, n+1);
    for i = 1 : n
        [Ke,Fe] = table (i,n,L)
        KG(i,i) = KG(i,i) + Ke(1,1) ;
        KG(i,i+1) = KG(i,i+1) + Ke(1,2) ;
        KG(i+1,i) = KG(i+1,i) + Ke(2,1) ;
        KG(i+1,i+1) = KG(i+1,i+1) + Ke (2,2) ;
    end
    KG(1,1) = 1 ;
    KG(n+1, n+1) = 1 ; // la condition limite    
    
    FG = zeros(n+1,1) ;
    for j = 1:n
        [Ke,Fe] = table (j ,n ,L)
        FG(j,1) = FG(j,1) + Fe(1,1);
        FG(j+1,1) = FG(j+1,1) + Fe(2,1);
    end
    FG(1,1) = 0 ;
    FG(n+1,n+1) = 0 ; // la condition limite
endfunction

function [U] = resolution(n , L)
    [K, F] = assemblage(n , L) ;
    U = inv(K)*F;
//    U = linsolve (K,F)
    x = 0:L/n:L;
    plot2d (x,U(x/(L/n)+1)) ;
endfunction

U = resolution(1000 ,10)
//x = 0:0.01:10 
//plot2d (x , cos (%pi * x)-1, style=2) ; 
