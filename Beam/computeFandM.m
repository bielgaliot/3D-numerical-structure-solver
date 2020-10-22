function [Fy, Mz] = computeFandM(n_d,nel,x,Tnod, Td, n_ne, n_i, u, Kel, R)

Fy=zeros(nel, n_ne);
Mz=zeros(nel, n_ne);
u_e=zeros(n_ne*n_i,1);

for e=1:nel %sacar long barra i llenar K_p
    
    if n_d==1
    x1=x(Tnod(e,1),1);
    x2=x(Tnod(e,2),1);
    l=x2-x1;
    R=1/l*[x2-x1 0 0 0; 
           0 l 0 0;
           0 0 x2-x1 0;
           0 0 0 l];
    end
    
    if n_d==2
        
       A=0; %edit if necessary 
       
    end
       
    for i=1:n_ne*n_i
    
    I=Td(e,i);
    u_e(i,1)=u(I,1);
  
    end
    
    F_int_el = Kel(:,:,e)*u_e;
    F_int_el_p = R*F_int_el;
    
    
    Fy(e,1) = -F_int_el_p(n_i-1);
    Fy(e,2) = F_int_el_p(2*n_i-1);
    
    Mz(e,1) = -F_int_el_p(n_i);
    Mz(e,2) = F_int_el_p(2*n_i);
    
end






end