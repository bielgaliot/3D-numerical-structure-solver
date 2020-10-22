function [pu, pt] = computepupt(n_d, nel, x, Tnod, Td, n_ne, n_i, u, R)


coefs=zeros(4,1);
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
    
    u_e_p=R*u_e;
    
    coefs = 1/(l^3)*[2 l -2 l;
                    -3*l -2*l^2 3*l -l^2;
                    0 l^3 0 0;
                    l^3 0 0 0;
                    ]*u_e;

            

    
    pu(e,[1,2,3,4])=[coefs(1,1), coefs(2,1), coefs(3,1), coefs(4,1)];
    
    pt(e,[1,2,3])=[3*coefs(1,1), 2*coefs(2,1), coefs(3,1)];
    
    
    
end



end
    