function Fe = computeF(nel,n_ne,n_i,n_d,Tnod,x,l,M,L1,L2,g)

Fe = zeros(n_ne*n_i,nel);

q1 = @(x) (0.85-0.15.*cos(pi.*x./L1)).*l;

q2 = @(y) -((1./L2.^2).*(L1-L2-y).*(L1+L2-y)).*l;

lambda1 = @(z) (3.*M./(2.*L1.^2)).*(L1-z)+M./(4.*(L1+L2));

lambda2 =  (M./(4.*(L1+L2)));



for el=1:nel
    
    if n_d==1
    x1=x(Tnod(el,1),1);
    x2=x(Tnod(el,2),1);
    l_e=x2-x1;
    R=(1/l_e)*[x2-x1 0 0 0; 
           0 l_e 0 0;
           0 0 x2-x1 0;
           0 0 0 l_e];
    end
    %calculo fuerza del segmento
    
    if el < nel/2
        
        Section_Lift = integral(q1,x1,x2);

        Section_Weight = integral(lambda1,x1,x2)*g;
    
    else
        Section_Lift = integral(q2,x1,x2);

        Section_Weight = lambda2*g*(x2-x1);
    end
    
    q_le = Section_Lift-Section_Weight;
    
    Fe_prima = q_le*(1/2)*[1;
                            l_e/6;
                            1;
                            -l_e/6];
    Fe (:,el) = R.'*Fe_prima;
    
    
end


end