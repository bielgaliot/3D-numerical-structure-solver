function Kel = computeKelBar(n_d,nel,x,Tnod,mat,Tmat)

R=zeros(2*(n_d+1));

for e=1:nel %sacar long barra i llenar K_p
    
    if n_d==1
    x1=x(Tnod(e,1),1);
    x2=x(Tnod(e,2),1);
    l=x2-x1;
    R=(1/l)*[x2-x1 0 0 0; 
           0 l 0 0;
           0 0 x2-x1 0;
           0 0 0 l];
       
    K_p =  (mat(Tmat(e),3)*mat(Tmat(e),1))/l^3*[12 6*l -12 6*l; 
                                            6*l 4*l^2 -6*l 2*l^2; 
                                            -12 -6*l 12 -6*l; 
                                            6*l 2*l^2 -6*l 4*l^2];
    Kel(:,:,e)=R.'*K_p*R;
    end
    
    if n_d==2 %%%%ADAPTAR PARA 2D
    y1=x(Tn(e,1),2);
    y2=x(Tn(e,2),2);
    l=sqrt((x2-x1)^2+(y2-y1)^2);
    R=1/l*[x2-x1 y2-y1 0 0; 0 0 x2-x1 y2-y1];
    K_p=mat(Tmat(e),2)*mat(Tmat(e),1)/l*[1 -1; -1 1];
    end
            
    

end


%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------

% It must provide as output:
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------



end