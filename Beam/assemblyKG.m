function [KG,Fext] = assemblyKG(n_ne,nel,n_i,n_dof,Td,Kel,Fe)

KG=zeros(n_dof);
Fext = zeros(n_dof,1);

for e=1:nel
    for i=1:n_ne*n_i
        I=Td(e,i);
        Fext(I) = Fext(I) + Fe(i,e);
        for j=1:n_ne*n_i
            J=Td(e,j);
            KG(I,J)=KG(I,J)+Kel(i,j,e);
        end
    end
end



%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_el       Total number of elements
%                  n_el_dof   Number of DOFs per element 
%                  n_dof      Total number of DOFs
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------
% It must provide as output:

%   - KG    Global stiffness matrix [n_dof x n_dof]
%            KG(I,J) - Term in (I,J) position of global stiffness matrix
%--------------------------------------------------------------------------
end