%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: 9/03/20
% Author/s: Biel Galiot, Javier Roset.
%

clear;
close all;

%% INPUT DATA

% Material properties
E = 71e9;

% Cross-section parameters
a = 30e-3;
b = 105e-3;
h = 900e-3;
t = 6e-3;

% Other data
g = 9.81;
L1 = 6;
L2 = 12;
Me = 2360;
M=36000;

% Number of elements for each part
% -- Total number of elements = nel_min*nel_part
nel_min = 2; % Modify the value
nel_part = [2,4,8,16,32,64,128];
%% PRECOMPUTATIONS

% Compute section: 
 A = b*(h+t)-2*((h-t)*((b/2)-(a/2)));
 Iy = ((h+t)^3*b)/12 - 1/6*(((b/2-a/2)*(h-t)^3));
% Compute parameter l: obtenerlo igualando fuerzas verticales, no momentos
 l = computeL(L1, L2, M, Me, g);

% Plot analytical solution
fig = plotBeamsInitialize(L1+L2);

% fixed nodes matrix creation

fixNod = [  1 1 0;
            1 2 0];
        
u_convergence=zeros(1,length(nel_part));

% Loop through each of the number of elements
for k = 1:length(nel_part)

    %% PREPROCESS
    
    % Number of elements
    nel = nel_min*nel_part(k);
    % Number of nodes
    nnod = nel+1;
    
    
 
    
    % Nodal coordinates
     %x(a,j) = %coordinate of node a in the dimension j
     
     x1=linspace(0,L1,nel_part(k)+1);
     x2=linspace(L1,L2+L1,nel_part(k)+1);
     x2(1)=[];
     
     x=transpose([x1 x2]);
    % Complete the coordinates
    
    % Nodal connectivities  
    %  Tnod(e,a) = 
    Tnod = zeros(nel,2);
    for e = 1:nel
        Tnod(e,1) = e;
        Tnod(e,2) = e+1;
    end

    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iy;  % Material (1)
    ];

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    Tmat = ones(nel,1);
        
    %% SOLVER
    
    n_d = size(x,2);                 % Number of dimensions
    n_ne = size(Tnod,2);            % Number of nodes for each bar
    n_i = 2;                       % Number of DOFs for each node
    n_dof = n_i*nnod;                  % Total number of degrees of freedom
    
    % Computation of the DOFs connectivities
    Td = connectDOFs(nel,n_ne,n_i,Tnod);
    
    % Computation of element stiffness matrices (PARA 2D ESTÁ AÚN POR
    % ADAPTAR)
    Kel = computeKelBar(n_d,nel,x,Tnod,mat,Tmat);
    
    %distributing force between elements  
    Fe = computeF(nel,n_ne,n_i,n_d,Tnod,x,l,M,L1,L2,g);
    
    %creating global stifness matrix and global force vector   
    [KG,Fext] = assemblyKG(n_ne,nel,n_i,n_dof,Td,Kel,Fe);
    
    %adding engine weight   
    Fext((nel_part(k)+1)*2-1) = Fext(1+nel_part(k)) - Me * g;
    
    total = 0;
    for force=1:size(Fext,1)  %comprobacion suma fuerzas Fext. Total debe ser 0
        if rem(Fext(force),2)==0
            total = total + Fext(force);
        end
    end
    
    
    %Applying boundry conditions
    [vL,vR,uR] = applyCond(n_i,n_dof,fixNod);
    
    %solving the system and finding displacements, rotations and reactions    
    [u,R] = solveSys(vL,vR,uR,KG,Fext);
    
    %save values for convergence plot
    u_convergence(k)=u(length(u)-1);
    
    %computation of the internal forces and bending moment
    [Fy, Mz] = computeFandM(n_d, nel, x, Tnod, Td, n_ne, n_i, u, Kel, R);
    
    %computation of the internal forces and bending moment
    [pu, pt] = computepupt(n_d, nel, x, Tnod, Td, n_ne, n_i, u, R);
    
    %% POSTPROCESS

    % Number of subdivisions and plots
    nsub = 10;
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
    
   
    
end



% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel_part))),'location','northeast');


%% Convergence Plot

for k = 1:length(nel_part)
u_for_plot(k)=abs((100*(u_convergence(k)-u(length(u)-1))/u(length(u)-1)));
end

nel_part(1)=[];
u_for_plot(1)=[];

figure
loglog(nel_part,u_for_plot)
title 'Convergence plot'
xlabel 'Number of divisions per part'
ylabel 'Difference between "analytical" and numerical solution (%)'
%xticklabels({'2','4','8','16','32','64','128'})
grid minor


