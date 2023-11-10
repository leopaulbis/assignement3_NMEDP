%PREPROCES
%Reference element
close all;
clear
degree = 1, typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);

nOfElementNodes = size(theReferenceElement.N,2);
figure(1), drawReferenceElement(theReferenceElement);

fileName=mesh_pyzo(5,'Q',degree,0); %creation of the associated mesh 
load(fileName);


figure(2), clf
PlotMesh(T,X,typeOfElement,'k-',1);

T_boun_h1 = Tb_h1(:,[1 end]);
X_boun_h1 = X(T_boun_h1(:),:);

T_boun_h2 = Tb_h2(:,[1 end]);
X_boun_h2 = X(T_boun_h2(:),:);

nodesCCD=find(ismember(X,X_boun_h1,'rows')|ismember(X,X_boun_h2,'rows'));%find the nodes on the dirichlet boundary
uCCD=DirichletValue2(X(nodesCCD,:));
% Definition of Dirichlet boundary conditions
% x = X(:,1); y = X(:,2); tol=1.e-10;
% nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary
% hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
% uCCD=DirichletValue(X(nodesCCD,:));

%System of equations (without BC)
[K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);
% 
% %Boundary conditions & SYSTEM SOLUTION
methodDBC = 2; % 1=System reduction or 2=Lagrange multipliers 
if methodDBC == 1 %System reduction
    unknowns= setdiff([1:size(X,1)],nodesCCD);
    %disp(unknowns)
    f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
    K=K(unknowns,unknowns);
    %System solution
    sol=K\f;
    %Nodal values: system solution and Dirichlet values
    u = zeros(size(X,1),1);
    u(unknowns) = sol; u(nodesCCD) = uCCD;
else %LagrangeMultipliers
    nOfDirichletDOF = length(uCCD); 
   
    nOfDOF = size(K,1);
    
    A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
    A(:,nodesCCD) = eye(nOfDirichletDOF);
    b = uCCD;
    Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
    ftot = [f;b];
    sol = Ktot\ftot;
    u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
end



%POSTPROCESS
figure(3)
PlotNodalField(u,X,T), title('FEM solution')
[Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement);
 %figure(4)
% PlotMesh(T,X,1,'k-'); hold on, quiver(Xpg(:,1),Xpg(:,2),Fgp(:,1),Fgp(:,1),'LineWidth',2), hold off
figure(5)

% Suppose you have a column vector u containing the values of u
% and a matrix X containing the coordinates (x, y)

% Extract x and y coordinates from the matrix X
x = X(:, 1);  % Assuming x-coordinates are in the first column
y = X(:, 2);  % Assuming y-coordinates are in the second column

% Create a scatter plot
scatter(x, y, [], u, 'filled');  % 'filled' option fills the markers based on u values

% Add labels and title
xlabel('x');
ylabel('y');
title('Piezometric Level u(x, y) Distribution');

% Add a color bar
colorbar,
% 
