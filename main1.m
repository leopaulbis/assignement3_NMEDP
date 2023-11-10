%PREPROCES
%Reference element
close all;
clear
degree = 2; typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
figure(1), drawReferenceElement(theReferenceElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
[X,T] = CreateMesh(0,nOfElementNodes,[0,1,0,1],41,41);
figure(2), clf
PlotMesh(T,X,typeOfElement,'k-',1);
%Definition of Dirichlet boundary conditions
x = X(:,1); y = X(:,2); tol=1.e-10;
nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary
hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
uCCD=DirichletValue(X(nodesCCD,:));

%System of equations (without BC)
[K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);

%Boundary conditions & SYSTEM SOLUTION
methodDBC = 1; % 1=System reduction or 2=Lagrange multipliers 
if methodDBC == 1 %System reduction
    unknowns= setdiff([1:size(X,1)],nodesCCD);
    f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
    K=K(unknowns,unknowns);
    %System solution
    sol=K\f;
    %Nodal values: system solution and Dirichlet values
    u = zeros(size(X,1),1);
    u(unknowns) = sol; u(nodesCCD) = uCCD;
else %LagrangeMultipliers
    nOfDirichletDOF = length(uCCD); nOfDOF = size(K,1);
    A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
    A(:,nodesCCD) = eye(nOfDirichletDOF);
    b = uCCD;
    Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
    ftot = [f;b];
    sol = Ktot\ftot;
    u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
end

sum(lambda);
%POSTPROCESS
figure(3)
PlotNodalField(u,X,T), title('FEM solution')
%[Xpg,Fgp] = gradientElementalField(u,X,T,theReferenceElement);
% figure(4)
% PlotMesh(T,X,1,'k-'); hold on, quiver(Xpg(:,1),Xpg(:,2),Fgp(:,1),Fgp(:,1),'LineWidth',2), hold off

%Comparison with analytical solution
figure(5)
[Xfine,Tfine] = CreateMesh(0,nOfElementNodes,[0,1,0,1],41,41);
PlotNodalField(analytical(Xfine),Xfine,Tfine), title('Analytical solution')

L2error=computeL2error(u,X,T,theReferenceElement);
H1error=computeH1error(u,X,T,theReferenceElement);
disp(L2error);

%PLOT OF THE L2 CONVERGENCE

%Reference element
degree = 1; typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
disp('H1');
disp(H1error);

[node,error_lin_L2]=convergence_l2(5,11,nOfElementNodes,theReferenceElement);
[node,error_lin_H1]=convergence_H1(5,11,nOfElementNodes,theReferenceElement);


% %Reference element
degree = 2; typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);

[node,error_quad_L2]=convergence_l2(5,11,nOfElementNodes,theReferenceElement);
[node,error_quad_H1]=convergence_H1(5,11,nOfElementNodes,theReferenceElement);



 figure(6)
% 
% Tracer la première série de données (erreur L2 pour les éléments linéaires)
loglog(node, error_quad_L2, '-o', 'DisplayName', 'Quadratic quadrilateral Element L2 Error');
hold on % Pour superposer le prochain tracé

% Tracer la deuxième série de données (erreur L2 pour les éléments quadratiques)
loglog(node, error_quad_H1, '-s', 'DisplayName', 'Quadratic quadrilateral Element H1 Error');

legend('Location', 'Best') % Afficher la légende
xlabel('log(element size)')
ylabel('log(L2 Error)')
% 
x = node; % Utiliser les mêmes valeurs x que pour les données
y =x.^2; % L2 Error = c * h^2, où c est une constante
loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 2)');

x = node; % Utiliser les mêmes valeurs x que pour les données
y =x.^3; % L2 Error = c * h^2, où c est une constante
loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 3)');
legend();
title('Convergence Plot')
% 
% 
 hold off 
% 
 figure(7)
% 
loglog(node, error_lin_L2, '-o', 'DisplayName', 'Linear quadrilateral Element L2 Error');
hold on
loglog(node, error_lin_H1, '-o', 'DisplayName', 'Linear quadrilateral element H1 Error');
x = node; % Utiliser les mêmes valeurs x que pour les données
y =x.^2; % L2 Error = c * h^2, où c est une constante
loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 2)');
x = node; % Utiliser les mêmes valeurs x que pour les données
y =x; % L2 Error = c * h^2, où c est une constante
loglog(x, y, '--', 'DisplayName', 'Reference Line (slope = 1)');
legend('Location', 'Best') % Afficher la légende
xlabel('log(element size)')
ylabel('log(L2 Error)')
title("Convergence Plot")
legend();






