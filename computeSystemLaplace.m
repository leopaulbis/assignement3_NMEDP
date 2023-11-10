function [K,f]=computeSystemLaplace(X,T,theReferenceElement,sourceTermFunction)

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;%function evaluated on the gauss point
Nxi=theReferenceElement.Nxi; %dérivative with respect to xi on the gauss point
Neta=theReferenceElement.Neta;%same with respect to eta

%number of nodes =number of row of the matrix which contains the coordinate of the nodes
nOfNodes = size(X,1); 
nOfElements = size(T,1);%number of element=nb of row in the connectivity matrix

%initialize a sparse matrix of size n0nodes*n0nodes whith room for up to 9*nonodes 
K=spalloc(nOfNodes,nOfNodes,9*nOfNodes);
f=zeros(nOfNodes,1);%initialise the vector of the second member 

%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element
    [Ke,fe]=elementalComputations(Xe,IPcoord,IPweights,N,Nxi,Neta,sourceTermFunction);%somme des intégrales sur omega e
    K(Te,Te)=K(Te,Te)+Ke; %assembly of elemental matrix
    f(Te) = f(Te) + fe;
    %figure(11), spy(K), disp('Press any key to continue'), pause
end

%_______________________________________
%Calcul de la matriu i vector elementals
function [Ke,fe]=elementalComputations(Xe,IPcoord,IPweights,N,Nxi,Neta,sourceTermFunction)

nnodes = size(Xe,1);
Ke=zeros(nnodes);
fe=zeros(nnodes,1);
xe = Xe(:,1); ye = Xe(:,2);
%Bucle en punts d integraci
for k=1:length(IPweights)
    Nk=N(k,:);
    Nkxi=Nxi(k,:);
    Nketa=Neta(k,:); 
    xk = Nk*Xe; %map the gauss point to the physical element=isoparametric transformation
    %Jacobia 
    J = [Nkxi*xe Nkxi*ye;Nketa*xe Nketa*ye];%jacobienne du changement isoparametric 
    % Derivadas de las funciones de forma respecto a (x,y)
    Nkxy = J\[Nkxi;Nketa];
    Nkx = Nkxy(1,:); Nky = Nkxy(2,:);
    %diferencial de volum
    dxy=IPweights(k)*det(J);
    Gk=Nkxy;
    Ke = Ke + Gk'*Gk*dxy;%fait le produit des deux gros trucs dans l'intégrale 
    fe = fe + sourceTermFunction(xk)*Nk'*dxy;
end
  