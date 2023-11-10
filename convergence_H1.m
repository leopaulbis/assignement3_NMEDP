function[node,error]=convergence_H1(i,num_dep,nOfElementNodes,theReferenceElement) %i=puissance de 10 max, num_dep=nombre de neouds pris au départ
    num=num_dep
    error=zeros(1,i)
    node=zeros(1,i)

    for k=1:i
        [X,T] = CreateMesh(0,nOfElementNodes,[0,1,0,1],num,num)
        node(1,k)=1/(num-1)
        num=num+10*k

        x = X(:,1); y = X(:,2); tol=1.e-10;
        nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); %Nodes on the boundary
        hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off
        uCCD=DirichletValue(X(nodesCCD,:));%on rentre les dirichlet conditions

        [K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);% on impose les conditions limites

        unknowns= setdiff([1:size(X,1)],nodesCCD); %ceux qui sont pqs sur le bord
        f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
        K=K(unknowns,unknowns);%K_uu
        %System solution
        sol=K\f;
        %Nodal values: system solution and Dirichlet values
        u = zeros(size(X,1),1);
        u(unknowns) = sol; u(nodesCCD) = uCCD;

        error(1,k)=computeH1error(u,X,T,theReferenceElement)
    end
end