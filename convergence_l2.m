function[node,error]=convergence_l2(i,num_dep,nOfElementNodes,theReferenceElement) %i=puissance de 10 max, num_dep=nombre de neouds pris au départ
    num=num_dep
    error=zeros(1,i)
    node=zeros(1,i)

    for k=1:i
        [X,T] = CreateMesh(0,nOfElementNodes,[0,1,0,1],num,num)
        node(1,k)=1/(num-1)
        num=num+10*k

        x = X(:,1); y = X(:,2); tol=1.e-10;
        nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol); 
        uCCD=DirichletValue(X(nodesCCD,:));
        [K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);

        unknowns= setdiff([1:size(X,1)],nodesCCD);
        f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
        K=K(unknowns,unknowns);%K_uu
        %System solution
        sol=K\f;
        %Nodal values: system solution and Dirichlet values
        u = zeros(size(X,1),1);
        u(unknowns) = sol; u(nodesCCD) = uCCD;

        error(1,k)=computeL2error(u,X,T,theReferenceElement)
    end
end