function[size_mesh,lambda_val]=mesh_influence_lambda(i,nOfElementNodes,theReferenceElement,degree,bool) %i=mesh fine max num_dep=nombre de neouds pris au départ
  
    lambda_val=zeros(1,i);
    size_mesh=zeros(1,i);

    for k=1:i

        fileName=mesh_pyzo(10,'Q',degree,k-1); %creation of the associated mesh 
        load(fileName);

        size_mesh(1,k)=k-1

        T_boun_h1 = Tb_h1(:,[1 end]);
        X_boun_h1 = X(T_boun_h1(:),:);

        T_boun_h2 = Tb_h2(:,[1 end]);
        X_boun_h2 = X(T_boun_h2(:),:);
        
        if bool

            T_boun = Tb_artificial(:,[1 end]);
            X_boun = X(T_boun(:),:);

            nodesCCD=find(ismember(X,X_boun_h1,'rows')|ismember(X,X_boun_h2,'rows')|ismember(X,X_boun,'rows'));%find the nodes on the dirichlet boundary
        else
            nodesCCD=find(ismember(X,X_boun_h1,'rows')|ismember(X,X_boun_h2,'rows'));
        end
        uCCD=DirichletValue2(X(nodesCCD,:));

        [K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);

        nOfDirichletDOF = length(uCCD); nOfDOF = size(K,1);
        A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
        A(:,nodesCCD) = eye(nOfDirichletDOF);
        b = uCCD;
        Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
        ftot = [f;b];
        sol = Ktot\ftot;
        u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
        
        nodesCCD2=find(ismember(X,X_boun_h2,'rows'));

        lambda2=lambda(1:numel(nodesCCD2));
        lambda_val(1,k)=sum(lambda2);

       
    end
end