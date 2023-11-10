degree = 2; typeOfElement=0; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
[size_mesh1,lambda_val1]=mesh_influence_lambda(4,nOfElementNodes,theReferenceElement,degree,false);
[size_mesh2,lambda_val2]=mesh_influence_lambda(4,nOfElementNodes,theReferenceElement,degree,true);


disp("size_mesh");
disp(size_mesh1);
disp("sum lambda");
disp(lambda_val1);
disp("size_mesh");
disp(size_mesh2);
disp("sum lambda");
disp(lambda_val2);
figure;

% Tracé du premier ensemble de données
plot(size_mesh1, lambda_val1, 'DisplayName', 'without Dirichlet on the artificial boundary');
hold on;

% Tracé du deuxième ensemble de données
plot(size_mesh2, lambda_val2, 'DisplayName', 'with Dirichlet on artificial boundary');

% Affichage de la légende
legend('Location', 'best');
xlabel('mesh size')
ylabel('flux')

hold off;


