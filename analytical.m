function u=analytical(X)

u = (1/2000)*h(X(:,1)).*h(X(:,2));


function y=h(x)
y = (x.^2).*((1-x).^2).*exp(10*x);



