
function y=analytic_deriv_x(X)
y=(1/2000)*h(X(:,2)).*v(X(:,1));

function y=v(x)
y=(2*x+4*x.^2-16*x.^3+10*x.^4).*exp(10*x);

function y=h(x)
y = (x.^2).*((1-x).^2).*exp(10*x);
