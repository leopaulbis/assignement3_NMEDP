function [H1err] = computeH1error(u, X, T, theReferenceElement)
   

    IPweights = theReferenceElement.IPweights;
    IPcoord = theReferenceElement.IPcoord;
    N = theReferenceElement.N;
    Nxi = theReferenceElement.Nxi;
    Neta = theReferenceElement.Neta;

    nOfElements = size(T, 1);
    L2err = 0;
    H1err = 0;

    for i = 1:nOfElements
        Te = T(i, :);
        Xe = X(Te, :);
        xe = Xe(:, 1);
        ye = Xe(:, 2);
        Ue = u(Te);

        for g = 1:length(IPweights)
            N_g = N(g, :);
            Nxi_g = Nxi(g, :);
            Neta_g = Neta(g, :);

            J = [Nxi_g * xe, Nxi_g * ye; Neta_g * xe, Neta_g * ye];

            grad_ref = [Nxi_g; Neta_g];
           
            grad = J \ grad_ref;

            Nx_g = grad(1, :);
            Ny_g = grad(2, :);

            dvolu = IPweights(g) * det(J);
            Xg = N_g * Xe;

            errorL2 = (analytical(Xg) - N_g * Ue) ^ 2;
            L2err = L2err + errorL2 * dvolu;
            
            errorH1 =(analytic_deriv_x(Xg) - Nx_g * Ue)^2 + (analytic_deriv_y(Xg) - Ny_g * Ue)^2;
            H1err = H1err +errorL2*dvolu+ errorH1 * dvolu;
        end
    end

    H1err = sqrt(H1err);
end
