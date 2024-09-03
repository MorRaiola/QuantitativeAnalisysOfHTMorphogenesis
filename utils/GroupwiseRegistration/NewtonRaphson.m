%% Newton-Raphson scheme for updating (nu): Degrees of freedom for TMM
% Written by Nishant Ravikumar 02/09/2015, The University of Sheffield
function nu = NewtonRaphson(nu_i,sPU,D)

tol = 1e-2;
iter=1;
nu_is=nu_i+tol;
while iter<100
  %  disp(['Iteration = ' num2str(iter)]);
    funX = -psi((nu_is)/2) + log((nu_is)/2)+ 1 + sPU + psi((nu_i+D)*0.5) - log((nu_i+D)*0.5);
    der1 = -psi(1,nu_is/2)/2 + 1/nu_is;
    step = funX/der1;
    nu_i=nu_is;
    nu_is = nu_i - 0.01*step;
    iter=iter+1;
end
nu = nu_is;