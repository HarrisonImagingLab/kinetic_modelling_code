function [k_1, k_2, k_3,k_4 V_T] = get2tcmparams(A)

%%A should be the vector of coefficients, A(1) and A(3) are amplitudes, A(3) and A(4) are parameters!

[theta_2, n] = min([A(2),A(4)]);
[theta_1, m]=max([A(2),A(4)]);

amps=[A(1), A(3)];

phi_1=amps(m);
phi_2=amps(n);

delta = theta_1 - theta_2;
theta_bar = theta_1 + theta_2;
psi = (theta_1 + (theta_2*phi_1/phi_2))/((phi_1/phi_2) +1 );

k_1 = phi_1* delta *(1/(theta_1-psi))
k_2 = 2*theta_1-delta-psi
k_4 = (delta^2 - (psi + k_2)^2)/(-4*k_2)
k_3= psi-k_4

V_T = (k_1/k_2)*(1 + (k_3/k_4))

end
