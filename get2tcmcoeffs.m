function [A] = get2tcmcoeffs(k_1,k_2,k_3,k_4)

deltat=sqrt( (k_2+k_3+k_4)^2 -4*k_2*k_4)

theta_1 = (k_2+k_3+k_4+deltat)/2;
theta_2 = (k_2+k_3+k_4-deltat)/2;

phi_1=k_1*(theta_1 - k_3 - k_4)/deltat;
phi_2=-1*k_1*(theta_2 - k_3 - k_4)/deltat;

A = [phi_1,theta_1,phi_2,theta_2];

end

