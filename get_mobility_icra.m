% Part of the code package related to the publication:
%
% H. O. Caldag & S. Yesilyurt (2020). Steering control of magnetic helical 
% swimmers in swirling flows due to confinement. ICRA 2020,
% https://doi.org/10.1109/ICRA40945.2020.9196521
%
% This function evaluates the mobility matrix for a slender helix.
% The formulation is a simplified version of the expressions in
%
% Y. Man, E. Lauga (2013). The wobbling-to-swimming transition of rotated 
% helices. Physics of Fluids, https://doi.org/10.1063/1.4812637.
%
% Inputs to the function:
%
% lam: Wavelength of the helix
% nlam: Number of rotations
% B: Major radius
% ct: Tangential drag coefficient, i.e. ct=2*pi*visc/(log(2*lam/b)-1/2)
%     where visc is the dynamic viscosity, lam is the wavelength and b is
%     the minor radius of the helix.
%
% The function outputs D1, the complete mobility matrix.

function D1 = get_mobility_icra(lam,nlam,B,ct)

kw= 2*pi/lam; C=ct*nlam/(sqrt(1+kw^2*B^2));
amob=1/(sqrt(1+kw^2*B^2));

A1=-C/lam; A2=3*pi^2*B^2; A3=A2+lam^2;
A4=kw^2*B^4 + 2/3*A2*nlam^2-1/4*B^2+2/3*lam^2*nlam^2;
A5=kw^2*B^4+2/3*A2*nlam^2-5/4*B^2+2/3*lam^2*nlam^2;

RL=A1*[2*A3 0 0; 0 2*A3 0; 0 0 8/3*A2+lam^2];
RC=A1*[lam*A2/(2*pi) nlam*lam*A3 0; -nlam*lam*A3 lam*A2/(6*pi) 0; 0 -lam^2*B -2*pi*lam*B^2];

RR=A1*lam^2*[A4 nlam*A2/6 0;...
    nlam*A2/6 A5 B/(pi*lam)*(2/3*A2+lam^2);...
    0 B/(pi*lam)*(2/3*A2+lam^2) 2*B/lam*(2/3*A2+lam^2)];

M=[RL RC; RC' RR];
D1=inv(M);

end