% Part of the code package related to the publication:
%
% H. O. Caldag, S. Yesilyurt, Steering Control of Magnetic Helical Swimmers
% in Swirling Flows due to Confinement, 2020 International Conference on
% Robotics and Automation. https://doi.org/10.1109/ICRA40945.2020.9196521
%
% swirl_compute is a function called from helix_control_with_swirl.m. The
% script evaluates the forces acting on the helix due to the swirling flow.
% The force is evaluated for small parts of helix and then added up.


function [F_flow]=swirl_compute(B, lam, a, ct, n, curpos, R,Gamma)

kw = 2*pi/lam; % Wave number
cn = 2*ct; % Normal drag coefficient
N=1e3; % Number of sections the helix is divided to for force evaluation
ds=1/N; % Infinitesimal piece size
s=0:ds:n; % Sections in the helix array

cks = cos(kw*s*a); 
sks = sin(kw*s*a); % Some shorthand notation

x  = B*cks; 
y  = B*sks;
z  = s*a; % Helix geometry in terms of s

dx = -B*kw*a*sks;
dy =  B*kw*a*cks; 
dz =  a; % Small pieces of the helix in terms of s

tmag = sqrt(a*a+dx.*dx+dy.*dy); % Local tangent
tx = dx./tmag; 
ty = dy./tmag;
tz = dz./tmag; % Components of the tangent
txyz=[tx;ty;tz]; % All together
omegaf=Gamma*2*pi; % Swirl flow amplitude
KU=[0 0 0; 0 0 -omegaf; 0 omegaf 0]; 
% Swirl components are placed in a matrix for the computations, see Eq. (24) 

scount=0; % Section counter
for s=0:ds:n*lam % As we sweep the sections
    scount=scount+1; % Increase the section counter
    spos=[x(scount);y(scount);z(scount)]; % Define the center coordinate of the section
    uf_body=R'*(KU*(curpos+R*spos)); % Express the swirling flow in body frame
    uf_body_t=(uf_body'*txyz(:,scount)).*txyz(:,scount); % Extract the tangential component
    uf_body_n=(uf_body-uf_body_t); % Extract the normal component
    
    dF_t=ct*uf_body_t*ds; % Evaluate the tangential and normal forces
    dF_n=cn*uf_body_n*ds;

    dF(:,scount)=dF_t+dF_n; % Add them up for the total force on a single section
end

F_flow=sum(dF'); % Add up the forces on the sections to find the total force
F_flow=F_flow';

end
