% Part of the code package related to the publication:
%
% H. O. Caldag, S. Yesilyurt, Steering Control of Magnetic Helical Swimmers
% in Swirling Flows due to Confinement, 2020 International Conference on
% Robotics and Automation. https://doi.org/10.1109/ICRA40945.2020.9196521
%
% This code simulates a slender helix with a swirling flow and applies
% steering control to tilt the helix towards the reference position. The
% reference position is (0,0) on the y-z plane and the swimmer propels in
% the x- direction. Swirling flow is added to replicate the flow field
% generated in confinement that leads to circular radial trajectories.
%
% Parameters of interest to the end user would be the geometric parameters
% (defined right below), control gain and swirl flow amplitudes (lines 41-42)

clearvars; close all;
set(groot,'defaultLineLineWidth',1.5);
set(groot,'DefaultAxesXGrid','on');
set(groot,'DefaultAxesYGrid','on');
visc = 1; % Non-dimensional viscosity
nlam = 3; % Number of rotations of the helix
lam_helix = 1; % Axial wavelength of the helix
B_helix   = 1; % Helix amplitude
htheta = atan((2*pi*B_helix)/lam_helix); % Helix angle
LAM_helix = lam_helix/cos(htheta); % Wavelength calculated along the helix centerline
bdim=0.0023; % Minor radius of the helix

% Non-dimensional geometric parameters
LAM = 1; % LAM_helix is the length scale
eta = B_helix/LAM_helix; % Non-dimensional helix amplitude
lam = lam_helix/LAM_helix; % Non-dimensional helix wavelength
b = bdim/ LAM_helix; % Non-dimensional minor radius of the helix
kw = 2*pi/lam; % Wave number
a = 1/sqrt(1+kw*kw*eta*eta); % See Eq. (28) for these definitions
phi   = 2*pi*nlam; % Number of rotations of the helix (in radians)

ct = 2*pi*visc/(log(2*lam/b)-1/2); % Tangential drag coefficient
cn = 2*ct; % Normal drag coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kp1 = 5; % Proportional control gain. Set to zero if you don't want any control input.
Gamma=pi/10; % Swirl flow amplitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B0 = 1; % Magnetic field amplitude
Vth = @(r) Gamma*r; % Swirling flow definition
U = 0; % External fluid flow amplitude (unused)
omegax = 2*pi; % Rotation frequency of the helix

Minv=get_mobility_icra(lam,nlam,eta,ct); % Resistance matrix evaluation
% This code piece contains the terms in Eqs. (29)-(32) in the article.

y_0 = 1;
z_0 = 1; % Initial radial position components

R0 = [0 0 1;...
      1 0 0;...
      0 1 0]; % Initial orientation of the particle
x0 = [0;y_0;z_0;R0(:,1);R0(:,2);R0(:,3)]; % Matrix containing the initial 
% position and orientation components

tf = 50; % Simulation duration
[t,x] = ode45(@(t,x) fun(t,x,U,Vth,omegax,B0,Kp1,Minv,B_helix, b, lam, a, ct, nlam,Gamma),[0,tf],x0);
% Main simulation command ^
% This outputs the position and orientation of the swimmer at times t.

% The additional run is to extract additional data from the simulation:
% [tmp,out] = fun(t,x,U,Vth,omegax,B0,Kp1,Minv,B_helix, b, lam, a, ct, nlam,Gamma);

% We utilize the ode solver to quickly solve the problem with variable time
% steps. Then we extract additional data on solved time points only with the 
% second call to fun. It is commented out by default.

%% Plot generation

figure(1) % displays the trajectory
subplot(411)
plot(t,x(:,1)) % Axial position
ylabel('x')
subplot(412)
plot(t,x(:,2)) % y- coordinate
ylabel('y')
subplot(413)
plot(t,x(:,3)) % z- coordinate
ylabel('z')
subplot(414)
plot(t,sqrt(x(:,2).^2+x(:,3).^2)) % Radial position
ylabel('r')
xlabel('t')

figure(2) % 3-dimensional and radial trajectory
subplot(121) % 3-dimensional trajectory
plot3(x(:,1),x(:,2),x(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
subplot(122) % Radial trajectory
plot(x(:,2),x(:,3))
xlabel('y');ylabel('z')
set(gcf,'Position',[680 503 806 375])

%% Function fun evaluates the trajectory of the helix.

% The function can be called in two ways. If it is called with a solver
% (ode45 used above), it solves the problem with variable time stepping.
% After solving with a solver, the function can be called on its own to
% extract additional data at the solved time steps. The first call is
% required to be able to make the second call.

function [dx,out] = fun(t,x,U,Vth,omegax,B0,Kp1,Minv,B_helix, b, lam, a, ct, n,Gamma)
eps = 1e-6; % Swirling flow will not have a component at very low radial positions

if nargout > 1 % If multiple outputs will be returned, we will do some additional evaluations
    r = sqrt(x(:,2).*x(:,2)+x(:,3).*x(:,3)); % Evaluate radial position
    th = atan2(x(:,3),x(:,2)); % Radial orientation angle
    
    % Swirling flow components
    Uy = -sin(th).*Vth(r).*(r>eps); 
    Uz =  cos(th).*Vth(r).*(r>eps);

    % Magnetic field components. Magnetic field rotates on y-z plane
    By = cos(omegax*t);
    Bz = sin(omegax*t);
    Bx = -Kp1*(By.*x(:,3) - Bz.*x(:,2));
    % Bx is the magnetic field in the propulsion direction and main input
    % for position control

    B = B0*[Bx';By';Bz']; % Complete expression of the magnetic field
    
    for k=1:size(t) % Going through time steps
        curpos=x(k,1:3)'; % Current position
       
        R(:,1) = x(k,4:6)';
        R(:,2) = x(k,7:9)';
        R(:,3) = x(k,10:12)'; % Rotation matrix components
        
        [F_flow]=swirl_compute(B_helix,lam, a, ct, n,curpos,R,Gamma);
        % This function evaluates the forces due to swirling flows
        
        Tm = cross([0;1;0],R'*B(:,k)); % Magnetic torque computation
        
        uw =  Minv*([zeros(3,1);Tm]-[F_flow;zeros(3,1)]); % Instantaneous velocity computation
        % Velocities = Mobility matrix * [Forces; Torques]

        out.u(k,:)  = (R*uw(1:3) + [0;Uy(k);Uz(k)])'; % Velocities with swirling flow added
        out.u2(k,:) = (R*uw(1:3))'; % Velocities without the swirling flow
        out.w(k,:)  = (R*uw(4:6))'; % Angular velocities
        out.eb(k)=norm(R'*[0;x(k,2);x(k,3)]); % Orientation vector
    end
    out.Vdot_y = x(:,2).*out.u(:,2);
    out.Vdot_z = x(:,3).*out.u(:,3);
    out.Vdot = out.Vdot_y + out.Vdot_z;
    out.Vdot2 = x(:,2).*out.u2(:,2) + x(:,3).*out.u2(:,3);
    out.r    = r; % Radial position
    out.th   = th; % Polar position
    dx = zeros(size(t)); % Reset at the end of simulation run
    return
end

% Below is what happens if ode45 solver is called:

RR = zeros(3);
dx = zeros(12,1);
for k = 1:3
    RR(:,k) = [x(k*3+1);x(k*3+2);x(k*3+3)];
end % Setting the initial orientation vector
[RL,RS,RV] = svd(RR);
R = RL*RV'; % Deriving the corresponding rotation matrix

r = sqrt(x(2)*x(2)+x(3)*x(3)); % Radial position
th = atan2(x(3),x(2)); % Polar position

curpos=x(1:3); % Current position

[F_flow]=swirl_compute(B_helix,lam, a, ct, n,curpos,R,Gamma);
% Compute the forcing by the swirling flow

By = sin(omegax*t); % Magnetic field components
Bz = cos(omegax*t);
Bx = -Kp1*(By*x(3) - Bz*x(2)); % Control input
B = B0*[Bx;By;Bz]; % Complete magnetic field

Tm = cross([0;1;0],R'*B);  % Magnetic torque in the swimmer frame

uw =  Minv*([zeros(3,1);Tm]-[F_flow; zeros(3,1)]); % Instantaneous velocity computation
% Velocities = Mobility matrix * [Forces; Torques]

dx(1:3,1) = R*uw(1:3);  % Linear velocities
w  = R*uw(4:6); % Angular velocities

for k=1:3
    dx(k*3+1:k*3+3,1) = cross(w,x(k*3+1:k*3+3));
end %dR/dt = w x R where R is the rotation matrix

disp([t,x(1:3)' (R*Tm)']) % Output message
end