function [F, M, trpy, drpy] = controller(qd, t, qn, params)
% CONTROLLER quadrotor controller
% The current states are:
%% Inputs:
%
% qd{qn}: state and desired state information for quadrotor #qn (qn
%         will be = 1 since we are only flying a single robot)
%
%  qd{qn}.pos, qd{qn}.vel   position and velocity
%  qd{qn}.euler = [roll;pitch;yaw]
%  qd{qn}.omega     angular velocity in body frame
% 
%  qd{qn}.pos_des, qd{qn}.vel_des, qd{qn}.acc_des  desired position, velocity, accel
%  qd{qn}.yaw_des, qd{qn}.yawdot_des
%
% t: current time
%    
% qn: quadrotor number, should always be 1
%    
% params: various parameters
%  params.I     moment of inertia
%  params.grav  gravitational constant g (9.8...m/s^2)
%  params.mass  mass of robot
%
%% Outputs:
%
% F: total thrust commanded (sum of forces from all rotors)
% M: total torque commanded
% trpy: thrust, roll, pitch, yaw (attitude you want to command!)
% drpy: time derivative of trpy
%
% Using these current and desired states, you have to compute the desired
% controls u, and from there F and M
%

% =================== Your code goes here ===================
% ...
% ==============================

%% extract state and paramaters
% state
r = qd{qn}.pos;
v = qd{qn}.vel;
euler = qd{qn}.euler;
omega = qd{qn}.omega;

% trajectory
r_t = qd{qn}.pos_des;
v_t = qd{qn}.vel_des;
a_t = qd{qn}.acc_des;
yaw_t = qd{qn}.yaw_des;
yawdot_t = qd{qn}.yawdot_des;

% parameters
I = params.I;
g = params.grav;
m = params.mass;

%% Initialize controls
% position control
% 20, 10, 15000, 1500
kp_xy = 30;
kp_z_xy = 3;
kp = kp_xy * [1; 1; kp_z_xy];

kd_xy = kp_xy/(sqrt(2)*1.05);
kd_z_xy = kp_z_xy;
kd = kd_xy * [1; 1; kd_z_xy];

% euler angle control
kpe = 30000 * ones(3,1);

kde = 3000 * ones(3,1);

%% calculate inputs

a_des = a_t - kd.*(v-v_t) - kp.*(r-r_t);

% Desired roll, pitch and yaw (in rad). In the simulator, those will be *ignored*.
% When you are flying in the lab, they *will* be used (because the platform
% has a built-in attitude controller). Best to fill them in already
% during simulation.

phi_des   = (a_des(1)*sin(yaw_t) - a_des(2)*cos(yaw_t))/g;
theta_des = (a_des(1)*cos(yaw_t) + a_des(2)*sin(yaw_t))/g;
psi_des   = yaw_t;
euler_des = [phi_des; theta_des; psi_des];

omega_des = [0; 0; yawdot_t];

%
u    = zeros(4,1); % control input u, you should fill this in
u(1) = m*(g + a_des(3));
u(2:4) = I*(-kpe.*(euler - euler_des) - kde.*(omega - omega_des));
                  
% Thrust
F    = u(1);       % This should be F = u(1) from the project handout

% Moment
M    = u(2:4);     % note: params.I has the moment of inertia

% =================== Your code ends here ===================

% Output trpy and drpy as in hardware
trpy = [F, phi_des, theta_des, psi_des];
drpy = [0, 0,       0,         0];

end

%
% ------------------------------------------------------------
%    should you decide to write a geometric controller,
%    the following functions should come in handy
%

function m = eulzxy2rotmat(ang)
    phi   = ang(1);
    theta = ang(2);
    psi   = ang(3);
    
    m = [[cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta), -cos(phi)*sin(psi), ...
          cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi)];
         [cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta),  cos(phi)*cos(psi), ...
          sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)];
         [-cos(phi)*sin(theta), sin(phi), cos(phi)*cos(theta)]];
end

function eul = rotmat2eulzxy(R)
    if R(3,2) < 1
        if R(3,2) > -1
            thetaX = asin(R(3,2));
            thetaZ = atan2(-R(1,2), R(2,2));
            thetaY = atan2(-R(3,1), R(3,3));
        else % R(3,2) == -1
            thetaX = -pi/2;
            thetaZ = -atan2(R(1,3),R(1,1));
            thetaY = 0;
        end
    else % R(3,2) == +1
        thetaX = pi/2;
        thetaZ = atan2(R(1,3),R(1,1));
        thetaY = 0;
    end
    eul = [thetaX, thetaY, thetaZ];
end

function w = veemap(R)
    w = [-R(2,3), R(1,3), -R(1,2)];
end
