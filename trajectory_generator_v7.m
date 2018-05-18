function [desired_state] = trajectory_generator(t, qn, varargin)
% TRAJECTORY_GENERATOR: Turn a Dijkstra or A* path into a trajectory
% t: A scalar, specifying inquiry time
%
% varargin: variable number of input arguments. In the framework,
% this function will first (and only once!) be called like this:
%
% trajectory_generator([],[], 0, path)
%
% i.e. map = varargin{1} and path = varargin{2}.
%
% path: A N x 3 matrix where each row is (x, y, z) coordinate of a
% point in the path. N is the total number of points in the path
%
% This is when you compute and store the trajectory.
%
% Later it will be called with only t and qn as an argument, at
% which point you generate the desired state for point t.
%
    
desired_state = [];
Vmax = 2.5;
Vmin = 0.3;
dmin = 0.2;
dmax = 4;

% use the "persistent" keyword to keep your trajectory around
% inbetween function calls

persistent traj Ts L Ldot Ldotdot

%
% When called without varargin (isempty(varargin) == true), compute
% and return the desired state here.
%
if isempty(varargin)
    i = find(Ts<=t, 1, 'last');
    if isempty(i)
        prev = traj(1,:);
        vec = zeros(1,3);
        t_nd = 0;
        T = 1;
    elseif i == length(Ts)
        prev = traj(end,:);
        vec = zeros(1,3);
        t_nd = 0;
        T = 1;
    else
        prev = traj(i,:);
        vec = traj(i+1,:)-traj(i,:);
        tf = Ts(i);
        T = Ts(i+1)-Ts(i);
        t_nd = (t-tf)/T;
    end
    desired_state.pos = prev' + vec'*L(t_nd);
    desired_state.vel = vec'*Ldot(t_nd)/T;
    desired_state.acc = vec'*Ldotdot(t_nd)/T^2;
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
else
    path = varargin{2};
    
    N = size(path,1);
    % create traj and length
    traj = path;
    d = zeros(N,1);
    for i = 2:N
        d(i) = d(i-1) + norm(traj(i,:)-traj(i-1,:));
    end
    
    diffd = diff(d);
    V = zeros(N-1,1);
    for i = 1:(N-1)
        d1 = diffd(i);
        if d1 < dmin
            V(i) = Vmin;
        elseif d1 < dmax
            V(i) = Vmin+(Vmax-Vmin)*(d1-dmin)/(dmax-dmin);
        else
            V(i) = Vmax;
        end
    end
    diffTs = diffd./V;
    Ts = [0; cumsum(diffTs)];
    
    % non-dimensionalize min acceleration formula
    Tf = 1;
    L0 = 0;
    Lf = 1;
    A = [1 0 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 0 1 0;
        -10/Tf^3 10/Tf^3 -6/Tf^2 -4/Tf^2 -3/Tf 1/(2*Tf);
        15/Tf^4 -15/Tf^4 8/Tf^3 7/Tf^3 3/Tf^2 -1/Tf^2;
        -6/Tf^5 6/Tf^5 -3/Tf^4 -3/Tf^4 -1/Tf^3 1/(2*Tf^3)];
    c = A*[L0; Lf; 0; 0; 0; 0];
    L = @(t) c'*[1; t; t^2; t^3; t^4; t^5];
    Ldot = @(t) c'*[0; 1; 2*t; 3*t^2; 4*t^3; 5*t^4];
    Ldotdot = @(t) c'*[0; 0; 2; 6*t; 12*t^2; 20*t^3];
end

end
