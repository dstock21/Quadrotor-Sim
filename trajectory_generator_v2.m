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
Vf = 1.3;
Vs = 0.5;
dmin = 0.7;

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
    
    % number of contractions
    N = size(path,1);
    % pathRem is Nx1 determining if waypoint should be removed
    pathRem = zeros(N,1);
    pathRem([1 end]) = [1;1];
    for i = 2:N-1
        vec1 = path(i+1,:)-path(i,:);
        vec2 = path(i,:)-path(i-1,:);
        pathRem(i) = norm(cross(vec1,vec2)) > 10^-6*norm(vec1)*norm(vec2);
    end
    ncon = sum(pathRem);
    
    % create traj and length
    traj = zeros(ncon,3);
    d = zeros(ncon,1);
    traj(1,:) = path(1,:);
    it = 2;
    for i = 2:N
        if pathRem(i)
            traj(it,:) = path(i,:);
            d(it) = d(it-1) + norm(traj(it,:)-traj(it-1,:));
            it = it+1;
        end
    end
    
    diffd = diff(d);
    V = Vf*ones(ncon-1,1);
    V(diffd<=dmin) = Vs;
    diffTs = diffd./V;
    Ts = [0; cumsum(diffTs)];
    
    % non-dimensionalize min acceleration formula
    Tf = 1;
    L0 = 0;
    Lf = 1;
    A = [1 0 0 0;
        0 0 1 0;
        -3/Tf^2 3/Tf^2 -2/Tf -1/Tf;
        2/Tf^3 -2/Tf^3 1/Tf^2 1/Tf^2];
    c = A*[L0; Lf; 0; 0];
    L = @(t) c'*[1; t; t^2; t^3];
    Ldot = @(t) c'*[0; 1; 2*t; 3*t^2];
    Ldotdot = @(t) c'*[0; 0; 2; 6*t];
end

end
