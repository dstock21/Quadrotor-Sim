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
Vmax = 1.8;
Vmin = 0.3;
dmin = 0.2;
dmax = 3;

% use the "persistent" keyword to keep your trajectory around
% inbetween function calls

persistent traj Ts C T Td Tdd Tddd Tdddd

%
% When called without varargin (isempty(varargin) == true), compute
% and return the desired state here.
%
if isempty(varargin)
    i = find(Ts>=t, 1, 'first');
    if t<0
        desired_state.pos = traj(0,:)';
        desired_state.vel = zeros(3,1);
        desired_state.acc = zeros(3,1);
    elseif isempty(i)
        desired_state.pos = traj(end,:)';
        desired_state.vel = zeros(3,1);
        desired_state.acc = zeros(3,1);
    else
        sz = [6,1];
        desired_state.pos = [dot(reshape(C(1,i,:),sz) , T(t));
                             dot(reshape(C(1,i,:),sz), T(t));
                             dot(reshape(C(1,i,:),sz), T(t))];
        desired_state.vel = [dot(reshape(C(1,i,:),sz), Td(t));
                             dot(reshape(C(1,i,:),sz), Td(t));
                             dot(reshape(C(1,i,:),sz), Td(t))];
        desired_state.acc = [dot(reshape(C(1,i,:),sz), Tdd(t));
                             dot(reshape(C(1,i,:),sz), Tdd(t));
                             dot(reshape(C(1,i,:),sz), Tdd(t))];
    end
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
else
    path = varargin{2};
    
    T = @(t) [t^5 t^4 t^3 t^2 t 1];
    Td = @(t) [5*t^4 4*t^3 3*t^2 2*t 1 0];
    Tdd = @(t) [20*t^3 12*t^2 6*t 2 0 0];
    Tddd = @(t) [60*t^2 24*t 6 0 0 0];
    Tdddd = @(t) [120*t 24 0 0 0 0];
    
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
    % number of contractions
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
    V = zeros(ncon-1,1);
    for i = 1:(ncon-1)
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
    Ts = cumsum(diffTs);
    
    % non-dimensionalize min acceleration formula
    C = zeros(3,ncon-1, 6);
    for i=1:3
        A = zeros(6*(ncon-1));
        b = zeros(6*(ncon-1),1);
        b(1:3) = [traj(1,i); 0; 0];
        b(4:6) = [traj(end,i); 0; 0];
        A(1:3,1:6) = [T(0); Td(0); Tdd(0)];
        A(4:6, 1:6) = [T(Ts(end)); Td(Ts(end)); Tdd(Ts(end))];
        
        % central waypoints
        for j = 2:(ncon-1)
            rmin = 6*j-5;
            rmax = 6*j;
            t = Ts(j-1);
            r = rmin:rmax;
            b(rmin) = traj(j,i);
            A(rmin, r-6) = T(t);
            
            A(r(2:6), [r-6, r]) = [T(t), -T(t);
                                Td(t), -Td(t);
                                Tdd(t), -Tdd(t);
                                Tddd(t), -Tddd(t);
                                Tdddd(t), -Tdddd(t)];
        end
        
        c = pinv(A)*b;
        C(i, :, :) = reshape(c, [6, ncon-1])';
    end
end

end
