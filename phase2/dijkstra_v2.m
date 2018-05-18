function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an M-by-3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path.  The first
%   row is start and the last row is goal.  If no path is found, PATH is a
%   0-by-3 matrix.  Consecutive points in PATH should not be farther apart than
%   neighboring cells in the map (e.g.., if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of points that were visited while performing the search.

if nargin < 4
    astar = false;
end

m = map.occgrid;
sz = size(m);
if length(sz) < 3
    sz = [sz, 1];
end

si = map.xyzToInd(start);
gi = map.xyzToInd(goal);

res = map.res_xyz;
[x, y, z] = meshgrid(res(1)*(1:sz(1)),res(2)*(1:sz(2)),res(3)*(1:sz(3)));

% nvis: not visited, d: distance, prev: previous node
nvis = true(sz);
d = inf*ones(sz);
prev = zeros(sz);

d(si) = 0;
num_expanded = 0;

if astar
    h = zeros(sz);
    for i=1:length(m(:))
        h(i) = norm([x(i)-x(gi) y(i)-y(gi) z(i)-z(gi)]);
    end
    f = d;
    f(si)=h(si);
    finf=f;
else
    dinf = d;
end

while nvis(gi)
    if astar
        [~, ui] = min(finf(:));
    else
        [~, ui] = min(dinf(:));
    end
    num_expanded = num_expanded+1;
    if ~m(ui)
        neighb = neighbors(sz, ui);
        vit = sub2ind(sz,neighb(:,1), neighb(:,2), neighb(:,3));
        for i=1:size(vit,1)
            vi = vit(i,:);
            if nvis(vi)
                dc = d(ui) + norm([x(ui)-x(vi) y(ui)-y(vi) z(ui)-z(vi)]);
                if dc < d(vi)
                    d(vi) = dc;
                    if astar
                        finf(vi) = dc+h(vi);
                    else
                        dinf(vi) = dc;
                    end
                    prev(vi) = ui;
                end
            end
        end
    end
    nvis(ui) = 0;
    if astar
        finf(ui) = inf;
    else
        dinf(ui) = inf;
    end
end
pathi = gi;
while pathi(1) ~= si
    previ = prev(pathi(1));
    if previ < 1
        pathi = [];
        break;
    end
    pathi = [previ; pathi];
end

if isempty(pathi)
    path = [];
else
    path = [start; map.indToXYZ(pathi); goal];
end

%% post-processing
i = 1;
while i<=size(path,1)-1
    x1 = path(i,:);
    x2 = path(i+1,:);
    if x1(3)<x2(3)
        n = ceil(norm(x1-x2)/0.3);
        traj_z = x1+(x2-x1).*linspace(0,1,n)';
        path = [path(1:i,:); traj_z; path(i+1:end,:)];
        i = i+n+1;
    else
        i = i+1;
    end
end

% number of contractions
N = size(path,1);
% pathRem is Nx1 determining if waypoint should be removed
pathRem = ones(N,1);
pathRem([1 end]) = [1;1];
cnt = 1;
while cnt >0
    cnt = 0;
    i = 2;
    while i<N && ~pathRem(i)
        i = i+1;
    end
    while ~isempty(i) && i < N
        if pathRem(i)
            j = find(pathRem(1:i-1),1,'last');
            x1 = path(j,:);
            x2 = path(i+1,:);
            x = x1+(x2-x1).*linspace(0,1,ceil(norm(x1-x2)/0.015))';
            xm = x(2:end-1,:) - [0 0 0.1];
            xp = x(2:end-1,:) + [0 0 0.1];
            if ~any(map.collide([x; xm; xp]))
                pathRem(i) = 0;
                cnt = cnt+1;
            end
        end
        i = i+2;
        while i<N && ~pathRem(i)
            i = i+1;
        end
    end
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
path = traj;
end


function v = neighbors(sz, ui)
[i, j, k] = ind2sub(sz, ui);

v = zeros(26,3);
v(1:9,1) = i-1;
v(10:17,1) = i;
v(18:26,1) = i+1;
v([1:3 10:12 18:20], 2) = j-1;
v([4:6 13:14 21:23], 2) = j;
v([7:9 15:17 24:26], 2) = j+1;
v([1:3:13 15:3:24],3) = k-1;
v([2:3:11 16:3:25],3) = k;
v([3:3:12 14:3:26],3) = k+1;

valid = (v(:,1)>0) & (v(:,2)>0) &( v(:,3)>0) & (v(:,1)<=sz(1)) &...
    (v(:,2)<=sz(2)) & (v(:,3)<=sz(3));
v = v(valid,:);
end

