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

