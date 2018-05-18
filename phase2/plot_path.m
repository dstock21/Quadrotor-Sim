function plot_path(map, path)
% PLOT_PATH Visualize a path through an environment
%   PLOT_PATH(map, path) creates a figure showing a path through the
%   environment.  path is an N-by-3 matrix where each row corresponds to the
%   (x, y, z) coordinates of one point along the path.
m = map.occgrid;

xyz = map.indToXYZ((1:length(m(:)))');
xyzo = xyz(logical(m(:)),:);
xyzno = xyz(~logical(m(:)),:);

figure(1);
hold on;
scatter3(xyzo(:,1), xyzo(:,2), xyzo(:,3),1,'r', 'filled');
scatter3(xyzno(:,1), xyzno(:,2), xyzno(:,3),1, 'b', 'filled');
plot3(path(:,1), path(:,2), path(:,3), 'k');
xlabel('x');
ylabel('y');
zlabel('z');
end