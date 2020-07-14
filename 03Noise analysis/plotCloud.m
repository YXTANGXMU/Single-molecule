function [Mat, x, y,Dist_index] = plotCloud(data, maxConduct, minConduct, maxDist, minDist)

n_bins = 500;
n_bins_1 = 800;
nanData = -1000;
s_C_bins = (maxConduct - minConduct)/(n_bins_1);
s_D_bins = (maxDist - minDist)/n_bins;
Conduct_index = linspace(maxConduct, minConduct, n_bins_1);
Dist_index = [nanData, linspace(minDist, maxDist, n_bins)];
% Conduct_index = roundn(Conduct_index, -3);
% Dist_index = roundn(Dist_index, -3);
Matrix = zeros(n_bins_1, n_bins);
Matrix = [Conduct_index', Matrix];
Matrix = [Dist_index; Matrix];
% data = roundn(data, -3);
[m, n] = size(data);
% for i = 1: m
%     
%     x_ind = find((data(i, 1) < Matrix(1, :)+s_D_bins/2) & (data(i, 1) >= Matrix(1, :)-s_D_bins/2));
%     y_ind = find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
%     
%     Matrix(y_ind, x_ind) = Matrix(y_ind, x_ind) + 1;
%     
% end
x_ind=cell(m,2);
y_ind=cell(m,2);
parfor i = 1: m
    x_ind{i}=find((data(i, 1) < Matrix(1, :)+s_D_bins/2) & (data(i, 1) >= Matrix(1, :)-s_D_bins/2));
%    y_ind{i}=find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
end
parfor i = 1: m
    y_ind{i}=find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
end
for i = 1: m
    Matrix(y_ind{i}, x_ind{i}) = Matrix(y_ind{i}, x_ind{i}) + 1;
end

Mat = Matrix(2:end, 2:end);
x = Matrix(1, 2:end);
y = Matrix(2:end, 1);

