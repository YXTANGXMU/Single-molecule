function [Mat, x, y,G_index,bins] = plotPSD(data, max_PSD, min_PSD, max_meanG, min_meanG,n_bins)

n_bins = n_bins;
n_bins2=n_bins;
nanData = -1000;
s_C_bins = (max_PSD - min_PSD)/(n_bins);
s_D_bins = (max_meanG - min_meanG)/n_bins2;
PSD_index = linspace(max_PSD, min_PSD, n_bins);
G_index = [nanData, linspace(min_meanG, max_meanG, n_bins2)];
% Conduct_index = roundn(Conduct_index, -3);
% Dist_index = roundn(Dist_index, -3);
Matrix = zeros(n_bins, n_bins);
Matrix = [PSD_index', Matrix];
Matrix = [G_index; Matrix];
% data = roundn(data, -3);
[m, n] = size(data);

for i = 1: m
    
    x_ind = find((data(i, 1) < Matrix(1, :)+s_D_bins/2) & (data(i, 1) >= Matrix(1, :)-s_D_bins/2));
    y_ind = find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
    
    Matrix(y_ind, x_ind) = Matrix(y_ind, x_ind) + 1;
    
    
end

Mat = Matrix(2:end, 2:end);
x = Matrix(1, 2:end);
y = Matrix(2:end, 1);

end