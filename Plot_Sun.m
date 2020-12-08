function [] = Plot_Sun
AU=1.495978707000000e+08;
C = imread('sun.jpg'); 
theta=0; %theta l'ho messo a caso, in realtà dovrebbe essere variabile
[x, y, z] = ellipsoid(0,0,0, 6963400/AU, 6963400/AU,6963400/AU,1E2);
%figure('Position', [0, 0, 1920*2, 1080*2]);
surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]), 'FaceColor', 'texturemap','EdgeColor','none');
axis equal;
hold on;
end