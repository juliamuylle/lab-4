function [] = Plot_Earth
C = imread('EarthTexture.jpg'); 
theta=0; 

[x, y, z] = ellipsoid(0,0,0, 6378.135, 6378.135,6356.750,1E2);
%figure('Position', [0, 0, 1920*2, 1080*2]);
surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]), 'FaceColor', 'texturemap','EdgeColor','none');
axis equal;
hold on;
end