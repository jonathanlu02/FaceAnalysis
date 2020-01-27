% Jonathan Lu
% AMATH 482
% HW1 - PCA for face recognition
% 1/19/17

close all; clc

N = 39; % 39 faces
A = []; % data matrix

%% Load images in %%%

cd('/Users/Jonathan/Documents/UW/AMATH:MATH/AMATH482/HW1/CroppedYale')
% cd('/Users/Jonathan/Documents/UW/AMATH:MATH/AMATH482/HW1/ExtendedYaleB')
currentDir = pwd;

for num = 11:N %for num = 1:N
    if num ~= 14 % for some reason YaleB14 doesn't exist?
        
        % process image file
        sampleDir = strcat(currentDir, '/', strcat('yaleB', sprintf('%02d', num)));
        cd(sampleDir)
        faceDir = strcat(sampleDir, '/*.pgm'); %% read all pgm files
        fullFile = dir(faceDir); % listing of files in directory
        
        for j = 1:length(fullFile)
            image = fullFile(j).name;
            u = imread(image);
            %[w,h] = size(u);
            u = double(u); % convert to double precision to make data matrix
            R = reshape(u, [], 1); % from matrix to col vector
            
            A = [A R];
        end
        cd('..') % move up to parent directory
    end
end
%% computing the SVD %%%

cd('/Users/Jonathan/Documents/UW/AMATH:MATH/AMATH482/HW1/CroppedYale')
load('dataMat1');
[U, S, V] = svd(A, 'econ');

%% exploring stuff %%%
cd('/Users/Jonathan/Documents/UW/AMATH:MATH/AMATH482/HW1/CroppedYale')
%load('svdData')


figure(2);
hold on

%rank = 2414
S_r = S;
B = U*S_r*V'; 
subplot(2,2,4), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 2414')
%rank = 1000
S_r(1001:end,1001:end) = 0;
B = U*S_r*V'; 
subplot(2,2,3), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 1000')
%rank = 500
S_r(501:end,501:end) = 0;
B = U*S_r*V'; 
subplot(2,2,2), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 500')
%rank = 250
S_r(251:end,251:end) = 0;
B = U*S_r*V'; 
subplot(2,2,1), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 250')
hold off

figure(1)
hold on

%rank = 100
S_r(101:end,101:end) = 0;
B = U*S_r*V'; 
subplot(2,2,4), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 100')
%rank = 50
S_r(51:end,51:end) = 0;
B = U*S_r*V'; 
subplot(2,2,3), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 50')
%rank = 25
S_r(26:end,26:end) = 0;
B = U*S_r*V'; 
subplot(2,2,2), face2 = reshape(B(:,1), 192, 168); pcolor(flipud(face2)), shading interp, colormap(gray);
xlabel('rank = 25')
%rank = 1
S_r(2:end,2:end) = 0;
B = U*S_r*V'; 
subplot(2,2,1), face1 = reshape(B(:,1), 192, 168); pcolor(flipud(face1)), shading interp, colormap(gray);
xlabel('rank = 1')

hold off

% The magnitudes of the singular values, where the latter sigma has greater
% importance
figure(3)
semilogy(diag(S), 'ko')
xlabel('\sigma_x')
ylabel('Magnitude')
title('Magnitude of \sigma_x')

figure(4)
plot(100*diag(S)/sum(diag(S)), 'ko'), xlim([0 500]);
xlabel('\sigma_x')
ylabel('Percentage (%)')

