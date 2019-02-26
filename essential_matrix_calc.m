%% Structure from Motion - 2 images
% Original by Alex
% Edited by Hyunggi Chang (changh95)

clear all;
close all;
clc;

%% 3D Space preparation + 2D conversion
% Create some points on the face of a cube
P_M = [
    0   0   0   0   0   0   0   0   0   1   2   1   2   1   2;
    0   1   2   0   1   2   0   1   2   0   0   0   0   0   0;
    0   0   0   -1  -1  -1  -2  -2  -2  0   0   -1  -1  -2  -2;
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1;
    ];
Num_ptS = length(P_M);

% Define image parameters for 2D conversion
L =1500;        % size of image in pixels

% Intrinsic parameter preparation
f = L;
u0 = L/2;
v0 = L/2;
 
K = [ f  0  u0;
      0  f  v0;
      0  0   1];

%% Define pose of model wrt to camera 1
euler_Rx = deg2rad(120);
euler_Ry = deg2rad(0);
euler_Rz = deg2rad(60);

Rx = [ 1    0                   0;
       0    cos(euler_Rx)     -sin(euler_Rx);
       0    sin(euler_Rx)      cos(euler_Rx) ];
Ry = [ cos(euler_Ry)    0     sin(euler_Ry);
       0                1     0;
      -sin(euler_Ry)    0     cos(euler_Ry) ];
Rz = [ cos(euler_Rz)    -sin(euler_Rz)   0;
       sin(euler_Rz)    cos(euler_Rz)   0;
       0                0               1 ];
   
R_m_c1 = Rx * Ry * Rz;   % XYZ Rotation order
T_m_c1 = [0; 0; 5];   % translation of model wrt camera

M_I1 = [ R_m_c1 T_m_c1 ];    % Extrinsic camera parameter matrix, camera 1 to model
%% Define pose of camera 2 wrt to camera 1 / model wrt to camera 2
euler_Rx = deg2rad(0);
euler_Ry = deg2rad(-25);
euler_Rz = deg2rad(0);

Rx = [ 1    0                   0;
       0    cos(euler_Rx)     -sin(euler_Rx);
       0    sin(euler_Rx)      cos(euler_Rx) ];
Ry = [ cos(euler_Ry)    0     sin(euler_Ry);
       0                1     0;
      -sin(euler_Ry)    0     cos(euler_Ry) ];
Rz = [ cos(euler_Rz)    -sin(euler_Rz)   0;
       sin(euler_Rz)     cos(euler_Rz)   0;
       0                0               1 ];
R_c2_c1 = Rx * Ry * Rz;   %XYZ rotation order
 
% Define translation of camera2 with respect to camera1
T_c2_c1 = [3; 0; 1];
 
% Figure out homogenous matrix of model wrt to camera 2
H_m_c1 = [ R_m_c1 T_m_c1 ;  0 0 0 1]; 
H_c2_c1 = [ R_c2_c1 T_c2_c1 ;   0 0 0 1]; 
H_c1_c2 = inv(H_c2_c1);
H_m_c2 = H_c1_c2 * H_m_c1;
 
% Extract rotation / translation from homogenous matrix
R_m_c2 = H_m_c2(1:3,1:3); 
T_m_c2 = H_m_c2(1:3,4);
 
% Extrinsic camera parameter matrix
M_I2 = [ R_m_c2 T_m_c2 ];


%% Render 2D images
I1 = zeros(L,L);
I2 = zeros(L,L);

% Transform points location wrt to cameras
p1 = M_I1 * P_M; % Extrinsic param * Points coordinates
p2 = M_I2 * P_M;

% 3D-2D transformation
p1(1,:) = p1(1,:) ./ p1(3,:);
p1(2,:) = p1(2,:) ./ p1(3,:);
p1(3,:) = p1(3,:) ./ p1(3,:); 

p2(1,:) = p2(1,:) ./ p2(3,:);
p2(2,:) = p2(2,:) ./ p2(3,:);
p2(3,:) = p2(3,:) ./ p2(3,:);

% Convert image points to normalized to unnormalized
u1 = K * p1;
u2 = K * p2;

% Plot images - each data point having a square dot size 4 for better
% visualisation
for i=1:length(u1)
    x = round(u1(1,i));   y = round(u1(2,i));
    I1(y-2:y+2, x-2:x+2) = 255; %Square size 4
end
figure(1), imshow(I1, []), title('View 1');

for i=1:length(u2)
    x = round(u2(1,i));   y = round(u2(2,i));
    I2(y-2:y+2, x-2:x+2) = 255; % Sqaure dot size 4
end
figure(2), imshow(I2, []), title('View 2');

disp('Points in image 1:');
disp(u1);
disp('Points in image 2:');
disp(u2);

% Save images
%imwrite(I1, 'I1.tif');
%imwrite(I2, 'I2.tif');
%% Calculate true essential matrix
% Load translation of camera 2 wrt to camera 1
t = T_c2_c1;

% Calculate skew-symmetric of translation, and multiply it by rotation of
% camera 2 wrt to camera 1.
E = [ 0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0] * R_c2_c1;
disp('True essential matrix:');
disp(E);

% Save files
%save('u1.mat', 'u1');   % Save points to files
%save('u2.mat', 'u2');
%save('E.mat', 'E');     % Save to file

%% Data point matching
% From this section, we proceed to estimate the transformation between
% camera 1 and camera 2.

% Display labelled points on the images for visualization
figure(3),imshow(I1, []);
for i=1:length(u1)
    x = round(u1(1,i));     y = round(u1(2,i));
    rectangle('Position', [x-4 y-4 8 8], 'EdgeColor', 'r');
    text(x+4, y+4, sprintf('%d', i), 'Color', 'r');
end
figure(4), imshow(I2, []);
for i=1:length(u2)
    x = round(u2(1,i));     y = round(u2(2,i));
    rectangle('Position', [x-4 y-4 8 8], 'EdgeColor', 'r');
    text(x+4, y+4, sprintf('%d', i), 'Color', 'r');
end

% Get normalized image points
p1 = inv(K)*u1;
p2 = inv(K)*u2;

%% Data pre-conditioning
% We scale and translate the image points such that the centroid of the model is at
% the origin, and the average distance of the points at the origin is
% sqrt(2).

%Extract x,y from normalized image (z=1 in normalized images).
xy = p1(1:2,:);             
N = size(xy,2);

% Calculate the centroid of points
t = sum(xy,2) / N;

% Translate all points such that the average distance of all points from
% the centroid is sqrt(2)
xnc = xy - t*ones(1,N); % vector distance of all points from the centroid
dc = sqrt(sum(xnc.^2)); % convert to scalar distance 
d_avg = sum(dc) / N;   % average distance to the origin
s = sqrt(2)/d_avg;     % the scale factor, so that avg dist is sqrt(2)
T1 = [s*eye(2), -s*t ; 0 0 1];
p1s = T1 * p1;

% Same goes for image 2
xy = p2(1:2,:);             
N = size(xy,2);
t = sum(xy,2) / N; 
xnc = xy - t*ones(1,N);  
dc = sqrt(sum(xnc.^2)); 
d_avg = sum(dc) / N;  
s = sqrt(2)/d_avg;   
T2 = [s*eye(2), -s*t ; 0 0 1];
p2s = T2 * p2;


%% Compute essential matrix E from point corresnpondences
% p1s' E p2s = 0, where p1s,p2s are the scaled image coords.
% We write out the equations in the unknowns E(i,j) A x = 0
A = [p1s(1,:)'.*p2s(1,:)'   p1s(1,:)'.*p2s(2,:)'  p1s(1,:)' ...
         p1s(2,:)'.*p2s(1,:)'   p1s(2,:)'.*p2s(2,:)'  p1s(2,:)' ...
         p2s(1,:)'             p2s(2,:)'  ones(length(p1s),1)];       
 
% Ax = 0, in which the singular vector of A corresponding to the 
%smallest sinuglar value is the last column of V in the A+UDV'. 

[U,D,V] = svd(A);
x = V(:,size(V,2));                  % get last column of V

Escale = reshape(x,3,3)'; % 3x3 matrix

% Force rank=2 and equal eigenvalues
[U,D,V] = svd(Escale);
Escale = U*diag([1 1 0])*V';

% Undo scaling
E = T1' * Escale * T2;

disp('Calculated essential matrix:');
disp(E);

%save('E.mat', 'E');     % Save to file

%% Extract motion parameters from essential matrix.
% We know that E = [tx] R, where  [tx] = [ 0 -t3 t2; t3 0 -t1; -t2 t1 0]
[U,D,V] = svd(E); % SVD of E

% From E = U diag(1,1,0) V', we extract t, which is the last column of U.
% There are 4 possible mathematical outcomes.
W = [0 -1 0; 1 0 0; 0 0 1];
Hresult_c2_c1(:,:,1) = [ U*W*V'   U(:,3) ; 0 0 0 1];
Hresult_c2_c1(:,:,2) = [ U*W*V'  -U(:,3) ; 0 0 0 1];
Hresult_c2_c1(:,:,3) = [ U*W'*V'  U(:,3) ; 0 0 0 1];
Hresult_c2_c1(:,:,4) = [ U*W'*V' -U(:,3) ; 0 0 0 1];

% make sure each rotation component is a legal rotation matrix
for k=1:4
    if det(Hresult_c2_c1(1:3,1:3,k)) < 0
        Hresult_c2_c1(1:3,1:3,k) = -Hresult_c2_c1(1:3,1:3,k);
    end
end

disp('Calculated possible poses, camera 2 to camera 1:');
disp(Hresult_c2_c1);

%% 3D Reconstruction
% Find the correct motion.  

% For matching image points p1 and p2, we know that p1 x M1 P = 0,
% p2 x M2 P = 0, where M1 and M2 are the projection matrices and P is the
% unknown 3D point.

% M1 is an identity matrix, and M2 is H_c1_c2.
M1 = [ 1 0 0 0;
       0 1 0 0;
       0 0 1 0];

% Get skew symmetric matrices for point number 1
p1x = [ 0        -p1(3,1)   p1(2,1);
        p1s(3,1)   0        -p1(1,1);
       -p1s(2,1)   p1s(1,1)   0  ];

p2x = [ 0        -p2(3,1)   p2(2,1);
        p2(3,1)   0        -p2(1,1);
       -p2(2,1)   p2(1,1)   0  ];

% Only one of the four solutions will yield a correct 3D point position that 
% is in front of both cameras (i.e. z>0 for both)
for i=1:4
    Hresult_c1_c2 = inv(Hresult_c2_c1(:,:,i));
    M2 = Hresult_c1_c2(1:3,1:4);
    
    A = [ p1x * M1; p2x * M2 ];
    % The solution to AP=0 is the singular vector of A corresponding to the
    % smallest singular value; that is, the last column of V in A=UDV'
    [U,D,V] = svd(A);
    P = V(:,4);                     % get last column of V
    P1est = P/P(4);                 % normalize

    P2est = Hresult_c1_c2 * P1est;
 
    if P1est(3) > 0 && P2est(3) > 0
        % We've found a good solution.
        Hest_c2_c1 = Hresult_c2_c1(:,:,i);
        break;      % break out of for loop; can stop searching
    end
end

% Now we have the transformation between the cameras (up to a scale factor)
fprintf('Reconstructed pose of camera2 wrt camera1:\n');
disp(Hest_c2_c1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hest_c1_c2 = inv(Hest_c2_c1);
M2est = Hest_c1_c2(1:3,:);
 
% Reconstruct point positions (these are good to the same scale factor)
fprintf('Reconstructed points wrt camera1:\n');
for i=1:length(p1)
    p1x = [ 0        -p1(3,i)   p1(2,i);
        p1(3,i)   0        -p1(1,i);
       -p1(2,i)   p1(1,i)   0  ];
    p2x = [ 0        -p2(3,i)   p2(2,i);
        p2(3,i)   0        -p2(1,i);
       -p2(2,i)   p2(1,i)   0  ];
 
    A = [ p1x * M1; p2x * M2est ];
    [U,D,V] = svd(A);
    P = V(:,4);                     % get last column of V  
    P1est(:,i) = P/P(4);                 % normalize
 
    fprintf('%f %f %f\n', P1est(1,i), P1est(2,i),P1est(3,i));
end
 
 
% Show the reconstruction result in 3D
figure;
plot3(P1est(1,:),P1est(2,:),P1est(3,:),'d');
axis equal;
axis vis3d