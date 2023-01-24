%% Init
close all
clear
clc
warning('off', 'images:initSize:adjustingMag');
%% Load the image
% turn the image to gray scal
% visualize the histogram for both 
FNT_SZ = 28;
img = imread('capture.png');
plotw = length(img)*10;
imshow(img)
%% 1-From 𝐶1, 𝐶2 find the horizon (vanishing) line ℎ of the plane orthogonal to the cone axis 
figure();
imshow(img);
title('Pick points for C1')
hold on;
[x, y]=getpts;
scatter(x,y,100,'filled');
%% Estimation of conic parameters

A1=[x.^2 x.*y y.^2 x y ones(size(x))];
[~,~,V1] = svd(A1);
N1 = V1(:,end);
cc1 = N1(:, 1);
[a1 b1 c1 d1 e1 f1] = deal(cc1(1),cc1(2),cc1(3),cc1(4),cc1(5),cc1(6));
C1=[a1 b1 c1 d1 e1 f1];

%% 
figure();
imshow(img);
title('Pick points for C2')
hold on;
[c, z]=getpts;
scatter(c,z,100,'filled');
%% Estimation of conic parameters

A2=[c.^2 c.*z z.^2 c z ones(size(c))];
[~,~,V2] = svd(A2);
N2 = V2(:,end);
cc2 = N2(:, 1);
[a2 b2 c2 d2 e2 f2] = deal(cc2(1),cc2(2),cc2(3),cc2(4),cc2(5),cc2(6));
C2=[a2 b2 c2 d2 e2 f2];

%% pick the conics minor axis
figure()
imshow(img);
title('Intersections and h')
[u, t]=getpts;
p1=[u(1) t(1) 1];
p2=[u(2) t(2) 1];
p3=[u(3) t(3) 1];
p4=[u(4) t(4) 1];

C1_minor_axis = cross(p1,p2);
C2_minor_axis = cross(p3,p4);
%% Find the dual conics, inverting the matrices of the conics
dual_C1_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(C1)));
dual_C2_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(C2)));

% Get the 4 tangent lines and plot them with different colors
tangent_lines = intersect(dual_C1_coeffs,dual_C2_coeffs);
% Plot the 4 tangent lines
figure()
imshow(img);
title('Intersections and h')
styles=['r-','m-','b-','y-']; widths=[2 1 2 1]; 
for k=1:length(tangent_lines)
    line = [tangent_lines(k,:), 1];
    plotLine(line,[-plotw plotw], styles(k), widths(k));
end
clear line; clear styles; clear widths; clear k;
%% 1.1.1 Intersect lines and conics
% Define the lines we found
L1 = [tangent_lines(3,:), 1];
L2 = [tangent_lines(1,:), 1];

intersection_line1 = [tangent_lines(2,:), 1];
intersection_line2 = [tangent_lines(4,:), 1];
% Intersects lines and conics
C2_L1_int = real(intersect(C2, [0 0 0 L1]));
C2_L2_int = real(intersect(C2, [0 0 0 L2]));
C1_L2_int = real(intersect(C1, [0 0 0 L2]));
C1_L1_int = real(intersect(C1, [0 0 0 L1]));

% Take only one of the 2 identically tangent points (molteplicity 2)
C2_L1_int = C2_L1_int(1,:);
C2_L2_int = C2_L2_int(1,:);
C1_L2_int = C1_L2_int(1,:);
C1_L1_int = C1_L1_int(1,:);

% Get vertical lines through tangency points of each conic
C1_major_axis = cross([C2_L2_int, 1],[C2_L1_int, 1]);
C1_major_axis = C1_major_axis/C1_major_axis(3);

C2_major_axis = cross([C1_L1_int, 1],[C1_L2_int, 1]);
C2_major_axis = C2_major_axis/C2_major_axis(3);



vertex = cross(L1, L2);
vertex = vertex / vertex(3);
vertex = vertex.';

%% Get the vanishing points
vanish_pointx = cross(C1_major_axis, C2_major_axis);
vanish_pointx = vanish_pointx / vanish_pointx(3);
vanish_pointx = vanish_pointx.';

vanish_pointz = cross(C1_minor_axis, C2_minor_axis);
vanish_pointz = vanish_pointz / vanish_pointz(3);
vanish_pointz = vanish_pointz.';

% center
center= cross(intersection_line1, intersection_line2);
center = center/center(3);

% Get the 4 intersection points
cone_axis = cross( vertex, center );
C2_points = intersect(C2, [0 0 0 cone_axis]);
C1_points = intersect(C1, [0 0 0 cone_axis]);

% Plot all tangency and vanishing points and lines
hold on;
plot(C2_L1_int(1),C2_L1_int(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(C2_L2_int(1),C2_L2_int(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(C1_L2_int(1),C1_L2_int(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(C1_L1_int(1),C1_L1_int(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(center(1),center(2),'r^','MarkerSize',6,'HandleVisibility','off');
plot(vertex(1),vertex(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(vanish_pointx(1),vanish_pointx(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(C1_points(1,1),C1_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(C1_points(2,1),C1_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(C2_points(1,1),C2_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(C2_points(2,1),C2_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plotLine(C1_major_axis , [-plotw plotw ],'-',2);
plotLine(C2_major_axis , [-plotw plotw ],'-',2);

%% Get the image of the line at the infinity
v1 = cross( C1_minor_axis, C2_minor_axis );
v2 = cross( C1_major_axis, C2_major_axis );
h = cross( v1, v2 );
plotLine(h , [-plotw plotw ], '-', 2);
% Get the image of the circular points by intersecting the ellipsis and the
% line at the infinity
circ_point1 = intersect(C2, [0 0 0 h]);
circ_point2 = intersect(C1, [0 0 0 h]);
circ_points = (circ_point1 + circ_point2) ./ 2;
legend('L2','line2','L1','line4', 'vertical1','vertical2', 'infinity');
hold off;

%% we still need to find the vanishing line in the y directions
% this can be achieved by intersecting the cone's axis with the line at
% infinity (h)

vanish_pointy = center;
vanish_pointy = vanish_pointy / vanish_pointy(3);
vanish_pointy = vanish_pointy.';

%% 2-From 𝑙1, 𝑙2, 𝐶1, 𝐶2 find the image projection 𝑎 of the cone axis
% we have already found l1, l2 and their intersections with the conics
% what is left is to find the centers of the conics
C1_center = 0.5*C1_L1_int + 0.5*C1_L2_int;
C1_center = [C1_center(1) C1_center(2) 1];
C2_center = 0.5*C2_L1_int + 0.5*C2_L2_int;
C2_center = [C2_center(1) C2_center(2) 1];
vertex2 = cross(C1_center,C2_center);
line_alpha = vertex;
% we find that the vertices are the same, just to check for the correctness
% of the solution


%% 3 Compute K from w
% Get the image of the absolute conic (w) from the constraints on
% orthogonal directions
% The image of the degenerate dual conic can be computed by multiplying the
% images of the circular points:
% 
% $\omega_{\infty} = I*J' + J*I'$
% 
I = [circ_points(1,:).'; 1];
J = [circ_points(2,:).'; 1];
l_inf = cross(I, J);
image_dual_conic = I*J.' + J*I.';

syms f, syms u0, syms v0;
K = [f, 0, u0; 0, f, v0; 0, 0, 1];

w_dual = K*K.';
w = inv(w_dual);

eqs = [ vanish_pointz.' * w * vanish_pointy;
        vanish_pointz.' * w * vanish_pointx;
        vanish_pointy.' * w * vanish_pointx;
      ];
res = solve(eqs);

f = real(double(res.f)); f = f(f > 0); f = f(1);
u0 = real(double(res.u0)); u0 = u0(1);
v0 = real(double(res.v0)); v0 = v0(1);

% Print the matrix K (normalized)
K = [f, 0, u0;0, f, v0;0, 0, 1]

%% 4  From ℎ and 𝐾, determine the orientation of the cone axis wrt to the camera reference
% this can be done easily by getting the inclination of h and since h is
% assumed to be perpindicular to the cone axis, we can get the orientation
% of the cone axis as well

axis_orientation = computeEuclideanAngles(h, C2_major_axis);

%% 5 How would you use 𝐾, ℎ, the axis orientation and the image 𝑉 of the cone vertex 
%in order to compute the cone semi-aperture angle α?
% we can find the homography from IAC
[~, DC, H] = svd(image_dual_conic);
normalization = sqrt(DC); 
normalization(3,3)=1;
H_n = normalization * H;
H = inv(H_n);

% then we can find the rectified lines of L1 or L2, wrt the cones axis

cone_axis_rectified = H'*cone_axis';
cone_axis_rectified = cone_axis_rectified/cone_axis_rectified(3);

L1_rectified = H'*L1';
L1_rectified = L1_rectified/L1_rectified(3);

L2_rectified = H'*L2';
L2_rectified = L2_rectified/L2_rectified(3);

cone_aperture_angle1 = computeEuclideanAngles(cone_axis_rectified, L1_rectified);
cone_aperture_angle2 = computeEuclideanAngles(cone_axis_rectified, L2_rectified);

figure()
hold on 
plotLine(cone_axis_rectified,[-plotw plotw ], styles(1), widths(1));
plotLine(L1_rectified,[-plotw plotw ], styles(1), widths(1));
plotLine(L2_rectified,[-plotw plotw ], styles(1), widths(1));
