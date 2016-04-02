%% Sandeep Manandhar
%% University of Brugundy
%% Visual perception and Structure from Motion

%% 3D scene - a cube with two pyramids
%x  y  z
X = [...
    10 10 10;
    20 10 10;
    10 20 10;
    20 20 10;
    %lower plate
    10,10,20;
    20 10 20;
    10 20 20;
    20 20 20;
    %upper plate
    15 15 25;
    15 15 5;
    15 15 15];
X = [X ones(size(X),1)];
h3D = plotShape(X);


f = 5; %focal length
cx = 0; cy = 0; %principal point
cam1_position = [16, 16, -35];  %%viewing up
cam2_position = [14, 16, -35];  %%viewing up

% K1 = [FX 0 0; 0 FY 0;0 0 1];

plot3(h3D,cam1_position(1),cam1_position(2),cam1_position(3) , '*g'); %%camera 1
plot3(h3D,cam2_position(1),cam2_position(2),cam2_position(3), '*r'); %%camera 2


line([cam1_position(1);cam1_position(1)+5], [cam1_position(2);cam1_position(2)], [cam1_position(3);cam1_position(3)], 'Color', [1, 0, 0]);
line([cam1_position(1);cam1_position(1)], [cam1_position(2);cam1_position(2)+5], [cam1_position(3);cam1_position(3)], 'Color', [0, 1, 0]);
line([cam1_position(1);cam1_position(1)], [cam1_position(2);cam1_position(2)], [cam1_position(3);cam1_position(3)+5], 'Color', [0, 0, 1]);

line([cam2_position(1);cam2_position(1)+5], [cam2_position(2);cam2_position(2)], [cam2_position(3);cam2_position(3)], 'Color', [1, 0, 0]);
line([cam2_position(1);cam2_position(1)], [cam2_position(2);cam2_position(2)+5], [cam2_position(3);cam2_position(3)], 'Color', [0, 1, 0]);
line([cam2_position(1);cam2_position(1)], [cam2_position(2);cam2_position(2)], [cam2_position(3);cam2_position(3)+5], 'Color', [0, 0, 1]);


% line([cam1_position(1) cam1_position(1) ]',[cam1_position(2) cam1_position(2)]',[cam1_position(3) cam1_position(3)+ f]');
% line([cam2_position(1) cam2_position(1) ]',[cam2_position(2) cam2_position(2)]',[cam2_position(3) cam2_position(3)+ f]');


% for i =1:length(X)
%     line([cam1_position(1), X(i, 1)]', [cam1_position(2), X(i, 2)]', [cam1_position(3), X(i, 3)]',  'LineWidth',1,...
%    'Color',[.8 .8 .8]); 
% line([cam2_position(1), X(i, 1)]', [cam2_position(2), X(i, 2)]', [cam2_position(3), X(i, 3)]',  'LineWidth',1,...
%    'Color',[.8 .8 .8]); 
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL TWO CAMERAs
%%

[K1 R1 R2] = getStereoCam(f, 0, 0, cam1_position, cam2_position);


P = eye(3,3);
P =[P zeros(3,1)];

T1 = K1*P*R1;
x1 = world2image(T1,X);

subplot(121)
plot(x1(1,:), x1(2,:), 'ob');
T2 = K1*P*R2;
x2 = world2image(T2,X);

subplot(122)
plot(x2(1,:), x2(2,:), 'ob');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%triangulate
%% Z = fB/D; X = uZ/f; Y = vZ/f;
% B = 2; %baseline
% f = 15; %focal length
% u = xcoords
% v = ycoords
% back projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xProjected = ones(3, length(x2));
B = cam1_position(1) - cam2_position(1); f = 15;
for i =1:length(x1)
    D = abs(x1(1,i) - x2(1,i));
    Z = f*B/D;
    xProjected(1,i) = x1(1,i)*Z/f;
    xProjected(2,i) = x1(2,i)*Z/f;
    xProjected(3,i) = Z;
end
xProjected = xProjected';

plotShape(xProjected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% F matrix With point correspondences
%%xT*F*x' = 0
%%
A = theAMatrix(x1', x2');
[au, as,av] = svd(A);
F_1 = av(:,end);
F_1 = reshape(F_1, [3,3])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% With intrinsics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
theT = inv(R1)*R2
R = theT(1:3,1:3)
t = theT(:,end)
tx = [0 -t(3) t(2);
    t(3) 0 -t(1);
    -t(2) t(1) 0];
E = tx*R;
F = inv(K1)'*E*inv(K1);

x2(:,2)'*F*x1(:,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% epilinesshow
%%
for tag = 1:length(x1)
    abc = F_1*x1(:,tag);  %use F_1 for svd estimation
    x = -10:10;
    y = -(abc(1)*x + abc(3))/abc(2);
    plot(x,y );
    hold on
    plot(x2(1,tag), x2(2,tag), 'Marker', 'o');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Canonic cameras
%first camera matrix
P1 = [...
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    ];
%second camera matrix
[u s v] = svd(F_1');
e = v(:,3)  %epipole e'
ex = [0 -e(3) e(2);
    e(3) 0 -e(1);
    -e(2) e(1) 0];
M = 1/norm(e) * ex*F_1;

P2 = [M e];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% constructing 3D upto projective space
p3D = [];
for i=1:length(x1)
u1 = x1(1,i); v1 = x1(2,i);
u2 = x2(1,i); v2 = x2(2,i);

Q = [...
    u1*P1(3,:) - P1(1,:);
    v1*P1(3,:) - P1(2,:);
    u2*P2(3,:) - P2(1,:);
    v2*P2(3,:) - P2(2,:);];

[U S V] = svd(Q);
px = V(:,end);
px = px./px(4);
p3D = [p3D ;px'];
plot3(px(1), px(2), px(3), 'ob');
grid on;
hold on;
end
%plotShape(p3D);
%%
%% Residual error
plot(x1(1,:), x1(2,:), 'ob');

px1 = world2image(P1,p3D);

px2 = world2image(P2, p3D);

plot(px1(1,:), px1(2,:), '*r');
hold on;
plot(x1(1,:), x1(2,:), 'ob');


d1=norm(px1 - x1); %error in first view
d2=norm(px2 - x2); %error in second view
%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Extract R and t
 E = E';
 ea =E(:,1);
 eb =E(:,2);
 ec =E(:,3);
 
[v idx]= max([norm(cross(ea,eb)), norm(cross(ea,ec)), norm(cross(eb,ec))]);
if idx == 1 a=ea; b=eb;
elseif idx == 2 a = ea; b = ec;
elseif idx == 3 a = eb; b = ec;
end

vc = cross(a,b)/norm(cross(a,b));
va = a/norm(a);
vb = cross(vc,va);

ua = E*va/norm(E*va);
ub = E*vb/norm(E*vb);
uc = cross(ua,ub);

D = [0 1 0; -1 0 0; 0 0 1];
V = [va vb vc];
U = [ua ub uc];

tu = [U(1,3) U(2,3) U(3,3)]';
Ra = U*D*V'
Rb = U*D'*V'

Pa = [Ra tu;];
Pb = [Ra -tu;];
Pc = [Rb tu;];
Pd = [Rb -tu;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cpa = Pa*p3D(1,:)'
cpb = Pb*p3D(1,:)'
cpc = Pc*p3D(1,:)'
cpd = Pd*p3D(1,:)'
cp2 = P*p3D(1,:)'
%%choosing cpa for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = eye(3,3);
P =[P zeros(3,1)];
CP1 = K1*P*[P1; 0 0 0 1];
CP2 = K1*P*[Pa; 0 0 0 1];
cpx1 = world2image(CP1,X);
cpx2 = world2image(CP2,X);
subplot(121);
plot(cpx1(1,:), cpx1(2,:), 'ob');
subplot(122);
plot(cpx2(1,:), cpx2(2,:), 'ob');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% triangulation
xProjected = ones(3, length(cpx2));
B = 2; f = 15;
for i =1:length(x1)
    D = abs(cpx1(1,i) - cpx2(1,i));
    Z = f*B/D;
    xProjected(1,i) = cpx1(1,i)*Z/f;
    xProjected(2,i) = cpx1(2,i)*Z/f;
    xProjected(3,i) = Z;
end
xProjected = xProjected';
figure(4);
plot3(xProjected(:,1), xProjected(:,2), xProjected(:,3), 'or');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Reconstruction by DLT
cp3D = [];
for i=1:length(x1)
u1 = cpx1(1,i); v1 = cpx1(2,i);
u2 = cpx2(1,i); v2 = cpx2(2,i);

Q = [...
    u1*P1(3,:) - P1(1,:);
    v1*P1(3,:) - P1(2,:);
    u2*Pa(3,:) - Pa(1,:);
    v2*Pa(3,:) - Pa(2,:);];

[U S V] = svd(Q);
zx = V(:,end);
zx = zx./zx(4);
cp3D = [cp3D ;zx'];
plot3(zx(1), zx(2), zx(3), 'ob');
hold on;
end
%%%%%%%%%%%%%%%
%%
plotShape(cp3D);
hold on;
plotShape(X);





