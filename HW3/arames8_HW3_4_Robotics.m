%%
clc
% 4)A rigid body is moving along a screw. The coordinate transformation between the body-fixed frame B and the global reference frame A is given by
syms eeta
assume(eeta,'real')
gab = [(1/9)*((-2*cos(eeta))+(6*sin(eeta))+2), (1/9)*((-5*cos(eeta))-4), (1/9)*((-4*cos(eeta))-(3*sin(eeta))+4), (1/9)*((-12*eeta)+(23*cos(eeta))+(6*sin(eeta))+13);
        (2/9)*(cos(eeta)+(3*sin(eeta))-1), (1/9)*((-4*cos(eeta))+(3*sin(eeta))+4), (1/9)*((-5*cos(eeta))-4), (1/9)*((12*eeta)+(22*cos(eeta))-(9*sin(eeta))-31);
        (1/9)*((8*cos(eeta))+1), (2/9)*(cos(eeta)+(3*sin(eeta))-1), (1/9)*((-2*cos(eeta))+(6*sin(eeta))+2), (-2/9)*((3*eeta)+cos(eeta)+(15*sin(eeta))-10);
        0,0,0,1];
gab_0 = [(1/9)*((-2*cos(0))+(6*sin(0))+2), (1/9)*((-5*cos(0))-4), (1/9)*((-4*cos(0))-(3*sin(0))+4), (1/9)*((-12*0)+(23*cos(0))+(6*sin(0))+13);
        (2/9)*(cos(0)+(3*sin(0))-1), (1/9)*((-4*cos(0))+(3*sin(0))+4), (1/9)*((-5*cos(0))-4), (1/9)*((12*0)+(22*cos(0))-(9*sin(0))-31);
        (1/9)*((8*cos(0))+1), (2/9)*(cos(0)+(3*sin(0))-1), (1/9)*((-2*cos(0))+(6*sin(0))+2), (-2/9)*((3*0)+cos(0)+(15*sin(0))-10);
        0,0,0,1];
%calculate RBT g
g = gab* inv(gab_0)
R = g(1:3,1:3);
p = g(1:3,4);
angle_R = acos((trace(R)-1)/2)
axis_R = (1/(2*sin(angle_R)))*[(R(3,2)-R(2,3)) ;(R(1,3)-R(3,1)) ;(R(2,1)-R(1,2))];
axis_R = axis_R/norm(axis_R)
a1 = axis_R(1); a2= axis_R(2); a3= axis_R(3);
axis_R_hat = [0, -a3, a2 ; a3, 0, -a1 ; -a2, a1, 0];
v = inv(((eye(3) - R)*axis_R_hat) + (axis_R*transpose(axis_R)*angle_R)) *p;
%The screw motions are
zheta = [v;axis_R]
zheta_hat = [axis_R_hat , v; 0,0,0,1];
h = dot(axis_R,v)
q = cross(axis_R,v)
