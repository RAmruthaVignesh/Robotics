function arames8_HW3_Robotics()
clc
%%
% 1)Compute the matrix exponential of the twist
twist_co = [-0.2; 0.5; 0.1; -0.4; 0.3; -0.1];
w_hat = hat(twist_co(4), twist_co(5), twist_co(6));
v = [twist_co(1); twist_co(2); twist_co(3)];
twist = [w_hat, v ; 0,0,0,];
% the matrix exponential is 
expm(twist)

%%
% 2)(a) Compute the exponential coordinates of g.
% (b) Find the screw motion to which the above rigid body motion corresponds.

la = sym('la');
g= [-0.573973  0.312409 -0.756938  1.57652; 
     0.520871 -0.573973 -0.631861 -1.28446;
    -0.631861 -0.756938  0.16672   0.875234;
       0         0         0         1];

R = g(1:3,1:3);
p = g(1:3,4);
angle_R = acos((trace(R)-1)/2)
axis_R = (1/(2*sin(angle_R)))*[(R(3,2)-R(2,3)) ;(R(1,3)-R(3,1)) ;(R(2,1)-R(1,2))];
axis_R_hat = hat(axis_R(1), axis_R(2), axis_R(3));
v = inv(((eye(3) - R)*axis_R_hat) + (axis_R*transpose(axis_R)*angle_R)) *p;
zheta = [v;axis_R]
zheta_hat = [axis_R_hat , v; 0,0,0,1];
%The exponential coordinates of g are
exp_co = zheta*angle_R
%The screw motions are
pitch = dot(axis_R,v)
q = cross(axis_R,v)
magnitude = angle_R;
%l = q + lambda* axis_R

%%
% 3) Compute the rigid-body transformation corresponding to the screw motion around the line l={[−2 3 1]T +λ[−12 20 9] T ,λ∈R} with the pitchh=1.4 for θ=140◦.
h=1.4;
theta = 140;
theta_rad = degtorad(theta);
w = [-12;20;9];
q=[-2; 3; 1]; 
w_unit = normc(w);
w_hat = hat(w_unit(1),w_unit(2),w_unit(3));
R = expm(w_hat*theta_rad)
%The Rigid body transformation is
RBT = [R ,((eye(3)-R)*q)+(h*theta_rad*w_unit);
     0,0,0, 1]


%% 
% 5) Give the geometric interpretation of the rigid body motion represented by the matrix
g=[0.327697, 0.229144, 0.916574, 0.632756;
  -0.229144, 0.960453, -0.158189, -1.26669;
 -0.916574, -0.158189, 0.367244, 1.23325;
 0,0,0,1];
R = g(1:3,1:3);
p = g(1:3,4);
angle_R = acos((trace(R)-1)/2);
axis_R = (1/(2*sin(angle_R)))*[(R(3,2)-R(2,3)) ;(R(1,3)-R(3,1)) ;(R(2,1)-R(1,2))];
axis_R_hat = hat(axis_R(1), axis_R(2), axis_R(3));
v = inv(((eye(3) - R)*axis_R_hat) + (axis_R*transpose(axis_R)*angle_R)) *p;
%The geometric interpretation twist and angle are
zheta = [v;axis_R]
angle_R = acos((trace(R)-1)/2)
zheta_hat = [axis_R_hat , v; 0,0,0,1];
h = dot(axis_R,v)
q = cross(axis_R,v)

%%
    function hat = hat(a1,a2,a3)
     hat = [0, -a3, a2 ; a3, 0, -a1 ; -a2, a1, 0];
    end
end
    