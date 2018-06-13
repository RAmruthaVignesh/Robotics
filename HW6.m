%% 1)
clear
clc
%Reference frmes same as your solution in HW5
syms l0 l1 theta1 theta2 theta3
w(7:9)= [0; l1/sqrt((l0*l0)+(l1*l1)); l0/sqrt((l0*l0)+(l1*l1))]
q(1:3) = [0;-l1;0];
q(4:6) =[0;-l1;0];
q(7:9)= [0;0;0];
w(1:3)= [0;0;1];w(4:6)= [-1;0;0]; 
z=[];
for i =1:3:9
    z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = z(1:6)';
z2 = transpose(z(7:12));
z3 = transpose(z(13:18));
%position
p_0 = [0;0;0]
%Reference Configuration
gst_0 = [eye(3), p_0; 0,0,0, 1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
gst_theta_1 = expm(z1_hat*theta1)*expm(z2_hat*theta2)*expm(z3_hat*theta3)*gst_0;
gst_theta_2 = expm(z2_hat*theta2)*expm(z3_hat*theta3)*gst_0;
gst_theta_3 = expm(z3_hat*theta3)*gst_0;
ad_inv_gst_1 = ad_inv(gst_theta_1);
ad_gst_2 = ad_inv(gst_theta_2);
ad_gst_3 = ad_inv(gst_theta_3);
z1_dag = ad_inv_gst_1*z1;
z2_dag = ad_gst_2*z2;
z3_dag = ad_gst_3*z3;
jst_b = [z1_dag, z2_dag, z3_dag]

%%
%2)
clear
syms l1 l2 theta1 theta2 theta3
%Reference frmes same as your solution in HW5
q(4:6) =[0;0;-l2];
q(7:9)= [0;0;0];
q(1:3) = [0,0,0];
w(1:3)= [0;0;1];w(4:6)= [0;-1;0]; w(7:9)=[0;-1;0];
z=[];
for i =1:3:9
    z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = transpose(z(1:6))
z2 = transpose(z(7:12))
z3 = transpose(z(13:18))
p_0 = [0;0;0];
gst_0 = [eye(3), p_0; 0,0,0, 1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
gst_theta_1 = expm(z1_hat*theta1)*gst_0;
gst_theta_2 = expm(z1_hat*theta1)*expm(z2_hat*theta2)*gst_0;
% gst_theta_3 = expm(z3_hat*theta3)*gst_0;
ad_gst_1 = ad(gst_theta_1);
ad_gst_2 = ad(gst_theta_2);
% ad_gst_3 = ad(gst_theta_3);
z1_dag = ad_gst_1*z1
z2_dag = ad_gst_2*z3
% z3_dag = ad_gst_3*z3;
jst_s = simplify([z1, z1_dag, z2_dag])

%%
%3)
jst_s_sq = simplify(transpose(jst_s) * jst_s)
x = det(jst_s_sq)
solx = solve(x, theta1, theta2, theta3)
% The singular configuration does not exist

%%
%4)
clear
syms l0 l1 d3 theta1 theta2 theta3
%Reference frames - S is at the base and T is at joint 3 and 2 intersection
q(4:6) =[0;0;l0];
% q(7:9)= [0;l1;l0];
q(1:3) = [0,0,0];
w(1:3)= [0;0;1];w(4:6)= [0;1;0]; %w(7:9)=[0;0;0];
v(1:3)=[0;0;1];
z=[];
for i =1:3:6
    z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = z(1:6)'
z2 = transpose(z(7:12))
z3 = [v(1:3),0,0,0]'
p_0 = [0;l1;l0]
gst_0 = [eye(3), p_0; 0,0,0, 1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
gst_theta_1 = expm(z1_hat*theta1);
gst_theta_2 = expm(z1_hat*theta1)*expm(z2_hat*theta2);
% gst_theta_3 = expm(z3_hat*theta3)*gst_0;
ad_gst_1 = ad(gst_theta_1);
ad_gst_2 = ad(gst_theta_2);
% ad_gst_3 = ad(gst_theta_3);
z1_dag = ad_gst_1*z2;
z2_dag = ad_gst_2*z3;
jst_s = simplify([z1, z1_dag, z2_dag])
%%
%5)
clear
syms l1 l2 theta1 theta2 theta3
% theta1=0; theta2=0; theta3=0;
q(4:6) =[0;0;0];
q(7:9)= [0;0;0];
q(1:3) = [0,0,0];
w(1:3)= [0;0;1];w(4:6)= [-1;0;0]; w(7:9)=[0;1;0];
z=[];
for i =1:3:9
     z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = z(1:6)'
z2 = transpose(z(7:12))
z3 = z(13:18)'
p_0 = [0;0;0]
gst_0 = [eye(3), p_0; 0,0,0, 1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
gst_theta_1 = expm(z1_hat*theta1)*gst_0;
gst_theta_2 = expm(z1_hat*theta1)*expm(z2_hat*theta2)*gst_0;
ad_gst_1 = ad(gst_theta_1);
ad_gst_2 = ad(gst_theta_2);
% ad_gst_3 = ad(gst_theta_3);
z1_dag = ad_gst_1*z2;
z2_dag = ad_gst_2*z3;
jst_s = [z1, z1_dag, z2_dag]
jst_s_sq = transpose(jst_s) * jst_s
x = det(jst_s_sq)
slr_conf = solve(x, theta2)


