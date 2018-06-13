clc
%%
%1a
clear
gab= [0.528383,  0.758351, -0.381725, 3.32152;
    -0.150539,  0.526175,  0.836945, 3.00205;
    0.835553, -0.384763,  0.392184, 2.63034; 
    0, 0, 0, 1];
Tb = [-3;5;-1;-4;2;3];
Rab = gab(1:3,1:3);
pab = gab(1:3,4);
pab_hat = hat(pab(1), pab(2), pab(3));
Adgab = [Rab, pab_hat*Rab ; zeros(3), Rab];
Ta = Adgab* Tb %twist in frame A
%1b
wb = [4; -3; 5; -1; 3; 2 ];
wa = inv(transpose(Adgab))* wb %Wrench in frame A

%%
%2
clear
syms l0 l1 theta1 theta2 theta3
w(7:9)= [0; l1/sqrt((l0*l0)+(l1*l1)); l0/sqrt((l0*l0)+(l1*l1))];
q(1:3) = [0;0;l0];
q(4:6) =[0;0;l0];
q(7:9)= [0;0;0];
w(1:3)= [0;0;1];w(4:6)= [-1;0;0]; 
z=[];
for i =1:3:9
    z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = z(1:6)'
z2 = transpose(z(7:12))
z3 = transpose(z(13:18))
%position
p_0 = [0;l1;l0]
%Reference Configuration
gst_0 = [eye(3), p_0; 0,0,0, 1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
%orientation
gst_theta = expm(z1_hat*theta1)*expm(z2_hat*theta2)*expm(z3_hat*theta3)*gst_0

%%
%3
clear
syms l1 l2
q(4:6) =[0;0;-l2];
q(7:9)= [0;0;0];
q(1:3) = [0,0,0];
w(1:3)= [0;0;1];w(4:6)= [0;-1;0]; w(7:9)=[0;-1;0];
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
gst_theta = expm(z1_hat)*expm(z2_hat)*expm(z3_hat)*gst_0

%%
%4)

%%
%5)
clear
syms l1 l2 l3 D th1 th2 th3 th4 th5 th6
q(7:9)= [0;l1;0];
q(10:12)= [0;0;D];
q(13:15)= [0;l1+l2;D];
q(16:18)= [0;0;D];
q(1:3) = [0,0,0];
q(4:6) =[0;0;0];
w(1:3) = [0;0;1];
w(4:6) = [1;0;0]; 
w(7:9) = [1;0;0];
w(10:12) = [0;1;0];
w(13:15) = [1;0;0];
w(16:18) = [0;1;0];
z=[];
for i =1:3:18
    z = [z, cross(-w(i: i+2),q(i: i+2)), w(i: i+2)];
end
z1 = z(1:6)'
z2 = transpose(z(7:12))
z3 = z(13:18)'
z4 = z(19:24)'
z5 = transpose(z(25:30))
z6 = z(31:36)'
 p_0 = [0;l1+l2+l3;D]
gst_0 = [eye(3), p_0; 0,0,0,1]
z1_hat = [hat(z1(4),z1(5),z1(6)), z1(1:3); 0,0,0,1];
z2_hat = [hat(z2(4),z2(5),z2(6)), z2(1:3); 0,0,0,1];
z3_hat = [hat(z3(4),z3(5),z3(6)), z3(1:3); 0,0,0,1];
z4_hat = [hat(z4(4),z4(5),z4(6)), z4(1:3); 0,0,0,1];
z5_hat = [hat(z5(4),z5(5),z5(6)), z5(1:3); 0,0,0,1];
z6_hat = [hat(z6(4),z6(5),z6(6)), z6(1:3); 0,0,0,1];
gst_theta = expm(z1_hat*th1)*expm(z2_hat*th2)*expm(z3_hat*th3)*expm(z4_hat*th4)*expm(z5_hat*th5)*expm(z6_hat*th6)*gst_0

        
