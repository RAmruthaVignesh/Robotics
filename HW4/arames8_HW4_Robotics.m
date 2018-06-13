clc
%%
%a) Compute the transformation matrix gab(t)
q = [2; -5; 3];
w = [20; -9; -12];
w = w/norm(w);
W = hat(w(1),w(2),w(3));
gab_0 = [0.510799,  0.147375,  0.846974, 0.606065;
    0.767512, -0.522024, -0.372044, 1.90451;
    0.38731,   0.840102, -0.379762, 1.12323;
    0, 0, 0, 1];
h = 1.3;
syms 't'
th = 3*(t^2) + 2*t;
Rab_t = expm(W*th);
pab_t = ((eye(3) - expm(W*th))*q) + h*th*w;
gab_t = [ Rab_t pab_t; 0 0 0 1]* gab_0

%b)Compute the spatial velocity
th = diff(th);
v = -cross(w,q)+(h*w);
z = [v;w];
vab_s = z*th

%c)
Rab = [gab_0(1:3,1:3)];
pab = [gab_0(1,4) ;gab_0(2,4) ;gab_0(3,4)];
Pab = [0,-pab(3),pab(2);
    pab(3),0,-pab(1);
    -pab(2),pab(1),0];
Ad_gab = [Rab,Pab*Rab;
    0,0,0,Rab(1,:);
    0,0,0,Rab(2,:);
    0,0,0,Rab(3,:)];
vab_b=Ad_gab*z*th

clear all;
%%
%% 2)Find a screw motion g(t) = eÎ¾tg(0) that passes through gab(t0) at time t0 and would result in the same rigid body velocity at t0

syms 't';
gab_t0 = [ 0.930114, -0.209659, -0.301549, 1.35887;
        -0.209659,  0.371024, -0.904646, -2.7234;
        0.301549,  0.904646,  0.301137,  0.738995;
        0, 0, 0, 1];
Vabs_t0 = [-1.2;0.2;0.4;-0.6;-1.0;0.8];
v = Vabs_t0(1:3); w = Vabs_t0(4:6);
W = hat(w(1),w(2),w(3));
Z = [W,v;0,0,0,0];
gt = expm(Z*t)*gab_t0

% Case when t = 0

t=0;
gt = expm(Z*t)*gab_t0

clear all;
%%
%3a)Find the spatial velocity as a function of time
syms t;
wab = [-1;-3;1];
Pab = [-t^2*cos(t) + t*sin(t) + 3;
    t^2*sin(t) + 2*t*cos(t) - 3*cos(t);
    -t^2 - 2*sin(t) + 3*t*cos(t) - 1];
Wab = hat(wab(1),wab(2), wab(3));
Rab = expm(Wab*t);
Vab_sp = (-Wab*Pab)+diff(Pab);
Vab_s = [Vab_sp;wab]

%b)Find the body velocity 
Vabb = [(Rab')*(diff(Pab));(Rab')*wab ]

%c)Find the transformation
t = 2;
Pab = [t*cos(t) - 2*t^2*sin(t);
    t^2*cos(t) - 3*t*sin(t) + 2*cos(t);
    t^2 - 3*sin(t) + t*cos(t)];
Rab = expm(Wab*1);
gab = [Rab,Pab;0,0,0,1]

clear all;
%%
% 4)What is the body velocity of the rigid body (in vector form) at t = 1s
vab_s = [0.5; -0.2; -0.3; -0.7; -0.4; 0.6];
gab_1 = [0.544224, 0.756286, -0.363114,  1.24676;
    -0.821625, 0.392997, -0.412899, -0.593889;
    -0.169567, 0.523052,  0.835262, -1.16116;
    0, 0, 0, 1];
Rab = [gab_1(1:3,1:3)];
pab = gab_1(13:15);
Pab = hat(pab(1),pab(2), pab(3));
x = zeros(3,3);
Adj_ab = [Rab' -Rab'*Pab; x Rab'];
vabb = Adj_ab*vab_s
%%
clear all;
%5)
%a)Write the formulas (but do not compute) for gab(t) and gbc(t)
%gab(t) = exp(zab^θ) gab(0)
%gbc(t) = exp(zbc^θ) gbc(0)

Vab_s = [-2; 7; 5; -1; -4; 8];
Vbc_b = [6; 9; -5; 4; 7; -1];
gab_0 = [-0.43432, -0.583051, -0.686599, 0.608676;
    -0.876071, 0.0961901, 0.472491, 1.37076;
    -0.209442, 0.806721, -0.552571, -1.28077;
    0, 0, 0, 1];
gbc_0 = [0.191896, 0.929215, 0.31581, 1.15862;
    0.623914, -0.363899, 0.691599, 0.385134;
    0.757567, 0.0643232, -0.649581, 1.30157;
    0, 0, 0, 1];
Rab0 = [gab_0(1:3,1:3)];
Rbc0 = [gbc_0(1:3,1:3)];

pab = [gab_0(1,4) ;gab_0(2,4) ;gab_0(3,4)];
pbc = [gbc_0(1,4) ;gbc_0(2,4) ;gbc_0(3,4)];
Pab = hat(pab(1),pab(2), pab(3));
Pbc = hat(pbc(1), pbc(2), pbc(3));

Adj_ab = [Rab0,Pab*Rab0;0,0,0,Rab0(1,:);0,0,0,Rab0(2,:);0,0,0,Rab0(3,:)];
Adj_bc = [Rbc0,Pbc*Rbc0; 0,0,0,Rbc0(1,:); 0,0,0,Rbc0(2,:); 0,0,0,Rbc0(3,:)];

Vbcs = (Adj_bc)*Vbc_b;
Vabb = inv(Adj_ab)*Vab_s;

%b) Compute the body velocity

gbc = inv(gbc_0);
Rbc = [gbc(1,1) gbc(1,2) gbc(1,3);
    gbc(2,1) gbc(2,2) gbc(2,3);
    gbc(3,1) gbc(3,2) gbc(3,3)];

pbc = [gbc(1,4) ;gbc(2,4) ;gbc(3,4)];
Pbc = [0,-pbc(3),pbc(2);pbc(3),0,-pbc(1);-pbc(2),pbc(1),0];
Adj_bc = [Rbc,Pbc*Rbc;
    0,0,0,Rbc(1,:);
    0,0,0,Rbc(2,:);
    0,0,0,Rbc(3,:)];

Vacb = (Adj_bc*Vabb)+Vbc_b

%c) Compute the spatial velocity
Vacs = Vab_s+(Adj_ab*Vbcs)