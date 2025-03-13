%% Q-learning based preview control for Pan-tilt
clear; close all; clc;
% Thuan, Huy, Tien, Nam
% by Khanh Tien Nguyen
%% Pan-tilt system control
%{
1. Model Pan-tilt control system
    phi_pan(k+1) = -a1_pan.phi_pan(k) -a2_pan.phi_pan(k-1)
                    +b1_pan.phi_h_pan(k)
    phi_tilt(k+1) = -a1_tilt.phi_tilt(k) -a2_tilt.phi_tilt(k-1)
                    +b1_tilt.phi_h_tilt(k)
%}
a1_pan = -1.91*10^0;
a2_pan = 9.14*10^-1;
b1_pan = 3.38*10^-3;
a1_tilt = -1.83*10^0;
a2_tilt = 8.39*10^-1;
b1_tilt = 6.09*10^-3;
% x(k+1) = Ax(k)+Bu(k)
% x(k) = [phi_ban(k) phi_pan(k-1) phi_tilt(k) phi_tilt(k-1)]'
% u(k) = [phi_h_pan phi_h_tilt]'
A = [-a1_pan -a2_pan        0        0
           1       0        0        0
           0       0 -a1_tilt -a2_tilt
           0       0        1        0];

B = [ b1_pan       0
           0       0
           0 b1_tilt
           0       0];

D = ones(4,1);
Q = 1*eye(4);
R = 10*eye(2);

% Size
n = size(A,2);
m = size(B,2);
q = size(D,2);
%% Variable
%% Parametter
T_s = 0.05;
t_end = 90;
t = 1:T_s:t_end;

%% Reference and disturbance signals
for i = 1:length(t)
    % reference signal
        if t(i) >= 10 && t(i) <= 40
            r{i} = 1*ones(n,1);
            %r{i} = sin(pi*t(i)/30)*ones(n,1);
        elseif t(i) >= 50 && t(i) <= 70
            r{i} = 2*ones(n,1);
        else
            r{i} = 0*ones(n,1);
        end
    % disturbance signal
        if t(i) >= 35 && t(i) <= 65
            d{i} = 0.003*sin(pi*t(i));
        else
            d{i} = 0;
        end
end

%% (1) Mr = 0, Md = 0;
Mr = 0;
Md = 0;

% Preview model control
% Preview model
A_r = zeros((Mr+1)*n);
A_d = zeros((Md+1)*q);
for o = 1:(Mr+1)*n-1
    A_r(o,o+1) = 1;
end
for o = 1:(Md+1)*q-1
    A_d(o,o+1) = 1;
end

if Mr == 0
    G_r = [zeros(n,n)];
else
    G_r = [zeros(n,n) -eye(n) zeros(n,(Mr-1)*n)];
end
G_d = [D zeros(n,(Md)*(q))];

A_n = [A G_r G_d
       zeros((Mr+1)*n,n) A_r zeros((Mr+1)*n,(Md+1)*q)
       zeros((Md+1)*q,n) zeros((Md+1)*q,(Mr+1)*n) A_d];

B_n = [B
       zeros((Mr+1)*n,m)
       zeros((Md+1)*q,m)];

F1 = [eye(n) zeros(n,(Mr+1)*n) zeros(n,(Md+1)*q)];
Q_n = blkdiag(Q,0*eye((Mr+1)*n),0*eye((Md+1)*q));

% The Optimal preview controller
[L1,P,~] = dlqr(A_n,B_n,Q_n,R);
K1 = -L1;
K_e = K1(:,1:n);
K_r = K1(:,n+1:(Mr+1)*n);
K_d = K1(:,n+(Mr+1)*n+1:end);

% Simulation (1)
x1{1} = zeros(n,1);
x1{2} = zeros(n,1);
u1{1} = zeros(m,1);
e1{1} = zeros(n,1);
e_s1{1} = zeros(n,1);

for i = 2:length(t)
    % error
    e1{i} = x1{i} - r{i};
    e_s1{i} = e_s1{i-1} + e1{i};
    %{
    r_s = zeros(m,1); d_s = zeros(m,1);
    if i+Mr <= length(t) && i+Md <= length(t)
        for o = 1:Mr
            r_s = r_s + K_r(:,((o-1)*n+1):((o-1)*n+n))*(r{i+o}-r{i+o-1});
        end
        for o = 1:Md
            d_s = d_s + K_d(:,((o-1)*q+1):((o-1)*q+q))*(d{i+o}-d{i+o-1});
        end
    end
    %}
    % input signal    
    u1{i} = K_e*e_s1{i};

        if i == length(t)
            break
        end
    % state
    x1{i+1} = A*x1{i} + B*u1{i} + D*d{i};
end

%% (2) Mr = 1, Md = 1;
Mr = 10;
Md = 10;

% Preview model control
% Preview model
A_r = zeros((Mr+1)*n);
A_d = zeros((Md+1)*q);
for o = 1:(Mr+1)*n-1
    A_r(o,o+1) = 1;
end
for o = 1:(Md+1)*q-1
    A_d(o,o+1) = 1;
end

if Mr == 0
    G_r = [zeros(n,n)];
else
    G_r = [zeros(n,n) -eye(n) zeros(n,(Mr-1)*n)];
end
G_d = [D zeros(n,(Md)*(q))];

A_n = [A G_r G_d
       zeros((Mr+1)*n,n) A_r zeros((Mr+1)*n,(Md+1)*q)
       zeros((Md+1)*q,n) zeros((Md+1)*q,(Mr+1)*n) A_d];

B_n = [B
       zeros((Mr+1)*n,m)
       zeros((Md+1)*q,m)];

F1 = [eye(n) zeros(n,(Mr+1)*n) zeros(n,(Md+1)*q)];
Q_n = blkdiag(Q,0*eye((Mr+1)*n),0*eye((Md+1)*q));

% The Optimal preview controller
[L2,P,~] = dlqr(A_n,B_n,Q_n,R);
K2 = -L2;
K_e = K2(:,1:n);
K_r = K2(:,n+1:(Mr+1)*n);
K_d = K2(:,n+(Mr+1)*n+1:end);

% Simulation (2)
x2{1} = zeros(n,1);
x2{2} = zeros(n,1);
u2{1} = zeros(m,1);
e2{1} = zeros(n,1);
e_s2{1} = zeros(n,1);

for i = 2:length(t)
    % error
    e2{i} = x2{i} - r{i};
    e_s2{i} = e_s2{i-1} + e2{i};
    
    r_s = zeros(m,1); d_s = zeros(m,1);
    if i+Mr <= length(t) && i+Md <= length(t)
        for o = 1:Mr
            r_s = r_s + K_r(:,((o-1)*n+1):((o-1)*n+n))*(r{i+o}-r{i+o-1});
        end
        for o = 1:Md
            d_s = d_s + K_d(:,((o-1)*q+1):((o-1)*q+q))*(d{i+o}-d{i+o-1});
        end
    end

    % input signal    
    u2{i} = K_e*e_s2{i} + r_s + d_s;

        if i == length(t)
            break
        end
    % state
    x2{i+1} = A*x2{i} + B*u2{i} + D*d{i};
end

%% (3) Mr = 20, Md = 20;
Mr = 20;
Md = 20;

% Preview model control
% Preview model
A_r = zeros((Mr+1)*n);
A_d = zeros((Md+1)*q);
for o = 1:(Mr+1)*n-1
    A_r(o,o+1) = 1;
end
for o = 1:(Md+1)*q-1
    A_d(o,o+1) = 1;
end

if Mr == 0
    G_r = [zeros(n,n)];
else
    G_r = [zeros(n,n) -eye(n) zeros(n,(Mr-1)*n)];
end
G_d = [D zeros(n,(Md)*(q))];

A_n = [A G_r G_d
       zeros((Mr+1)*n,n) A_r zeros((Mr+1)*n,(Md+1)*q)
       zeros((Md+1)*q,n) zeros((Md+1)*q,(Mr+1)*n) A_d];

B_n = [B
       zeros((Mr+1)*n,m)
       zeros((Md+1)*q,m)];

F1 = [eye(n) zeros(n,(Mr+1)*n) zeros(n,(Md+1)*q)];
Q_n = blkdiag(Q,0*eye((Mr+1)*n),0*eye((Md+1)*q));

% The Optimal preview controller
[L3,P,~] = dlqr(A_n,B_n,Q_n,R);
K3 = -L3;
K_e = K3(:,1:n);
K_r = K3(:,n+1:(Mr+1)*n);
K_d = K3(:,n+(Mr+1)*n+1:end);

% Simulation (3)
x3{1} = zeros(n,1);
x3{2} = zeros(n,1);
u3{1} = zeros(m,1);
e3{1} = zeros(n,1);
e_s3{1} = zeros(n,1);

for i = 2:length(t)
    % error
    e3{i} = x3{i} - r{i};
    e_s3{i} = e_s3{i-1} + e3{i};
    
    r_s = zeros(m,1); d_s = zeros(m,1);
    if i+Mr <= length(t) && i+Md <= length(t)
        for o = 1:Mr
            r_s = r_s + K_r(:,((o-1)*n+1):((o-1)*n+n))*(r{i+o}-r{i+o-1});
        end
        for o = 1:Md
            d_s = d_s + K_d(:,((o-1)*q+1):((o-1)*q+q))*(d{i+o}-d{i+o-1});
        end
    end

    % input signal    
    u3{i} = K_e*e_s3{i} + r_s + d_s;

        if i == length(t)
            break
        end
    % state
    x3{i+1} = A*x3{i} + B*u3{i} + D*d{i};
end

%% Plot
x1_plot = cell2mat(x1);
x2_plot = cell2mat(x2);
x3_plot = cell2mat(x3);
r_plot = cell2mat(r);
d_plot = cell2mat(d);

u1_plot = cell2mat(u1);
u2_plot = cell2mat(u2);
u3_plot = cell2mat(u3);

e1_plot = cell2mat(e1);
e2_plot = cell2mat(e2);
e3_plot = cell2mat(e3);

figure(1)
subplot(1, 2, 1);
plot(t,r_plot(1,:));
hold on
plot(t,x1_plot(1,:));
plot(t,x2_plot(1,:));
plot(t,x3_plot(1,:));
ylabel('$\theta_{pan}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$r$','$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);

subplot(1, 2, 2);
plot(t,r_plot(3,:));
hold on
plot(t,x1_plot(2,:));
plot(t,x2_plot(3,:));
plot(t,x3_plot(3,:));
ylabel('$\theta_{tilt}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$r$','$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);


figure(2)
subplot(2, 2, 1);
plot(t,u1_plot(1,:));
hold on
plot(t,u2_plot(1,:));
plot(t,u3_plot(1,:));
ylabel('$\hat{\theta}_{pan}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);

subplot(2, 2, 2);
plot(t,u1_plot(2,:));
hold on
plot(t,u2_plot(2,:));
plot(t,u3_plot(2,:));
ylabel('$\hat{\theta}_{tilt}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);

subplot(2, 2, 3);
plot(t,e1_plot(1,:));
hold on
plot(t,e2_plot(1,:));
plot(t,e3_plot(1,:));
ylabel('$e_{tilt}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);

subplot(2, 2, 4);
plot(t,e1_plot(2,:));
hold on
plot(t,e2_plot(2,:));
plot(t,e3_plot(2,:));
ylabel('$e_{tilt}$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$M_r=0,M_d=0$','$M_r=1,M_d=1$','$M_r=20,M_d=20$','Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 6);