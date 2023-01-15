clear
clc


% system
clc
numerators = {0 -39200 0; [50/9 0 0] 0 [-1 0 0]};
denominators = [1 0 0 0 0];
system = tf(numerators,denominators);
step(system)
pause



% 1. realization
clc
state_space = canon(system, "modal");
A = state_space.A
B = state_space.B
C = state_space.C
D = state_space.D
pause



% 2. controllability
clc
phi_c = ctrb(A,B)
rank(phi_c)    % => two uncontrollable modes
T = [   0       0       0       0       0       65536   0       0;
        0       0       0       0       2048    0       0       0;
        0       0       0       16      0       0       0       0;
        0       0.5     0       0       0       0       0       0;
        0       0       22.2222 0       0       0       0       0;
        0.3472  0       0       0       0       0       0       0;
        0       0       0       0       0       0       1       0;
        0       0       0       0       0       0       0       1];
T = T';
A_prime = T^(-1) * A * T;
A_11 = A_prime(1:6,1:6);
B_prime = T^(-1) * B;
B_1 = B_prime(1:6,1:3);
C_prime = C * T;
C_1 = C_prime(1:2,1:6);
sys = ss(A_11, B_1, C_1, D)
step(system,sys);
legend('System','Controlable Part','Location','NorthEast');
pause



% 3. observability
clc
phi_o = obsv(A,C)
rank(phi_o)
pause



% 4. stability
clc
%       4.a. lyapunov
eigenvalues = eig(A) % => oscillating stability

%       4.b. asymptotic
eigenvalues % => unstable

%       4.c. BIBO
system % => oscillating stability

%       4.d. T
eigenvalues
system         % => unstable
pause


% 5.1. state feedback
clc
desired_poles = [-1 -2 -3 -4 -5 -6];
K_prime = place(A_11,B_1,desired_poles);
K_prime = [K_prime zeros(3,2)];
K = K_prime * T^(-1)
A_cl = A - (B * K);
eig(A_cl)
sys_cl = ss(A_cl, B, C, D)
step(system,sys_cl);
legend('System','Close-Loop System','Location','NorthEast');
pause



% 5.2. integral state
clc
M = [A_11 B_1;-C_1 zeros(2,3)];
size(M)
rank(M)
A_int = [A_11 zeros(6,2); -C_1 zeros(2,2)];
B_int = [B_1 ; zeros(2,3)];
C_int = [C_1 zeros(2,2)];
rank(ctrb(A_int, B_int))
desired_poles_int = [-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8];
K_int = place(A_int, B_int, desired_poles_int);
K_int_1 = K_int(1:3, 1:6);
K_int_2 = -K_int(1:3, 7:8);
A_int_cl = [A_11-(B_1 * K_int_1)     (B_1 * K_int_2);
            -C_1                      zeros(2,2)];
B_int_cl = [zeros(7,3);[1 1 1]];
C_int_cl = [C_1 zeros(2)];
sys_int = ss(A_int_cl, B_int_cl, C_int_cl, 0)
step(system,sys_int);
legend('System','Integral State Feedback','Location','NorthEast');
pause



% 6. linear observer
clc
desired_poles_obs = [-0.7 -0.8 -0.9 -1.0 -1.1 -1.2];
L = place(A_11', C_1', desired_poles_obs);
A_hat = [A_11   -B_1*K(1:3,1:6); L'*C_1   A_11-L'*C_1-B_1*K(1:3,1:6)];
B_hat = [B_1; zeros(6,3)];
C_hat = eye(12);
sys_hat = ss(A_hat, B_hat, C_hat, 0)
X0 = [1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0] * 10^6;
t = 0:0.0001:0.5;
u = 0 * t;
U = [u; u; u];
[Y,T,X] = lsim(sys_hat, U, t, X0);

subplot(2,3,1)
plot(t, X(:,1),'b')
hold on
plot(t, X(:,7),'k--')

subplot(2,3,2)
plot(t, X(:,2),'b')
hold on
plot(t, X(:,8),'k--')

subplot(2,3,3)
plot(t, X(:,3),'b')
hold on
plot(t, X(:,9),'k--')

subplot(2,3,4)
plot(t, X(:,4),'b')
hold on
plot(t, X(:,10),'k--')

subplot(2,3,5)
plot(t, X(:,5),'b')
hold on
plot(t, X(:,11),'k--')

subplot(2,3,6)
plot(t, X(:,6),'b')
hold on
plot(t, X(:,12),'k--')