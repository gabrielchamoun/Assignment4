% Question 4

% Linear RLC Circuit
% Low Pass Filter with ~46dB of gain at DC and a corner frequency at 8Hz

clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')
global G C b

R1 = 1;
R2 = 2;
R3 = 5.89;  % extracted from sim in part 1
R4 = 0.1;
R0 = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;

%-----------------------------------------------------------
% MNA Formulation:
%-----------------------------------------------------------
b = [0;0;0;0;0;1;0;0];  % Stimuli

G = [   % Conductance
    1/R1,-1/R1,0,0,0,1,0,0;
    -1/R1,1/R1+1/R2,0,0,0,0,1,0;
    0,0,1/R3,0,0,0,-1,0;
    0,0,0,1/R4,-1/R4,0,0,1;
    0,0,0,-1/R4,1/R4+1/R0,0,0,0;
    1,0,0,0,0,0,0,0;
    0,1,-1,0,0,0,0,0;
    0,0,-alpha/R3,1,0,0,0,0;
    ];

C = [   % Energy Storage
    C1,-C1,0,0,0,0,0,0;
    -C1,C1,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-L1,0;
    0,0,0,0,0,0,0,0;
    ];

%-----------------------------------------------------------
% Time Settings:
%-----------------------------------------------------------
h = 1/1000;
t_vec = 0:h:1;

%-----------------------------------------------------------
% Backward Euler Part A: Step Response
%-----------------------------------------------------------
xn = [0;0;0;0;0;0;0;0];                 % initial values
voutA = zeros(size(G,1),numel(t_vec));  % Vector to hold outputs
vinA = zeros(1,numel(t_vec));           % Vector to hold inputs
vinA(1,1) = inputA(t_vec(1));
for n=1:numel(t_vec)-1
    voutA(:,n) = xn;                % Storing output
    vinA(n+1) = inputA(t_vec(n+1)); % Storing input
    bnn = b * vinA(n+1);            % Capturing next input
    xn = ((C./h) + G)\((C./h)*xn + bnn);    % BE Calculation
end
voutA(:,n+1) = xn;

%-----------------------------------------------------------
% Backward Euler Part B higher frequency: Sinusoidal Input
%-----------------------------------------------------------
f = 1/0.03; % Frequency in Hz
xn = [0;0;0;0;0;0;0;0];                 % initial values
voutB1 = zeros(size(G,1),numel(t_vec)); % Vector to hold outputs
vinB1 = zeros(1,numel(t_vec));          % Vector to hold inputs
vinB1(1,1) = inputB(f, t_vec(1));
for n=1:numel(t_vec)-1
    voutB1(:,n) = xn;
    vinB1(1,n+1) = inputB(f, t_vec(n+1));
    bnn = b * vinB1(1,n+1);
    xn = ((C./h) + G)\((C./h)*xn + bnn);
end
voutB1(:,n+1) = xn;

%-----------------------------------------------------------
% Backward Euler Part B lower frequency: Sinusoidal Input
%-----------------------------------------------------------
f = 1/0.3; % Frequency in Hz
xn = [0;0;0;0;0;0;0;0];                 % initial values
voutB2 = zeros(size(G,1),numel(t_vec)); % Vector to hold outputs
vinB2 = zeros(1,numel(t_vec));          % Vector to hold inputs
vinB2(1,1) = inputB(f, t_vec(1));
for n=1:numel(t_vec)-1
    voutB2(:,n) = xn;
    vinB2(1,n+1) = inputB(f, t_vec(n+1));
    bnn = b * vinB2(1,n+1);
    xn = ((C./h) + G)\((C./h)*xn + bnn);
end
voutB2(:,n+1) = xn;

%-----------------------------------------------------------
% Backward Euler Part C: Gaussian Pulse
%-----------------------------------------------------------
xn = [0;0;0;0;0;0;0;0];                 % initial values
voutC = zeros(size(G,1),numel(t_vec));  % Vector to hold outputs
vinC = zeros(1,numel(t_vec));           % Vector to hold inputs
vinC(1,1) = inputC(t_vec(1));
for n=1:numel(t_vec)-1
    voutC(:,n) = xn;
    vinC(1,n+1) = inputC(t_vec(n+1));
    bnn = b * vinC(1,n+1);
    xn = ((C./h) + G)\((C./h)*xn + bnn);
end
voutC(:,n+1) = xn;


%-----------------------------------------------------------
% Plotting A:
%-----------------------------------------------------------
% plotting the step Response of the circuit
figure("Name","Part A Response")
plot(t_vec, vinA,'red',"LineWidth",1.75);
hold on
plot(t_vec, voutA(5,:),'k',"LineWidth",1.75);
title('Part A Response', 'FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Magnitude (V)','FontSize',12);
legend({'Input','Output'});

% plotting the frequency contents of the reponse
figure("Name","Part A Spectrum")
semilogy(abs(fftshift(fft(voutA(5,:)))))
xlabel('Frequency (F)','FontSize',12);
ylabel('Magnitude','FontSize',12);
xlim([250 750])


%-----------------------------------------------------------
% Plotting B:
%-----------------------------------------------------------
% plotting the Sinusoidal Response of the circuit
figure("Name","Part B highF Response")
plot(t_vec, vinB1,'red',"LineWidth",1.75);
hold on
plot(t_vec, voutB1(5,:),'k',"LineWidth",1.75);
title('Part B highF Response', 'FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Magnitude (V)','FontSize',12);
legend({'Input','Output'});

% plotting the frequency contents of the reponse
figure("Name","Part B highF Spectrum")
semilogy(abs(fftshift(fft(voutB1(5,:)))))
xlabel('Frequency (F)','FontSize',12);
ylabel('Magnitude','FontSize',12);
xlim([250 750])


%-----------------------------------------------------------
% Plotting B:
%-----------------------------------------------------------
% plotting the Sinusoidal Response of the circuit w/ lower Frequency
figure("Name","Part B lowF Response")
plot(t_vec, vinB2,'red',"LineWidth",1.75);
hold on
plot(t_vec, voutB2(5,:),'k',"LineWidth",1.75);
title('Part B Response', 'FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Magnitude (V)','FontSize',12);
legend({'Input','Output'});

% plotting the frequency contents of the reponse
figure("Name","Part B lowF Spectrum")
semilogy(abs(fftshift(fft(voutB2(5,:)))))
xlabel('Frequency (F)','FontSize',12);
ylabel('Magnitude','FontSize',12);
xlim([250 750])


%-----------------------------------------------------------
% Plotting C:
%-----------------------------------------------------------
% plotting the Pulse Response of the circuit
figure("Name","Part C Response")
plot(t_vec, vinC,'red',"LineWidth",1.75);
hold on
plot(t_vec, voutC(5,:),'k',"LineWidth",1.75);
title('Part C Response', 'FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Magnitude (V)','FontSize',12);
legend({'Input','Output'});


% plotting the frequency contents of the reponse
figure("Name","Part C Spectrum")
semilogy(abs(fftshift(fft(voutC(5,:)))))
xlabel('Frequency (F)','FontSize',12);
ylabel('Magnitude','FontSize',12);
xlim([250 750])