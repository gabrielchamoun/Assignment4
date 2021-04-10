% Question 5

clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')
global G C b

R1 = 1;
R2 = 2;
R3 = 5.89;  % extracted from sim in part 1
R4 = 0.1;
R0 = 1000;
C1 = 0.25;
Cn = 0.00001;
L1 = 0.2;
alpha = 100;

%-----------------------------------------------------------
% MNA Formulation:
%-----------------------------------------------------------
b = [0;0;1;0;0;0;0;0];  % Stimuli noise source on N3

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
    0,0,Cn,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-L1,0;
    0,0,-alpha*Cn,0,0,0,0,0;
    ];

%-----------------------------------------------------------
% Time Settings:
%-----------------------------------------------------------
h = 1/1000;
t_vec = 0:h:1;

%-----------------------------------------------------------
% Backward Euler 
%-----------------------------------------------------------
xn = [0;0;0;0;0;0;0;0];
In = 0.001*randn(1,numel(t_vec));
vout = zeros(size(G,1),numel(t_vec));
vin = zeros(1,numel(t_vec));
vin(1,1) = In(1);
for n=1:numel(t_vec)-1
    vout(:,n) = xn;
    vin(n+1) = inputC(t_vec(n+1));
    bnn = b * In(n+1);
    bnn(6) = vin(n+1);
    xn = ((C./h) + G)\((C./h)*xn + bnn);
end
vout(:,n+1) = xn;


%-----------------------------------------------------------
% Plotting 
%-----------------------------------------------------------
% plotting the thermal noise
figure("Name","Thermal Noise Histogram")
histogram(In);
title('Thermal Noise Histogram', 'FontSize',14);

% plotting the input and Transient Response of the output
figure("Name","Transient Response")
plot(t_vec, vin,'r',"LineWidth",1);
hold on
plot(t_vec, vout(5,:),'k',"LineWidth",1.75);
title('Part A Response', 'FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Magnitude (V)','FontSize',12);
legend({'Input Pulse','Output w/noise'});

% plotting the frequency contents of the reponse
figure("Name","Response Spectrum")
semilogy(abs(fftshift(fft(vout(5,:)))))
xlabel('Frequency (F)','FontSize',12);
ylabel('Magnitude','FontSize',12);
xlim([250 750])