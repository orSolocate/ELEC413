% Laser rate equation model
% by Lukas Chrostowski, 2016

function main

global I;   % current

time = [0 20e-9];   % simulation time

options = odeset('RelTol',1e-6,'AbsTol',1e-6);
Ith = 3.3e-3;   % threshold current

I=5 * Ith;
neff=3.7;
L=0.05;%5cm
c=3e8;
tau=2*L/c;
tau_in=2*neff*L/c;
R3=0.1/100;
y0 = [0 0 0]; % initial conditions for the laser (off)
[t, y] = ode23(@RateEqs, time, y0, options);
figure;
plot (t*1e9, y(:,1)*4.27e-5);
xlabel ('Time [ns]')
ylabel ('Power [mwatts]')
hold on;

%feedback
time_period = 20e-9; options = ddeset('MaxStep', 1e-4*time_period);
Z = [0 0 0]; 
sol = dde23(@(t,y,Z) ratedde(t, y, Z, tau, R3, tau_in, I),  [tau], [1e6; 0; 0.1e6], [0 time_period], options);
plot (sol.x.*1e9, sol.y(1,:)*4.27e-5);

legend ('3 x Ith no feedback', '3 x Ith with feedback');
hold off;
end

function dy = RateEqs (t, y)

global I;   % current
S = y(1);   % Photon Number
phi=y(2);
N = y(3);   % Carrier Number

% laser parameters
tp = 3e-12; % photon lifetime, s
ts = 2e-9;  % carrier lifetime, s
G0 = 1e4;   % Modal gain, s-1
Ntr = 4e6;  % Transparency carrier number
eta = 0.9;  % Quantum efficiency
Rsp = 100e9;% Spontaneous emission rate

% constants
q = 1.6e-19;
alpha=4.0;

dy = zeros(2,1);
G = G0 * ( N - Ntr);
Nth=Ntr+1/(G0*tp);
dy(1) = (G - 1/tp) * S + Rsp;
dy(2)=0.5*alpha*G0*(N-Nth);
dy(3) = eta * I / q - N/ts - G * S;

end     


function dydt = ratedde(t, y, Z, tau, R3, tau_in, I)

S = y(1);   % Photon Number
phi=y(2);   %
N = y(3);   % Carrier Number

% laser parameters
tp = 3e-12; % photon lifetime, s
ts = 2e-9;  % carrier lifetime, s
G0 = 1e4;   % Modal gain, s-1
Ntr = 4e6;  % Transparency carrier number
eta = 0.9;  % Quantum efficiency
Rsp = 100e9;% Spontaneous emission rate
% constants
q = 1.6e-19;
alpha=4.0;
omega0   = 193.414e12*2*pi; % Frequency

% Set-up the phase term and feedback:
ylag1 = Z(:,1);
Slag1 = ylag1(1);
philag1 = ylag1(2);

R2= 2.51e-5; 
kappa = (1-R2)*sqrt(R3)/(tau_in*sqrt(R2));

dydt = zeros(2,1);
G = G0 * ( N - Ntr);
Nth=Ntr+1/(G0*tp);
dydt(1) = (G - 1/tp) * S + Rsp+ 2*kappa*sqrt(Slag1*S)*cos(omega0*tau + phi - philag1);
dydt(2)=0.5*alpha*G0*(N-Nth)- kappa*sqrt(Slag1/S)*sin(omega0*tau + phi - philag1);
dydt(3) = eta * I / q - N/ts - G * S;

end     
      