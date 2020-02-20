% Laser rate equation model
% by Lukas Chrostowski, 2016
%edited by Or Bahari

function main

global I;   % current

y = [0 0; 0 0]; % initial conditions for the laser (off)
time = [0 10e-9];   % simulation time

I_sweep = 0:0.1e-3:10e-3;
S_sweep = zeros (length(I_sweep),1);
N_sweep = zeros (length(I_sweep),1);
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
%%
for i = 1:length(I_sweep)
    I = I_sweep(i);
    [t, y] = ode23(@RateEqs, time, y(end,:), options);
    S_sweep(i) = y(end,1);
    N_sweep(i) = y(end,2);
end
figure(30); 
plot (I_sweep*1e3, S_sweep*4.27e-8*1e3*.5);
xlabel ('Current [mA]')
ylabel ('Optical output power [mW]')
hold on;
plot(I_sweep*1e3,power*1e3)
legend('Matlab*0.5','Lumrical')
title("LI curve Matalb Vs. Lumerical");
hold off;
%%
Ith = 3.3e-3;   % threshold current

I=2 * Ith;
load('6.6mA.mat');
y0 = [0 0]; % initial conditions for the laser (off)
[t, y] = ode23(@RateEqs, time, y0, options);
figure(31);
plot (t*1e9, y(:,1)*4.27e-8*1e3);
hold on;
%plot from lumerical I=6.6mA
%power1 = zeros(3970,1);                          
%power1 = pwr(1:3970,:);
plot (osc.time*1e9,pwr*1e3)
xlabel ('Time [ns]')
ylabel ('Optical Power [mWatts]')
title("photon density I=6.6mA Matlab Vs. Lumerical");
legend ('Matlab 6.6mA', 'Lumerical 6.6mA');

hold off;

figure(32);
I=5 * Ith;
load('16.5mA.mat');
[t, y] = ode23(@RateEqs, time, y0, options);
plot (t*1e9, y(:,1)*4.27e-8*1e3);
hold on;
plot (osc.time*1e9,pwr*1e3)
xlabel ('Time [ns]')
ylabel ('Optical Power [mWatts]')
title("photon density I=16.5mA Matlab Vs. Lumerical");
legend ('Matlab 16.5mA', 'Lumerical 16.5mA');
hold off;
end

function dy = RateEqs (t, y)

global I;   % current
S = y(1);   % Photon Number
N = y(2);   % Carrier Number

% laser parameters
tp = 3e-12; % photon lifetime, s
ts = 2e-9;  % carrier lifetime, s
G0 = 1e4;   % Modal gain, s-1
Ntr = 4e6;  % Transparency carrier number
eta = 0.9;  % Quantum efficiency
Rsp = 100e9;% Spontaneous emission rate

% constants
q = 1.6e-19;

dy = zeros(2,1);
G = G0 * ( N - Ntr);
dy(1) = (G - 1/tp) * S + Rsp;
dy(2) = eta * I / q - N/ts - G * S;

end     
      