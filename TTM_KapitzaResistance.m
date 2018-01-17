clearvars
close all
set(0,'Units', 'pixels', 'DefaultAxesFontSize',12,'DefaultAxesFontName'...
    ,'CMU Serif','DefaultAxesFontWeight','normal'); % normal bold
%========== Definition of the temperatures ==========%
Tlattice = 300;                 % lattice temperature [K]
Te = 300;                         % electron temperature [K]

% %========== gold constants ==========%
% g1=2.3*10^(16);             % e-ph coupling [W m^-3 K^-1)]
% k1=317 * Te / Tlattice;    % thermal conductivity [W m^-1 K^-1]
% Ce1=71 * Te;                   % electron heat capacity [Jm^-3 K^-3]
% Ci1=2.5*10^6;                % lattice heat capacity [J m^-3 K^-1]
% c_s_1 = 3450;                % sound speed [m/s]
% D1=k1/Ce1;                     % electron diffusivity
% tau1=Ce1/g1;                  % e-ph relaxation time
% %========== cobalt constants ==========%
% g2=4.3*10^17;               % e-ph coupling [W m^-3 K^-1)]
% k2=69 * Te / Tlattice;     % thermal conductivity [W m^-1 K^-1]
% Ce2=660 * Te;                % electron heat capacity [Jm^-3 K^-3]
% Ci2=3.5*10^6;               % lattice heat capacity [J m^-3 K^-1]
% c_s_2 = 6300;                % sound speed [m/s]
% D2=k2/Ce2;                    % electron diffusivity
% tau2=(Ce2/g2);               % e-ph relaxation time

%========== nickel constants ==========%
g1=10.9*10^(17);             % 2008 - Zhibin Lin (e-ph coupling [W m^-3 K^-1)]
k1=90.9 * Te / Tlattice;    % HANDBOOK E-10 thermal conductivity [W m^-1 K^-1]
Ce1=1077.4 * Te;                   % 2008 - Zhibin Lin     electron heat capacity [Jm^-3 K^-3]
Ci1=3.955*10^6;                % WIKIPEDIA lattice heat capacity [J m^-3 K^-1]
c_s_1 = 6040;                % HANDBOOK E-47 sound speed [m/s]
D1=k1/Ce1;                     % electron diffusivity
tau1=Ce1/g1;                  % e-ph relaxation time

%========== nickel constants ==========%
g2=10.9*10^(17);             % 2008 - Zhibin Lin (e-ph coupling [W m^-3 K^-1)]
k2=90.9 * Te / Tlattice;    % HANDBOOK E-10 thermal conductivity [W m^-1 K^-1]
Ce2=1077.4 * Te;                   % 2008 - Zhibin Lin     electron heat capacity [Jm^-3 K^-3]
Ci2=3.955*10^6;                % WIKIPEDIA lattice heat capacity [J m^-3 K^-1]
c_s_2 = 6040;                % HANDBOOK E-47 sound speed [m/s]
D2=k2/Ce2;                     % electron diffusivity
tau2=Ce2/g2;                  % e-ph relaxation time


% %========== SLG constants ==========%
% g2=4.3*10^17;               % !!!e-ph coupling [W m^-3 K^-1)]
% k2=0.96 * Te / Tlattice;     % HANDBOOK E-6 thermal conductivity [W m^-1 K^-1]
% Ce2=60 * Te;                % !!!electron heat capacity [Jm^-3 K^-3]
% Ci2=2.175*10^6;               % HANDBOOK OF GLASS PROPERTIES BANSAL&DOREMUS lattice heat capacity [J m^-3 K^-1]
% c_s_2 = 3980;                % HANDBOOK E-47 (longitudinal, 'heavy silicate flint')   sound speed [m/s]
% D2=k2/Ce2;                    % !!!electron diffusivity
% tau2=(Ce2/g2);               % !!!e-ph relaxation time

%========== other definitions ==========%

initial_dt=0.03e-12;       % time step [s]
t_final = 20e-12;              % maximum time
time_temp = 0: initial_dt : t_final;

% thicknessAu = 150e-9;      % thickness of gold [nm]
% thicknessCo = 30e-9;        % thickness of cobalt [nm]

thicknessAu = 39e-9;      % thickness of nickel [nm]
thicknessCo = 1e-9;        % thickness of !!also nickel!! SLG [nm]

R_Kap = 0;%e-10;               % Kapitza resistance [m^2 K W^-1]

%========== special acoustic mesh ==========%
delta_z1=c_s_1*initial_dt;  % mesh in gold
delta_z2=c_s_2*initial_dt;  % mesh in cobalt

%========== technical definitions ==========%
N1= fix(thicknessAu / delta_z1);     % calc. of number of points
N2= fix (thicknessCo / delta_z2);    % calc. of number of points
N_total=N1+N2;                            % total size of structure
d1=(N1-1)*delta_z1;

%========== time definition for the laser pulse ==========%
N_time=2^9 ;                               % number of points
delta_t=0.05e-12;                          % time step [ps]
tau_laser=100 * 1e-15;                  % FWHM definition
t0=0.15e-12;                                  % time shift
Te_omega=zeros(N_total,N_time);
Te_final=zeros(N_total);
Ti_omega=zeros(N_total,N_time);
Ti_final=zeros(N_total);
z_show = zeros(1, N_total);

%========== laser pulse creation ==========%
time_show = zeros(1,N_time);
f_time = zeros(1,N_time);

for index_time=1:N_time
    time=delta_t*(index_time-N_time/2-0.5);
    time_show(index_time)=time;
    f_time(index_time)=10^8*exp(-(time - 2 * t0)^2/tau_laser^2);
end

%========== Fourier transform laser pulse ==========%
delta_omega = (2*pi/delta_t)/(N_time-1);
omega_show = zeros(1,N_time);
f_omega = zeros(1,N_time);

for index_omega=1:N_time
    omega=-pi/delta_t+delta_omega*(index_omega-1);
    omega_show(index_omega)=omega;
    int1=0;
    for index_time=1:N_time
        time=delta_t*(index_time-N_time/2-0.5);
        int1=int1+f_time(index_time)*exp(i*omega*time);
    end
    f_omega(index_omega)=int1;
end

%========== inverse Fourier transform laser pulse ==========%
f_reconstr = zeros(1,N_time);

for index_time=1:N_time
    time=delta_t*(index_time-N_time/2-0.5);
    int1=0;
    for index_omega=1:N_time
        omega=-pi/delta_t+delta_omega*(index_omega-1);
        int1=int1+f_omega(index_omega)*exp(-i*omega*time);
    end
    f_reconstr(index_time)=int1/N_time;
end

figure(1)
subplot(121)
plot(time_show * 1e12,f_time,'-k',time_show * 1e12,abs(f_reconstr),'.r','LineWidth',2)
title('Laser pulse & reconstructed laser pulse');
xlabel('Time (ps)'); grid on; xlim([0 0.7]);

subplot(122)
plot(omega_show/(2*pi) * 1e-12,real(f_omega),omega_show/(2*pi) * 1e-12,imag(f_omega));
title('Real[Fourier transform of the laser pulse]');
xlabel('Frequency (THz)');

%========== Calculation of the Te ==========%
for index_omega=1:N_time
    omega=delta_omega*(index_omega-1)-pi/delta_t;
    % calculation of the Maxwell wave numbers p
    p1=sqrt(-1i*(omega/D1)*(1 - (Ci1/Ce1)/(0 + 1i*omega*tau1*Ci1/Ce1)));
    if (real(p1)<0) 
        p1 = -p1; 
    end;
    
    p2=sqrt(-1i*(omega/D2)*(1 - (Ci2/Ce2)/(0 + 1i*omega*tau2*Ci2/Ce2)));
    if (real(p2)<0) 
        p2 = -p2; 
    end;
    
    exp_plus=exp(p1*d1); exp_minus=exp(-p1*d1);
    
    discrim = (-k1 * p1 + k2 * p2 - R_Kap * k1 * k2 * p1 * p2) * k1 * p1 * exp_minus + ...
        ( k1 * p1 + k2 * p2 + R_Kap * k1 * k2 * p1 * p2) * k1 * p1 * exp_plus;
    
    C11 = f_omega(index_omega) * (k1 * p1 + k2 * p2 + k2 * p2 * R_Kap * k1 * p1) / discrim;
    C12 = f_omega(index_omega) * (k1 * p1 - k2 * p2 + k2 * p2 * R_Kap * k1 * p1) / discrim;
    C21 = C11 + C12 - R_Kap * k1 * p1 * (C11 - C12);
    % Te in gold layer
    for index_z=1:N1-1
        z=-d1+(index_z-1)*delta_z1;
        z_show(index_z)=z;
        Te_omega(index_z,index_omega)=C11*exp(-p1*z)+C12*exp(p1*z);
    end
    % Te in cobalt layer
    for index_z=N1+1:N_total
        z=(index_z-N1)*delta_z2; z_show(index_z)=z;
        Te_omega(index_z,index_omega)=C21*exp(-p2*z);
    end
    % intermediate point between two materials
    Te_omega(N1,index_omega)=0.5*(Te_omega(N1-1,index_omega)+Te_omega(N1+1,index_omega));
end

%========== Inverse Fourier transform for Te ==========%
N_pnts=length(time_temp);
T_e_total=zeros(N_pnts,N_total);            % matrix of Te(z, t)
% T_i_total=zeros(N_pnts,N_total);             % matrix of Ti(z, t)

index_time = 0;
for idTime=time_temp
    index_time=index_time + 1;
    Te_final=zeros(1,N_total);Ti_final=zeros(1,N_total);
    for index_z=1:N_total
        int1=0;
        for index_omega=1:N_time
            omega=delta_omega*(index_omega-1)-pi/delta_t;
            int1=int1+exp(-1i*omega*idTime)*Te_omega(index_z,index_omega);
        end
        Te_final(index_z)=int1;
    end
    T_e_total(index_time,:)=Te_final ./sqrt(N_time);
end
% we neglect the imaginary component
T_e_total=real(T_e_total);

%========== Calculation of the Ti ==========%
tau1_p=Ci1/g1*10^12;
tau2_p=Ci2/g2*10^12;

E1=exp(-time_temp/tau1_p);
E2=exp(-time_temp/tau2_p);
% lattice temperature in gold
for ii=1:N1-1
    buf(:,ii)=conv(T_e_total(:, ii)', E1) ./tau1_p;
end
% lattice temperature in cobalt
for ii=N1+1:N_total
    buf(:,ii)=conv(T_e_total(:,ii)', E2) ./tau2_p;
end
% intermediate point between cobalt and gold
buf(:,N1)=0.5*(buf(:,N1-1)+buf(:,N1+1));
T_i_total = buf(1: size(T_e_total, 1), :) .* delta_t;
T_i_total=real(T_i_total);

Te1_surf1=T_e_total(:,1);           % Te at air-metal interface
Te1_surf2=T_e_total(:,N1-1);     % Te at metal1-metal2 interface (in metal1)
Te2_surf2=T_e_total(:,N1+1);    % Te at metal1-metal2 interface (in metal2)

Ti1_surf1=T_i_total(:,1);             % Ti at air-metal1 interface
Ti1_surf2=T_i_total(:,N1-1);       % Ti at metal1-metal2 interface (in metal1)
Ti2_surf2=T_i_total(:,N1+1);      % Ti at metal1-metal2 interface (in metal2)

% Plotting the profiles of the T_e
hfig = figure(2);
set(hfig, 'position', [0 0 600 450])
plot(time_show * 1e12,f_time.*max(Te1_surf1)/max(f_time),...
    time_temp* 1e12,Te1_surf1,time_temp* 1e12,Te1_surf2,...
    time_temp* 1e12,Te2_surf2)
legend('Laser pulse','T_{e,surface 1}^{Au}',...
    'T_{e,surface 2}^{Au}','T_{e,surface 2}^{Co}','Fontsize', 14);
title('T_e at air-metal1 and metal1-metal2 interfaces');
xlabel('Delay time, ps','FontSize',16);
ylabel('Temperature, a. u.' , 'FontSize',16)
xlim([0 t_final * 1e12]);
box on;

% Plotting the profiles of the T_i
hfig = figure(3);
set(hfig, 'position', [0 0 600 450])
plot(time_show* 1e12,f_time.*max(Ti2_surf2)/max(f_time),...
    time_temp* 1e12,Ti1_surf1,time_temp* 1e12,Ti1_surf2,...
    time_temp* 1e12,Ti2_surf2)
legend('Laser pulse','T_{i,surface 1}^{Au}','T_{i,surface 2}^{Au}',...
    'T_{i,surface 2}^{Co}', 'Fontsize', 14);
title('T_i at air-metal1 and metal1-metal2 interfaces');
ylabel('Temperature, a. u.' , 'FontSize',16)
xlabel('Delay time, ps','FontSize',16);
xlim([0 t_final* 1e12]);
box on;

% Plotting the 2D distribution of the T_e
hfig = figure(4);
set(hfig, 'position', [0 0 600 450])
mesh((z_show + thicknessAu) * 1e9, time_temp * 1e12, T_e_total)
title('Evolution of the T_e(z, t)');
ylabel('Delay time, ps' , 'FontSize',16)
xlabel('Thickness, nm','FontSize',16);
xlim([0 (thicknessAu + thicknessCo) * 1e9])
view([0 -90])

% Plotting the 2D distribution of the T_i
hfig = figure(5);
set(hfig, 'position', [0 0 600 450])
mesh((z_show + thicknessAu) * 1e9, time_temp * 1e12, T_i_total)
title('Evolution of the T_i(z, t)');
ylabel('Delay time, ps' , 'FontSize',16)
xlabel('Thickness, nm','FontSize',16);
xlim([0 (thicknessAu + thicknessCo) * 1e9])
view([0 -90])

% Plotting sinusoidal pattern

lambda = 10                              %grating period in um
resolution=500                              %resolution
number_of_periods=3                         %number of periods
for i=1:(resolution*number_of_periods)
    x1_temp(i)=(sin(i/resolution*pi))^2;             %sinus pattern
end

Sinus_temp=transpose(x1_temp)*T_i_total(end,:);

hfig = figure(6);
set(hfig, 'position', [0 0 600 450])
mesh(z_show * 1e9, (1:(resolution*number_of_periods))/resolution*lambda,Sinus_temp)
title('Final T_i(z, t) of Grating');
ylabel('Position, um' , 'FontSize',16);
xlabel('Thickness, nm','FontSize',16);
view([0 -90])

% Calculating heatDepthProfile

number_of_z_steps=round((thicknessAu)/c_s_1/initial_dt);
this_one(1:number_of_z_steps,1)=1e9*transpose(z_show(1,1:number_of_z_steps));
this_one(1:number_of_z_steps,2)=transpose(T_i_total(end,1:number_of_z_steps));

hfig = figure(7);
plot(this_one(:,1),this_one(:,2))

ExpFit = fit(this_one(:,1),this_one(:,2),'exp1')

hfig = figure(8);
plot(ExpFit,this_one(:,1),this_one(:,2))
