% Evaluates the wavenumber-frequency spectra as measured by an array of 80 sub-resonant sensors
% as described in Damani et al. 2024 (https://arc.aiaa.org/doi/abs/10.2514/6.2024-1904)
% Uses the incompressible Chase (1980) model
% Assumes constant convection velocity (can support frequency dependent form - Line 39)
% Inputs: Boundary layer conditions (Defined in Section 1)
%         Sensor sensitivity definition (Done using discretized points in space - Section 3)
% Colormap package used
%   cubehelix - https://www.mathworks.com/matlabcentral/fileexchange/43700-cubehelix-colormap-generator-beautiful-and-versatile
% Creator: Shishir Damani - 03/06/2025
% Last updated: 03/17/25

clear;
set(0,'defaulttextInterpreter','latex') % latex axis labels
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)
set(0,'DefaultLineLineWidth',1)

%% Input Parameters (Change based on requirements)
% Flow condition (Units can be changed. Also change in lines 24-26)
Temp = 83; % Air temperature [F]
Patm = 27.88; % Air pressure [inches Hg]

% Calculating air properties
AirDensity = AirDensityCalc(Patm,'inHg',Temp,'F'); % Density of air
c0 = AirSpeedOfSoundCalc(Temp,'F'); % Speed of sound in air
viscosity = AirViscosityCalc(Temp,'F'); % Dynamic viscosity of air
pref = 20e-6; % Reference pressure for air

freq = linspace(10,3000,500); % Frequency vector

% Define BL parameters (Sample Wind Tunnel Data available here: https://doi.org/10.7294/28395098.v1)
% Change according to tunnel data or simulations
Uinfty = 38.19;  % Freestream velocity [m/s]
Ue = 38.51;  % Edge velocity [m/s]
utau = 1.342;  % Friction velocity [m/s]
delta = 66.3e-3;  % Boundary layer thickness [m]
deltas = 9.6e-3;  % Displacement thickness [m]
tauw = AirDensity*utau^2; % Wall shear stress [Pa]
Uc_Ue = 0.7; % Ratio between convective and edge velocity (Needed for model)
Ut_Ue = utau/Ue; % Ratio between friction and edge velocity (Needed for model)

kc = 2*pi*freq./Uc_Ue./Ue; % Convective wavenumber
ka = 2*pi*freq./c0; % Acoustic wavenumber

%% Defining sensor area and distinct sensor separations for CSD evaluation   
% The sensor comprises of 80 pores arranged over a rectangle
% A simple representation shown here (Circle - pores and # - boundary)
% ############   ###########    ############   ###########    
% #  ⭘  ⭘  #   #  ⭘  ⭘  #    #  ⭘  ⭘  #   #  ⭘  ⭘  # --> 
% #  ⭘  ⭘  #   #  ⭘  ⭘  #    #  ⭘  ⭘  #   #  ⭘  ⭘  #
% #  ⭘  ⭘  #   #  ⭘  ⭘  #    #  ⭘  ⭘  #   #  ⭘  ⭘  # --> Flow
% #  ⭘  ⭘  #   #  ⭘  ⭘  #    #  ⭘  ⭘  #   #  ⭘  ⭘  #
% #  ⭘  ⭘  #   #  ⭘  ⭘  #    #  ⭘  ⭘  #   #  ⭘  ⭘  # --> 
% ############   ###########    ############   ###########
L1 = 5e-3; % Dimension along flow direction
L3 = 50e-3; % Dimension along spanwise direction (across flow)
z1 = -L1/2+0.4375e-3:1.375e-3:L1/2;   % Pore distribution in stream direction for single sensor [m]
z3 = -L3/2+1.25e-3:2.5e-3:L3/2;   % Pore distribution in spanwise direction for single sensor [m]
[z1A,z3A] = ndgrid(z1,z3); % Creating a mesh of pore distribution

% Defining sensitivity over each pore
Sn = ones(size(z1A))/length(z1)/length(z3); % Assumed uniform hence summation over all pores = 1
Sm = Sn; % Sensitivity of another sensor (can be changed if assuming a different sensitivity)

% Defining sensor position in space (Array comprised of 80 sensors)
y1 = (-39:40)*5.5*1e-3; % Streamwise location of sensor centroid [m] (Assumed linear arrangement)
y3 = ones(size(y1))*0e-3; % Spanwise location of sensor centroid [m] (Assumed no spanwise separation)

% Finding distinct sensor separation for CSD evaluation
cy=y1+1i*y3; % Creates an array denoting the sensor position as complex (real - stream, imag - span)
dy=cy-cy'; % Evaluates all possible combinations of separations in space
dyu=unique(dy(real(dy)>=0)); % Condition eliminates sensor separations that are just the reverse of each other.
dy1u=unique(round(real(dyu)*1e10))/1e10; % Get unique y1 spacings, 1e10 products are to set tolerance of unique
dy1f = [dy1u;-dy1u(end-1:-1:2)+max(dy1u*2)]; % Creating an array with -ve and +ve separations
npts=length(dy1f); % Total number of discrete positions at which CSD will be evalauted
deltay1=dy1u(2)-dy1u(1); % Minimum separation between sensor [m]
dk1=2*pi/deltay1/npts; % Defining wavenumber resolution based on total number of positions in space
k1stream = (0:npts-1)*dk1; % Defining wavenumber vector [rad/m] for evaluating the model

%% Estimating the spectrum using Incompressible Chase (1980) model
% For further details on this process refer section 6.2 in PhD dissertation
om = 2*pi*freq*deltas/Ue; % Defining normalized frequency for model function
RT = utau^2*delta/Ue/(viscosity/AirDensity); % ratio of inner to outer scales

fprintf('Evaluating Chase Model for the array....\n')
%%% Cross-spectral matrix
Chase_CppArray = zeros([length(dyu) length(freq)]);
Coh_ChaseAvg = zeros([length(dyu) length(freq)]);
for n = 1:length(dyu)
    disp(n)
    % Uses defined model function (see section at end of script) to evaluate the cross-spectra between
    % indvidual pores and performs a selective summation based on sensor position
    Chase_CppArray(n,:) = ChaseCpp_Sensor(om,real(dyu(n)),imag(dyu(n)),z1A(:),z3A(:),z1A(:),z3A(:),Sn(:),Sm(:),Uc_Ue,Ut_Ue,deltas,deltas/delta,Ue,tauw); 
end

% Plot and process cross-spectral matrix
% Average (flatten) cross-spectra for duplicate y1 separations (even if y3 separation is different)
Chase_CppArrayFlat = zeros([length(dy1u) length(freq)]);
Chase_CoherenceFlat = zeros([length(dy1u) length(freq)]);
for m = 1:length(dy1u)
    a = find(abs(real(dyu)-dy1u(m))<1e-9);
    Chase_CppArrayFlat(m,:) = mean(Chase_CppArray(a,:),1); % Normalized CSD function
end

% Conjugate flipping for Fourier transform (Appending negative separations)
Chase_CppArrayFull = [conj(Chase_CppArrayFlat);Chase_CppArrayFlat(end-1:-1:2,:)];
% Applying a Hanning window prior to Fourier transform and scaling the
% normalized CSD function to yield dimensional form
Chase_CppArrayWin = Chase_CppArrayFull.*fftshift(hanning(length(Chase_CppArrayFull(:,1)),'periodic'))./sqrt(3/8);

% Fourier transform in space to yield wavenumber-frequnecy spectrum
Phipp_Chase = fft(Chase_CppArrayWin)*deltay1/2/pi; % Pa^2/(rad/s)/(rad/m)

fprintf('Done.\n')

%% Plotting autospectrum as measured by spatially averaging over sensor
figure
semilogx(freq,10*log10(2*pi*abs(Chase_CppArrayFlat(1,:))./pref^2),'k','LineWidth',2);
xlabel('Frequency [Hz]');
ylabel('$10log_{10}(\phi_{pp}(\omega)/p_{ref}^2)$ [dB]')
grid on

%% Plotting Cross-spectral density function
[freq_grid,strmloc_grid] = ndgrid(freq,dy1u);
nlevels = 41;
clvls = linspace(10,75,nlevels);

figure
set(gcf,'position',[500 500 400 300])
contourf(strmloc_grid,freq_grid,10*log10(2*pi*abs(Chase_CppArrayFlat')./pref^2),clvls,'EdgeColor','none');hold on
set(gca,'YScale','log')
cmap = cubehelix(nlevels,0.3,-0.5,1,1,[0.0,1.0],[0.0,1.0]); % CSM colorscale
colormap(cmap);colorbar
xlabel('$\Delta x_1 [m]$'); ylabel('Frequency [Hz]');
title('Cross-Spectrum Estimate')
set(gca,'TickDir','out');

%% Plotting k1-f spectrum

figure
set(gcf,'position',[500 500 500 400])
contourf(k1stream,freq,10*log10(abs(Phipp_Chase')./pref.^2./2),nlevels,'linestyle','none') % /2 because single-sided
hold on
plot(ka,freq,'--w','Linewidth',2)
plot(kc,freq,'--r','Linewidth',2)
xlim([k1stream(1) k1stream(end)/2])
ylim([freq(1) freq(end)])
cmap = cubehelix(nlevels,1.5,-1.5,1.5,1,[0.0,1.0],[0.0,1.0]);colormap(cmap);colorbar
clim([0 60]); 
xlabel('$k_1 [rad/m]$'); ylabel('Frequency [Hz]');
title('$k1-f$ Spectrum Estimate')
set(gca,'TickDir','out');

%% Line plot at interested frequency
omega_dstarU_interest = 2; % Input normalized frequency (can be changed)
[valf, nn] = min(abs(om - omega_dstarU_interest));

freq_interest = omega_dstarU_interest*Ue/deltas/2/pi;
[~,freq_ind] = min(abs(freq - freq_interest));

figure
semilogx(k1stream./kc(nn),10*log10(abs(Phipp_Chase(:,nn)).*kc(nn)./abs(Chase_CppArrayFlat(1,nn))./2),'k'); grid on
xlabel('$k_1/k_c$');
ylabel('$10log_{10}(\phi_{pp}(k_1,\omega).k_c/\phi_{pp}(\omega))$')
title(strcat('$\omega\delta^*/U_{\infty}=$',num2str(round(om(nn),3)),' / f=', num2str(round(freq(freq_ind),1)),' Hz'));

%% Line plot at interested wavenumber

k1dstar_interest = 1.12; % Input normalized wavenumber (can be changed)
[valk1, mm] = min(abs(k1stream.*deltas - k1dstar_interest));

figure
parm = Uinfty/utau^4/AirDensity^2/deltas^2;
semilogx(freq*deltas*2*pi/Uinfty,10*log10(abs(parm.*Phipp_Chase(mm,:))./2),'k'); grid on
xlabel('$\omega\delta^*/U_\infty$');
ylabel('$10log_{10}(\phi_{pp}(k_1,\omega) U_{\infty}/\tau_w^2(\delta^*)^2)$')
title(strcat('$k_1 \delta^*=$',num2str(round(k1stream(mm).*deltas,3)),' / $k_1$=', num2str(round(k1stream(mm),1)),' rad/m'))

%% Functions
%%% Sensor cross-spectrum using Chase
function res = ChaseCpp_Sensor(om,y1,y3,z1n,z3n,z1m,z3m,Sn,Sm,UcUe,UtUe,deltas,dsd,Ue,tauw)
% res - spatial crosspectrum between 2 sub-resonant sensors in dimensional form  
% y1, y3 - difference in sensor origin locations y(n)-y(m) normalized on deltastar, scalars
% z1n,z3n - coordinates of pore locations on sensor n (should be 1D arrays of same size, 1xN), relative to sensor origin normalized on deltastar
% z1m,z3m - coordinates of pore locations on sensor m (should be 1D arrays of same size, 1xM), relative to sensor origin normalized on deltastar
% Sn,Sm - pore sensitivities on sensor n and sensor m (should be 1D arrays of size 1xN and 1xM)
% om - angular frequency normalized as omega*deltastar/Ue
% dsd - deltastar/delta
% UcUe - Uc/Ue
% UtUe - utau/Ue
% Ue - Edge velocity [m/s]
% tauw - Wall shear stress [Pa]

res=zeros(size(om));
for n=1:length(om)
    % Based on Eq.13 in Damani et al. 2024 (https://arc.aiaa.org/doi/abs/10.2514/6.2024-1904)
    summand = conj(Sn*Sm').*ChaseCpp_Point(om(n),y1+z1n'-z1m,y3+z3n'-z3m,UcUe,UtUe,deltas,dsd,Ue,tauw)'; % Chase_Cpp transposed for consistency
    res(n) = sum(summand(:)); % Summation over sensor pores
end

end

%%% Pointwise cross-spectrum from Chase spectrum
function res = ChaseCpp_Point(om,x1,x3,UcUe,UtauUe,deltas,dsd,Ue,tauw)
% res - spatial crosspectrum in dimensional form 
% x1, x3 - spatial separation normalized as x/deltastar
% om - angular frequency normalized as omega*deltastar/Ue
% dsd - deltastar/delta
% UcUe - Uc/Ue
% UtauUe - utau/Ue
% deltas - Displacement thickness [m]
% Ue - Edge velocity [m/s]
% tauw - Wall shear stress [Pa]
% Based on Eq.5 in Damani et al. 2024 (https://arc.aiaa.org/doi/abs/10.2514/6.2024-3048)

h = 3;
CT = 0.0047;
CM = 0.155;
b = 0.75;
aplus = 2*pi/3.*(CT*h + CM*h);
rT = CT*h./(CT*h + CM*h);
rM = 1 - rT;
alpha = sqrt(1 + (b.*om./UcUe./dsd).^(-2));
mu = h*UtauUe/UcUe;
z1 = mu.*alpha.*om./UcUe./deltas.*x1;
z3 = alpha.*om./UcUe./deltas.*x3;
zeta = sqrt(z1.^2 + z3.^2);
fM = alpha.^(-3).*(zeta + 1 + alpha.^2.*mu.^2*(1-z1.^2)./zeta + 1i.*2.*alpha.*mu.*z1);
fT = alpha.^(-3).*(zeta + 1 + alpha.^2.*(1-z3.^2)./zeta + 1i.*2.*alpha.*mu.*z1);
omega = om./deltas.*Ue;
GammaM = aplus*rM.*tauw^2.*omega.^(-1).*exp(1i.*om./deltas./UcUe.*x1).*exp(-zeta).*fM;
GammaT = aplus*rT.*tauw^2.*omega.^(-1).*exp(1i.*om./deltas./UcUe.*x1).*exp(-zeta).*fT;

PM = rM.*aplus.*tauw^2.*omega.^(-1).*alpha.^(-3).*(1 + mu^2.*alpha.^2);
PT = rT.*aplus.*tauw^2.*omega.^(-1).*alpha.^(-3).*(1 + alpha.^2);

res = GammaM + GammaT;
res(res == inf) = PT + PM; % Set the diagonal elements to autospectrum value

end

