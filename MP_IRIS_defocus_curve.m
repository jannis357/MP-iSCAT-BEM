%% BEM simulation of a defocus curve of a nanoparticle in MP-IRIS with variable protective layer thickness
% Authors: Jan Pac, Thomas Juffmann
% Code adapted from https://github.com/derinsevenler/SP-IRIS-BEM/

%% Define illumination
lambda = 517;

% Illumination parameters
illum.enei = lambda; % free-space wavelength, nm
illum.pol = [1 0 0]; % e.g. linearly polarized along x-axis
illum.dir = [0 0 -1]; % unit vector, i.e. k/norm(k)

%% MP-IRIS substrate parameters
subst.nList = [1.33 1.4616 0.05+3.2i];   % Water, SiO2, Silver

subst.thickness = 60; % Protective layer thickness

%% Collection parameters
micr.collectNA = 0.8;

micr.x = 0;%linspace(-2e3, 2e3, 41); 
micr.y = 0;%linspace(-2e3, 2e3, 41);
micr.z = linspace(-5e3, 5e3, 41); %linspace(-8e3, 8e3, 41);


%% Define nanoparticle
np.material=1.44;    % Protein as mentioned in Piliarik and Sandoghdar "Direct optical sensing of single unlabelled proteins and super-resolution imaging of their binding sites", Nature Comm. 5 (2014) p. 4495
%np.material = 'gold';
np.type = 'sphere';
np.size = 10;
np.gapToSurface = 2; %nm between top of film and bottom of np
   
%% Attenuator ring
attenuation = sqrt(1.5*10^-4);
phase_val = pi/2;

%% Run BEM Simulation
scattered = np_bemsim(illum,subst,np);      % Calculate the scattered far fields

%% Calculate reflected reference field
% Uses jreftran (https://www.mathworks.com/matlabcentral/fileexchange/50923-jreftran-a-layered-thin-film-transmission-and-reflection-coefficient-calculator)
% to be compatible with complex refractive indices
[r,t,R,T,A] = jreftran_rt(lambda,[NaN,subst.thickness,NaN],[subst.nList(1), subst.nList(2), subst.nList(3)],0,0);

reflected.E = [r;0;0]*attenuation*exp(1i*phase_val);
reflected.dir = [0;0;1];

%% Calculate image
[I_ref, I_tot, I_scat] = mp_iris_image(scattered, reflected, micr, subst.nList(1));

%% Calculate normalized intensity 
I_norm=I_tot./I_ref;

%% Line plot of defocus curve
figure(2);
plot(micr.z/1000,squeeze(I_norm(1,1,:)));
xlabel('Defocus [µm]');
ylabel('Normalized Intensity');
title(np.size + " nm protein");

%% Calculate normalized intensity range (NI-range)
ni_range = max(I_norm(1,1,:)) - min(I_norm(1,1,:));

