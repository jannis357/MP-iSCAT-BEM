%% BEM simulation of a defocus curve of a nanoparticle in MP-iSCAT
% Authors: Jan Pac, Thomas Juffmann
% Code adapted from https://github.com/derinsevenler/SP-IRIS-BEM/

%% Define illumination
illum.enei=517;
illum.pol=[1 0 0];
illum.dir=[0 0 1];   %come in with positive kz to be aligned with most optics textbooks

%% Define nanoparticle

%np.material='gold';
np.material=1.44;    % Protein as mentioned in Piliarik and Sandoghdar "Direct optical sensing of single unlabelled proteins and super-resolution imaging of their binding sites", Nature Comm. 5 (2014) p. 4495
np.type='sphere';
np.size=10;         % can be a vector for non spherical!
np.gapToSurface=2;  % JP: keep above 1, else there is an error
%np.orientation= ???;    % needs to be specified for non-spherical!

%% Define substrate indices
subst.nList(1)=1.33; %immersion medium of sample particle! (Mostly water in our case)
subst.nList(2)=1.52; %cover glass
subst.nList(3)=1.52; %cover glass 
subst.thickness=1; %thickness of coating

%% Define microscope stage:
%N=41;
micr.x=0;%linspace(-2e3,2e3,N);
micr.y=micr.x;
N2=50;
micr.z=linspace(-5e3,5e3,N2);
micr.collectNA=0.8;

%% Attenuator ring
micr.Phase = pi/2;     % phase shift of reference light
micr.attenuation = sqrt(1.5*10^-4);         % attenuation at phase ring! Atten is of E field --> 0.01 is an intensity attenuation of 10^-4!
mediumIndex=subst.nList(3);    % JP: Specifys (and is only used for) the collection NA of the objective in mp_iscat_image_trans!
                               % With Snell's law, mediumIndex is the refractive index of the medium that comes before the medium where the objective is.
                               % (in our case: coverslip glass)


%% Run BEM Simulation
far_field = np_bemsim_trans(illum,subst,np);          % calculate scattered fields

%% Calculate transmitted reference field
transmittedfield = mp_iscat_transmission(illum, subst,micr.Phase,micr.attenuation);

%% Calculate image
[I_ref, I_tot, I_scat] = mp_iscat_image_trans(far_field, transmittedfield, micr, mediumIndex);

%% Calculate normalized intensity 
I_norm=I_tot./I_ref; 

%% Line plot of defocus curve
figure(1);
plot(micr.z/1000,squeeze(I_norm(1,1,:)));
xlabel('Defocus [µm]');
ylabel('Normalized Intensity');
title(np.size + " nm protein");

%% Calculate normalized intensity range (NI-range)
ni_range = max(I_norm(1,1,:)) - min(I_norm(1,1,:));
