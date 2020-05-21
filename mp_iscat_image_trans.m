function [I_ref, I_tot, I_scat] = mp_iscat_image_trans(scattered, reflected, micr, mediumIndex)
% Calculate the resulting image from the scattered and reference fields in
% MP-iSCAT
% Code adapted from https://github.com/derinsevenler/SP-IRIS-BEM/

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
k = 2*pi/scattered.enei;

% 1. Collect all the rays within the NA
theta = acos(scattered.p.nvec(:,3)); % theta should be positive for all these rays
collectedRays = find(theta < asin(micr.collectNA/mediumIndex));
Esca = scattered.e(collectedRays,:);
kvecs = k.* scattered.p.nvec(collectedRays,:);

% 2. Get the scattered fields at the sample points
N = 2^7; % mesh density parameter
[kx,ky] = meshgrid(linspace(min(kvecs(:,1)),max(kvecs(:,1)),N), linspace(min(kvecs(:,2)),max(kvecs(:,2)),N));
kz = sqrt(k^2 - kx.^2 - ky.^2);
dk = kx(1,2)-kx(1,1);
Ex_inf = griddata(kvecs(:,1), kvecs(:,2), Esca(:,1), kx, ky);
Ey_inf = griddata(kvecs(:,1), kvecs(:,2), Esca(:,2), kx, ky);
Ez_inf = griddata(kvecs(:,1), kvecs(:,2), Esca(:,3), kx, ky);
% we have to mesh over a rectangular range so there are some elements out of bounds in the corners
Ex_inf(isnan(Ex_inf)) = 0;
Ey_inf(isnan(Ey_inf)) = 0;
Ez_inf(isnan(Ez_inf)) = 0;

% Novotny Eq. 3.33
sca_prefactor = 1i/(2*pi); % JP: r*exp(-i*k*r) is already included in the far field solution (the far field solution comes from the BEM solver from the toolbox!)
for zi = 1:length(micr.z)
	for xi = 1:length(micr.x)
		for yi = 1:length(micr.y)
			sca_phase = exp(1i*(kx.*micr.x(xi) + ky.*micr.y(yi) + kz.*(micr.z(zi))));
			E_x_sca(xi,yi,zi) = sca_prefactor .* sum(sum(Ex_inf.*sca_phase./kz)).*dk^2;
			E_y_sca(xi,yi,zi) = sca_prefactor .* sum(sum(Ey_inf.*sca_phase./kz)).*dk^2;
		end
	end
end

% 3. Get the reference fields at the sample points. It's a single plane
% wave, so we don't need to grid and integrate.
[xx, yy,zz] = meshgrid(micr.x,micr.y,micr.z);
krefl = k*reflected.dir;
ref_phase = exp(1i* (krefl(1).*xx + krefl(2).*yy + krefl(3)*(zz)) );
E_x_ref = ref_phase.*reflected.E(1);
E_y_ref = ref_phase.*reflected.E(2);

% 4. Get the measurable intensity
I_ref = abs(E_x_ref).^2 + abs(E_y_ref).^2;
I_scat = abs(E_x_sca).^2 + abs(E_y_sca).^2;
I_tot = abs(E_x_sca + E_x_ref).^2 + abs(E_y_sca + E_y_ref).^2;
end