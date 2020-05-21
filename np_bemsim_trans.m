function far_field = np_bemsim_trans(illum,subst,np)
% Calculate the scattered far fields of an nanoparticle using the MNPBEM
% toolbox in MP-iSCAT
% Code adapted from https://github.com/derinsevenler/SP-IRIS-BEM/

% 1. Set substrate
% ----------------
immersionMaterial 	= epsconst(subst.nList(1)^2);
filmMaterial 		= epsconst(subst.nList(2)^2);
substrateMaterial 	= epsconst(subst.nList(3)^2);
if ischar(np.material) & ( strcmp(np.material,'gold') | strcmp(lower(np.material),'au') )
	particleMaterial = epstable('gold.dat');
else
	particleMaterial = epsconst(np.material^2);
end
materials = {immersionMaterial, filmMaterial, substrateMaterial, particleMaterial};

layerZCoordinates = [0 ,-subst.thickness];  %JP: tried to change it but doesn't work then
substrate = layerstructure(materials, [1,2,3], layerZCoordinates, layerstructure.options);


% 2. Set particle
% ------------------------
if strcmp(np.type,'sphere')
	% set the np mesh intelligently
	npMeshSize = round(pi*np.size^2/32);
	particle = trisphere( npMeshSize, np.size);
elseif strcmp(np.type,'rod')
	% set the np mesh intelligently
	n1 = round(max([np.size(1), 20]));
	n2 = round(max([np.size(1)/3, 20]));
	n3 = round(max([np.size(2), 15]));
    
	particle = trirod(np.size(1), np.size(2), [n1,n2,n3]);	particle = rot(particle, 180/pi*acos(np.orientation(3)/norm(np.orientation)), [[0 -1; 1 0]*np.orientation(1:2)'; 0]);
     figure(18); 
     plot(particle, 'EdgeColor', 'b' );
     grid on
     
else
	error('np.type must be either ''sphere'' or ''rod''');
end
% np is placed so lowest point is at z = + np.gapToSurface
particle = shift(particle, [ 0, 0, - min( particle.pos( :, 3 ) ) + np.gapToSurface ] );

% 3. Run simulation
% -----------------
op = bemoptions( 'sim', 'ret',  'interp', 'curv','layer', substrate, 'waitbar', 0);
particle = comparticle( materials, { particle }, [ 4, 1 ], 1, op ); %JP: [4,1] --> means that the particle is immersed in the first entry of the materials vector!

if ~exist( 'greentab', 'var' ) || ~greentab.ismember( substrate, illum.enei, particle )
	tab = tabspace( substrate, particle );
	greentab = compgreentablayer( substrate, tab );
end
op.greentab = greentab;
bem = bemsolver( particle, op );
exc = planewave( illum.pol, illum.dir, op );
sig = bem \ exc( particle, illum.enei );


% 4. Measure scattered far fields
% -------------------------------

unitSphere = trisphere(2^10,2);
dirVecs = unitSphere.verts;  %TJ: X,Y,Z components of unit sphere
dirVecs(dirVecs(:,3)<0,:) = []; % We are only interested in forward-scattering TJ: make sure to come in with pos. k vector (JP: for a transmission setup)

spec = spectrum( dirVecs, op );  %spec.pinfty.nvec saves dirVecs   TJ: !? Do we need to specify the medium? And if so, which one?

f = farfield( spec, sig );

far_field.p = f.p;
far_field.e = f.e;              %TJ: electric field
far_field.h = f.h;              %TJ: magnetic field
far_field.enei = f.enei;        %TJ: wavelength
