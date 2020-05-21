function transmittedfield = mp_iscat_transmission(illum, subst,Phase,attenuation)
% Calculates the transmitted reference field in MP-iSCAT

transmittedfield.E = illum.pol*exp(1i*Phase)*attenuation;
transmittedfield.dir = illum.dir;