# MP-iSCAT-BEM
BEM simulations with the MNPBEM toolbox for my master thesis in the group of Thomas Juffmann, University of Vienna
Group homepage: https://imaging.univie.ac.at/

The code is adapted from BEM simulations of the SP-IRIS community, which is published on https://github.com/derinsevenler/SP-IRIS-BEM/
The SP-IRIS code was used for simulations in: D. Sevenler, O. Avci, and M. S. Ünlü. "Quantitative interferometric reflectance imaging for the detection and measurement of biological nanoparticles". In: Biomedical Optics Express 8.6 (2017), pp. 2976-2989. URL: https://doi.org/10.1364/BOE.8.002976
We especially thank Derin Sevenler for sharing his code with us and a helpful discussion!

The BEM simulations are performed with the MNPBEM toolbox (https://physik.uni-graz.at/~uxh/mnpbem/html/mnpbem_product_page.html).
The toolbox is published in: U. Hohenester and A. Trügler. "MNPBEM - A Matlab toolbox for the simulation of plasmonic nanoparticles". In: Computer Physics Communications 183.2 (2012), pp. 370-381. URL: https://doi.org/10.1016/j.cpc.2011.09.009
The code needs an installed MNPBEM toolbox to be executable!

The simulation for the MP-IRIS uses 'jreftran - A layered thin film transmission and reflection coefficient calculator' by Shawn Divitt (https://de.mathworks.com/matlabcentral/fileexchange/50923-jreftran-a-layered-thin-film-transmission-and-reflection-coefficient-calculator)

-----------------------

MP_iSCAT_defocus_curve.m --> Script to simulate a defocus curve of a nanoparticle in MP-iSCAT
MP_IRIS_defocus_curve.m --> Script to simulate a defocus curve of a nanoparticle in MP-IRIS

The other files are functions needed to run the two scripts.
