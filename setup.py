from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os

sdens = "1e8"
dens = float(sdens)
lum = 1e5
ri = 500
ro = 15000
radius = 5000
nradial=100
nlam=1000
vin=1 # scaling factor for infall
# step function for abundance of molecule
Tlow = 10
Thigh = 2000
abunch3cn=1e-7

prho=4.0 #power law index
test=model(dens,prho,ri=ri,ro=ro,nphot=100000,radius=radius,nradial=nradial,lum=lum)
#test.add_gaussian_variations(0.2)

#write input filepwd
mrw=True

test.write_input(mrw=mrw)
#run model
try:
    out = test.calculate_model(ncores=20)
    print (f'stdout: {out.stdout}')
    print (f'stderr: {out.stderr}')
    test.add_lines(vin=vin, Tlow=Tlow, Thigh=Thigh, abunch3cn=abunch3cn)
    test.make_vtk()
    test.make_cube(nlam=nlam, ncores=20)
    dir = f'Dens={dens:5.2e}_Lum={lum:5.2e}_ri={ri}_ro={ro}_radius={radius}_prho={prho}_nradial={nradial}_nlam={nlam}_vin={vin}_MRW={mrw}_Tlow={Tlow}_Thigh={Thigh}_NCH3CN={abunch3cn}'
    fitsfile = f'{dir}.fits'
    test.make_fits(fitsfile)
    with open(f'{dir}.log', 'w') as logfile:
        logfile.write(f'stdout:\n {out.stdout}')
        logfile.write(f'stderr:\n {out.stderr}')
 
    cmd = f'mkdir {dir}; cp *.inp *.dat *.out *vtk {dir}; mv *.fits {dir}'
    os.system(cmd)
except Exception as e:
    print (e)

# write code to read a file

#sed
#test.sed()

#wls=[450,850,1000,2000,3000] #micron
#these functions take some time and could benefit from multithreading too.
#test.make_synth_maps(wls)
#test.make_tau_surface(wls)
