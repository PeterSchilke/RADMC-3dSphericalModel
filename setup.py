from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os

sdens = "1e10"
dens = float(sdens)
lum = 1e5
ri = 100
ro = 5000
radius = 500
nradial=100
nlam=1000
vin=1 # scaling factor for infall
# step function for abundance of molecule
Tlow = 10
Thigh = 2000
abunch3cn=1e-7
wls=[1360, 870, 450, 350]

prho=4.0 #power law index
test=model(dens,prho,ri=ri,ro=ro,nphot=100000,radius=radius,nradial=nradial,lum=lum)
r_dev = 2
phi_dev = 4
test.add_gaussian_variations(r_dev, phi_dev)

#write input filepwd
mrw=False

incl = 0
phi = 0

test.write_input(mrw=mrw)
#run model
try:
    dir = f'Dens={dens:5.2e}_Lum={lum:5.2e}_ri={ri}_ro={ro}_radius={radius}_prho={prho}_nradial={nradial}_nlam={nlam}_vin={vin}_MRW={mrw}_Tlow={Tlow}_Thigh={Thigh}_NCH3CN={abunch3cn}_rvar={r_dev}_phivar={phi_dev}'
    out = test.calculate_model(ncores=20)
    print (f'stdout: {out.stdout}')
    print (f'stderr: {out.stderr}')
    test.add_lines(vin=vin, Tlow=Tlow, Thigh=Thigh, abunch3cn=abunch3cn)
    test.make_vtk()
    test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=20)
    fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
    test.make_fits(fitsfile)
#   these functions take some time and could benefit from multithreading too.
#   test.make_synth_maps(wls)
#    test.make_tau_surface(wls)

    with open(f'{dir}.log', 'w') as logfile:
        logfile.write(f'stdout:\n {out.stdout}')
        logfile.write(f'stderr:\n {out.stderr}')
 
    cmd = f'mkdir {dir}; cp *.inp *.dat *.out *vtk {dir}; mv *.fits *log {dir}'
    os.system(cmd)
except Exception as e:
    print (e)

# write code to read a file

#sed
#test.sed()

