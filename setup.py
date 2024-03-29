from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os
import sys

sdens = "1e9"
dens = float(sdens)
lum = 1e4
ri = 50
ro = 2000
radius = 100
nradial=100
nlam=2000
vin=1 # scaling factor for infall
# step function for abundance of molecule
Tlow = 50
Thigh = 800
abunch3cn=1e-8
wls=[1360, 870, 450, 350]

prho=2.0 #power law index
test=model(dens,prho,ri=ri,ro=ro,nphot=100000,radius=radius,nradial=nradial,lum=lum)
r_dev = 0.0
phi_dev = 0.1
test.add_gaussian_variations(r_dev, phi_dev)

#write input file
mrw=False


test.write_input(mrw=mrw)
#run model
try:
    dir = f'Dens={dens:5.2e}_Lum={lum:5.2e}_ri={ri}_ro={ro}_radius={radius}_prho={prho}_nradial={nradial}_nlam={nlam}_vin={vin}_MRW={mrw}_Tlow={Tlow}_Thigh={Thigh}_NCH3CN={abunch3cn}_rvar={r_dev}_phivar={phi_dev}'
    out = test.calculate_model(ncores=20)
    print (f'stdout: {out.stdout}')
    print (f'stderr: {out.stderr}')
    with open(f'{dir}.log', 'w') as logfile:
        logfile.write(f'stdout:\n {out.stdout}')
        logfile.write(f'stderr:\n {out.stderr}')

except Exception as e:
    print (f'stdout: {out.stdout}')
    print (f'stderr: {out.stderr}')
    with open(f'{dir}.log', 'w') as logfile:
        logfile.write(f'stdout:\n {out.stdout}')
        logfile.write(f'stderr:\n {out.stderr}')
        logfile.write(f'Exception: {e}')
    print (e)
    sys.exit(1)


try:
    test.add_lines(vin=vin, Tlow=Tlow, Thigh=Thigh, abunch3cn=abunch3cn)
    test.make_vtk()
 
    incl = 0
    phi = 0
    test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=20)
    fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
    test.make_fits(fitsfile)

    incl = 0
    phi = 90
    test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=20)
    fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
    test.make_fits(fitsfile)

    incl = 90
    phi = 0
    test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=20)
    fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
    test.make_fits(fitsfile)

    incl = 90
    phi = 90
    test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=20)
    fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
    test.make_fits(fitsfile)


#   these functions take some time and could benefit from multithreading too.
#   test.make_synth_maps(wls)
#    test.make_tau_surface(wls)

    with open(f'{dir}.log', 'a+') as logfile:
        logfile.write(f'stdout:\n {out.stdout}')
        logfile.write(f'stderr:\n {out.stderr}')
 
    cmd = f'mkdir {dir}; cp *.inp  *.out {dir}; mv *.fits *log *vtk *dat {dir}'
    os.system(cmd)
except Exception as e:
    print (e)
    sys.exit(1)

# write code to read a file

#sed
#test.sed()

