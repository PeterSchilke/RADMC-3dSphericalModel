from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os
import sys
from multiprocessing import Pool

dens = 1e9
lum = 1e4
ri = 50
ro = 10000
radius = 200
nradial = 20
nlam=2000
vin=1 # scaling factor for infall
# step function for abundance of molecule
Tlow = 50
Thigh = 800
abunch3cn=1e-7
wls=[1360, 870, 450, 350]
ncores = 16

prho=2.0 #power law index
test=model(dens,prho,ri=ri,ro=ro,nphot=100000,radius=radius,nradial=nradial,lum=lum)
r_dev = 0.3
phi_dev = 0.3
test.add_gaussian_variations(r_dev, phi_dev)

#write input file
mrw=False


test.write_input(mrw=mrw)
#run model
try:
    dir = f'Dens={dens:5.2e}_Lum={lum:5.2e}_ri={ri}_ro={ro}_radius={radius}_prho={prho}_nradial={nradial}_nlam={nlam}_vin={vin}_MRW={mrw}_Tlow={Tlow}_Thigh={Thigh}_NCH3CN={abunch3cn}_rvar={r_dev}_phivar={phi_dev}'
    out = test.calculate_model(ncores=ncores)
    test.make_vtk()

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
 
    argument = []

    def make_cube(argument):
        nlam = argument[0]
        incl = argument[1]
        phi = argument[2] 
        ncores = argument[3]        
        dir = argument[4]
        test.make_cube(nlam=nlam, incl = incl, phi = phi, ncores=ncores)
        fitsfile = f'{dir}_incl={incl}_phi={phi}.fits'
        test.make_fits(fitsfile)
        cmd = f'mv image.out {dir}_incl={incl}_phi={phi}.out'
        os.system(cmd)
    
    for incl in [0, 45, 90]:    
        for phi in [0, 45, 90]:
            argument.append([nlam, incl, phi, ncores, dir])
            
    with Pool() as pool:
        result = pool.map_async(make_cube, argument)
        result.wait()
  


#   these functions take some time and could benefit from multithreading too.
#   test.make_synth_maps(wls)
#    test.make_tau_surface(wls)
    cmd = f'mkdir {dir}; cp *.inp  *.out {dir}; mv *.fits *log *vtk *dat {dir}'
    os.system(cmd)
#    with open(f'{dir}.log', 'a+') as logfile:
#        logfile.write(f'stdout:\n {out.stdout}')
#        logfile.write(f'stderr:\n {out.stderr}')
 
    
except Exception as e:
    print (e)
    sys.exit(1)

# write code to read a file

#sed
#test.sed()

