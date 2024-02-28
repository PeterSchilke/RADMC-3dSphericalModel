from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os

sdens = "1e7"
dens = float(sdens)
lum = 1e5
ri = 150
ro = 15000
radius = 5000
nradial=100
nlam=30



prho=4.0 #power law index
test=model(dens,prho,ri=ri,ro=ro,nphot=100000,radius=radius,nradial=nradial,lum=lum)
#test.add_gaussian_variations(0.2)

#write input filepwd

test.write_input(mrw=True)
test.add_lines()
#run model
test.calculate_model(ncores=20)

test.make_vtk()

test.make_cube(nlam=nlam, ncores=20)

dir = f'Dens={dens:5.2e}_Lum={lum:5.2e}_ri={ri}_ro={ro}_radius={radius}_prho={prho}_nradial={nradial}_nlam={nlam}'
fitsfile = f'{dir}.fits'
test.make_fits(fitsfile)

cmd = f'mkdir {dir}; cp *.inp *.dat *.out *vtk {dir}; mv *.fits {dir}'

os.system(cmd)

# write code to read a file

#sed
#test.sed()

#wls=[450,850,1000,2000,3000] #micron
#these functions take some time and could benefit from multithreading too.
#test.make_synth_maps(wls)
#test.make_tau_surface(wls)
