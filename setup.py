from model_setup_lines import model_setup_lines as model
import astropy.constants as const
import os

sdens = "1e6"

dens = float(sdens)


prho=4.0 #power law index
test=model(dens,prho,size=10000,nphot=100000,radius=5000)
#test.add_gaussian_variations(1.0)

#write input file
test.write_input(mrw=True)
test.add_lines()
#run model
test.calculate_model(ncores=20)

test.make_vtk()

test.make_cube(ncores=20)

fitsfile = f'{sdens}.fits'
test.make_fits(fitsfile)

cmd = f'mkdir {sdens}; cp *.inp *.fits *.dat *.out *.fits *vtk {sdens}'

os.system(cmd)

# write code to read a file

#sed
#test.sed()

#wls=[450,850,1000,2000,3000] #micron
#these functions take some time and could benefit from multithreading too.
#test.make_synth_maps(wls)
#test.make_tau_surface(wls)
