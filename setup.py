from model_setup_lines import model_setup_lines as model
import astropy.constants as const

rho0=1e6 *1.6735575e-24   #reference density cgs units
prho=4.0 #power law index
test=model(rho0,prho,size=5000,nphot=100000,radius=5000)
#test.add_gaussian_variations(1.0)

#write input file
test.write_input(mrw=True)
test.add_lines()
#run model
test.calculate_model(ncores=20)

test.make_vtk()

test.make_cube(ncores=20)

test.make_fits('1e6')
# write code to read a file

#sed
#test.sed()

#wls=[450,850,1000,2000,3000] #micron
#these functions take some time and could benefit from multithreading too.
#test.make_synth_maps(wls)
#test.make_tau_surface(wls)
