# Based on work from Tom Carron, modified by Peter Schilke


import astropy.constants as const
import astropy.units as u
import numpy as np
import math
import os
import subprocess
import matplotlib.pyplot as plt
import time
import matplotlib.pylab as plb
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from spectral_cube import SpectralCube
from spectral_cube import Projection
import radio_beam
from radio_beam import Beam
from astropy.io import fits
from radmc3dPy import image, analyze, natconst    
#from radmc3dPy.analyze import *  
#from radmc3dPy.natconst import * 
from astropy.wcs import WCS
import sys
#sys.path.append('../../SED_sgrb2/SED_fit')
#from fit import extract_dimensions

'''3D model with spherical coordinates'''

class model_setup_lines:
    def __init__(self,dens,prho,ri=50,ro=5000,nphot=100000,radius=3000,nradial=100,lum=1e4):

        au=(const.au).to("cm").value

        self.nphot=nphot

        self.ri = ri * au
        self.ro = ro * au
 
        sigma = const.sigma_sb.value
        Lsun = const.L_sun.value
        Rsun = const.R_sun.value
        Msun = const.M_sun.value
        G = const.G.value
        #
        # Star parameters
        #
        radius_star = 13.4
        mass_star = 30
        self.mstar    = mass_star*(const.M_sun).to("g").value # in kg, check if correct units
        self.rstar    = radius_star*(const.R_sun).to("cm").value  # in m , check units
        self.tstar    = (lum*Lsun/(4*np.pi*sigma*(radius_star*Rsun)**2))**0.25 
        self.ls       = (const.L_sun).to("erg/s").value  #solar luminosity
        self.pstar    = np.array([0.,0.,0.]) #position in cartesian coords
        #
        Lum = sigma*4*np.pi*(radius_star*Rsun)**2*self.tstar**4/Lsun 

        print (f'Luminosity: {Lum:5.2e}, Temperature: {self.tstar:5f}' )
        # Wavelengths - this eventually needs a function to calculate it based of start and endpoint and maybe number of intervals.
        #
        lam1     = 0.01e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        lam5	 = 9.0e4
        n12      = 100
        n23      = 100
        n34      = 100
        n45	     = 100
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=False)
        lam45    = np.logspace(np.log10(lam4),np.log10(lam5),n45,endpoint=True)
        self.lam      = np.concatenate([lam12,lam23,lam34,lam45])
        self.nlam     = self.lam.size

        # grid parameters - eventually replace this with a function or put as argument.
        self.nr = nradial  # 
        self.ntheta = 90  # 
        self.nphi = 180  # 
        #

        #
  
        self.ri       = np.logspace(np.log10(self.ri),np.log10(self.ro), self.nr+1)
        self.thetai   = np.linspace(0, np.pi,self.ntheta+1)
        self.phii     = np.linspace(0, 2*np.pi,self.nphi+1)
        self.rc       = 0.5e0 * ( self.ri[:-1] + self.ri[1:] )
        rdiff = np.diff(self.ri)
        rtheta       = 0.5e0 * ( self.thetai[:-1] + self.thetai[1:] )
        rphi       = 0.5e0 * ( self.phii[:-1] + self.phii[1:] )

        self.rho0 = dens * 2.3 * natconst.mp
        self.prho     = prho   


        #
        # Make the dust density model
        #
        self.radius = radius * au
        rr, tt, pp      = np.meshgrid(self.rc, rtheta, rphi,indexing='ij')
        self.rhod     = self.rho0 /(1.+rr**2/self.radius**2)**(self.prho/2)
        densr = dens/(1.+self.rc**2/self.radius**2)**(self.prho/2)
    #    print (self.rhod.shape, rr.shape)
#        for i in range(rr.shape[0]):
#            print (self.rhod[i], rr[i])
       #
        # Make the dust density model in the layer
        cold = rdiff*densr
        print (f'total column density {cold.sum(): 5.2e}')
   
    def add_gaussian_variations(self, std_deviation):
        # Generate Gaussian variations for each voxel
        gaussian_variations = np.random.normal(0, std_deviation, size=(self.nr, self.ntheta, self.nphi))

        # Apply variations to the density
        self.rhod *= 10**gaussian_variations

    def add_lines(self, vin=0, Tlow=50, Thigh=1000, abunch3cn=1e-8):
        
        # Model parameters
        #
        dusttogas= 0.01
        vturb0   = 3.*1e5

        vturb   = np.zeros((self.nr,self.ntheta,self.nphi)) + vturb0
        data = analyze.readData(dtemp=True)
        dt = data.dusttemp[:,:,:,0]
        print (dt.shape)
 
        factch3cn = abunch3cn/(2.3*natconst.mp)
        nch3cn    = self.rhod*factch3cn
        nch3cn[dt>Thigh] = 0
        nch3cn[dt<Tlow] = 0

        with open('numberdens_ch3cn.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(self.nr*self.ntheta*self.nphi))           # Nr of cells
            data = nch3cn.ravel(order='F')          # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
        #
        # Write the gas velocity field  - right now no velocity
        #
        with open('gas_velocity.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(self.nr*self.ntheta*self.nphi))           # Nr of cells
            for iphi in range(self.nphi):
                for itheta in range(self.ntheta):
                    for ir in range(self.nr):
                        ri = self.rc[ir]
                        vr = -vin*np.sqrt(2*const.G.cgs.value*self.mstar/ri)
                        f.write('%13.6e %13.6e %13.6e\n'%(vr,0,0))
        #
        # Write the microturbulence file
        #
        with open('microturbulence.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(self.nr*self.ntheta*self.nphi))           # Nr of cells
            data = vturb.ravel(order='F')          # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

        #
        # Write the lines.inp control file
        #
        with open('lines.inp','w') as f:
            f.write('1\n')
            f.write('1\n')
            f.write('ch3cn    leiden    0    0\n')
    

    def write_input(self,amr=False,mrw=False):
        # Write the wavelength file
    
        with open('wavelength_micron.inp','w+') as f:
            f.write('%d\n'%(self.nlam))
            for value in self.lam:
                f.write('%13.6e\n'%(value))
        #
        #
        # Write the stars.inp file
        #
        with open('stars.inp','w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n'%(self.nlam))
            f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(self.rstar,self.mstar,self.pstar[0],self.pstar[1],self.pstar[2]))
            for value in self.lam:
                f.write('%13.6e\n'%(value))
            f.write('\n%13.6e\n'%(-self.tstar))
        #
        # Write the grid file
        #
        with open('amr_grid.inp','w+') as f:
            f.write('1\n')                       # iformat
            f.write('0\n')                      # AMR grid style  (10=layer-style AMR)
            f.write('100\n')                       # Coordinate system
            f.write('0\n')                       # gridinfo
            f.write('1 1 1\n')                   # Include x,y,z coordinate
            f.write('%d %d %d\n'%(self.nr,self.ntheta,self.nphi))     # Size of grid/usr/bin/time 
            for value in self.ri:
                f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
            for value in self.thetai:
                f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
            for value in self.phii:
                f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)


        #
        # Write the density file
        #
        with open('dust_density.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(self.nr*self.ntheta*self.nphi)) # Nr of spatial data points (incl redundant ones)
            f.write('1\n')                       # Nr of dust species
            data = self.rhod.ravel(order='F')         # Create a 1-D view of rhod
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

        #
        # Dust opacity control file
        #
        with open('dustopac.inp','w+') as f:
            f.write('2               Format number of this file\n')
            f.write('1               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('1               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('silicate        Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')
        #
        # Write the radmc3d.inp control file
        #
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(self.nphot))
            f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
            f.write('iranfreqmode = 1\n')
            f.write('tgas_eq_tdust = 1\n')

            if mrw:
                f.write('modified_random_walk = 1')

    def calculate_model(self,ncores=None):
        t0 = time.time()
        if ncores is None:
            out = subprocess.run(["radmc3d",  "mctherm"], capture_output=True, text=True)
        else:
            out = subprocess.run(['radmc3d', 'mctherm',  'setthreads',  f'{ncores}'], capture_output=True, text=True)
        #os.system('radmc3d sed')
        t1 = time.time()

        total = t1-t0
        print(f'Calculating the model cost: {total}')
        with open('cost.out','w+') as f:
            f.write(f'Calculating the model cost: {total}')
        #Make the necessary calls to run radmc3d
        return out

    def make_vtk(self):
        os.system('radmc3d vtk_dust_temperature 1 vtk_dust_density 1')
   
    def make_cube(self, lambda1 = 1358.0, lambda2=1361.5, nlam=300, ncores=20):
        cmd = f'radmc3d image lambdarange {lambda1} {lambda2} nlam {nlam} setthreads {ncores}'
        os.system(cmd)

    def make_circular_image(self, wavel=10):
        cmd = f'radmc3d image circ lambda {wavel}'
        os.system(cmd)
    

    def make_fits(self, filename):
        im=image.readImage()
        im.writeFits(filename)
        
    def sed(self):
        #plot sed
        os.system('radmc3d sed')
        s    = readSpectrum()
        plt.figure()
        lammic = s[:,0]
        flux   = s[:,1]
        nu     = 1e4*const.c.cgs/lammic
        nufnu  = nu*flux
        nulnu  = nufnu*4*math.pi*(const.pc.cgs)*(const.pc.cgs)
        plt.plot(lammic,nulnu/self.ls)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\lambda$ [$\mu$m]')
        plt.ylabel(r'$\nu L_{\nu}$ [$L_{\odot}$]')
        plt.axis()
        plt.savefig('sed.png',dpi=200,bbox_inches='tight')
        #plt.show()
    
    def make_synth_maps(self,wls):
        t0 = time.time()
        for wl in wls:
            os.system('radmc3d image lambda '+str(wl))
            im   = image.readImage()
            im.writeFits('model'+str(wl)+'.fits')
            process_radmc_image('model'+str(wl)+'.fits','model'+str(wl)+'_smooth.fits',0.4,overwrite=True)
            cim = im.imConv(fwhm=[0.4, 0.4], pa=0., dpc=8340.) 
            #plt.figure()
            #plotImage(cim, arcsec=True, dpc=8340., cmap=plb.cm.gist_heat)
            #cim.writeFits('model'+str(wl)+'.fits', dpc=8340.)
            #plt.show()
        t1 = time.time()
        total=t1-t0
        print("Calculating the images cost: "+str(total)+" s")
        #plot synthetic maps at each wavelength specified.

    def make_tau_surface(self,wls):
        for wl in wls:
            os.system("radmc3d tausurf 1.0 lambda "+str(wl))
            im   = image.readImage()
            data = np.squeeze(im.image[:, ::-1, 0].T)
            #plotImage(im,arcsec=True,dpc=8340.)
            im.writeFits('tau_surf_'+str(wl)+'.fits')
            wcs = WCS(fits.getheader('tau_surf_'+str(wl)+'.fits'))
            newhdu = fits.PrimaryHDU(data,wcs.to_header())
            newhdu.writeto('tau_surf_'+str(wl)+'.fits',overwrite=True)

    def make_tau_map(self,wls):
        for wl in wls:
            os.system("radmc3d image lambda "+str(wl)+" tracetau")
            im   = image.readImage()
            data = np.squeeze(im.image[:, ::-1, 0].T)
            #plotImage(im,arcsec=True,dpc=8340.)
            im.writeFits('tau'+str(wl)+'.fits')
            wcs = WCS(fits.getheader('tau'+str(wl)+'.fits')).celestial
            newhdu = fits.PrimaryHDU(data,wcs.to_header())
            newhdu.writeto('tau'+str(wl)+'.fits',overwrite=True)

    def make_column_density_map(self,wls):
        for wl in wls:
            os.system("radmc3d image lambda "+str(wl)+" tracecolumn")
            im   = image.readImage()
            data = np.squeeze(im.image[:, ::-1, 0].T)
            #plotImage(im,arcsec=True,dpc=8340.)
            im.writeFits('column_density_'+str(wl)+'.fits')
            #open fits file, get wcs, overwrite
            wcs = WCS(fits.getheader('column_density_'+str(wl)+'.fits')).celestial
            newhdu = fits.PrimaryHDU(data,wcs.to_header())
            newhdu.writeto('column_density_'+str(wl)+'.fits',overwrite=True)

    
    def density_profile(self):
        #plot density vs radius
        au=const.au.cgs.value # AU [cm]
        a    = readData(ddens=True,binary=False)
        r    = a.grid.x[:]
        density = a.rhodust[:,0,0,0]
        plt.figure(1)
        plt.plot(r/au,density)
        plt.xlabel('r [au]')
        plt.ylabel(r'$\rho_{dust}$ [$g/cm^3$]')
        plt.show()
        plt.savefig('density.png')

'''
Support functions below
'''
def smooth_fits_image(input_file, output_file, target_resolution_major, target_resolution_minor):
    '''
    smooths an image to a target resolution.
    '''
    # Open the input FITS file
    hdul = fits.open(input_file)
    data = hdul[0].data

    # Get the current beam resolution
    current_resolution_major = hdul[0].header['BMAJ']  # Assumes major axis beam resolution is stored in BMAJ keyword
    current_resolution_minor = hdul[0].header['BMIN']  # Assumes minor axis beam resolution is stored in BMIN keyword

    # Compute the kernel width ratios
    sigma_ratio_major = target_resolution_major / current_resolution_major
    sigma_ratio_minor = target_resolution_minor / current_resolution_minor

    # Create a Gaussian kernel for smoothing
    kernel = Gaussian2DKernel([sigma_ratio_major, sigma_ratio_minor])

    # Convolve the data with the kernel
    smoothed_data = convolve(data, kernel)

    # Update the header with the new beam resolution
    hdul[0].header['BMAJ'] = target_resolution_major
    hdul[0].header['BMIN'] = target_resolution_minor

    # Save the smoothed data to a new FITS file
    hdul[0].data = smoothed_data
    hdul.writeto(output_file, overwrite=True)

    # Close the input FITS file
    hdul.close()


def process_radmc_image(input_fits, output_fits,beam_size_arcsec,overwrite=False):
    hdulist = fits.open(input_fits)
    data = extract_dimensions(hdulist[0].data)
    header = hdulist[0].header
    wcs = WCS(header)

    # Extract pixel size information from the WCS. Pixel sizes are given in degrees
    pixel_size_x, pixel_size_y = wcs.pixel_scale_matrix[1, 1], wcs.pixel_scale_matrix[0, 0]

    # Convert beam size from arcseconds to pixels
    beam_size_x = abs(int(beam_size_arcsec / 3600 / pixel_size_x))
    beam_size_y = abs(int(beam_size_arcsec / 3600 / pixel_size_y))

    #check if beam sizes are odd *must be odd for kernel
    if beam_size_x % 2 == 0:
        beam_size_x+=1
    if beam_size_y % 2 == 0:
        beam_size_y+=1

    # Calculate the standard deviation of the Gaussian kernel
    beam_stddev_x = beam_size_x / (2 * np.sqrt(2 * np.log(2)))

    print(beam_size_x,beam_size_y)
    # Create a 2D Gaussian kernel
    kernel = Gaussian2DKernel(beam_stddev_x, x_size=int(beam_size_x), y_size=int(beam_size_y))

    # Convolve the data with the Gaussian kernel
    smoothed_data = convolve(data, kernel, normalize_kernel=True)

    # Scale the entire image to convert pixel values from Jy/pixel to Jy/beam
    #conversion_factor = ((pixel_size_x/3600) * (pixel_size_y/3600) ) / (np.pi * beam_size_arcsec**2)
    #smoothed_data *= conversion_factor

    # Update the header to reflect the new units and beam information
    header['BUNIT'] = 'Jy/beam'
    header['BMAJ'] = beam_size_arcsec
    header['BMIN'] = beam_size_arcsec

    # Save the smoothed data to a new FITS file
    fits.writeto(output_fits, smoothed_data, header, overwrite=overwrite)

    print(f"Smoothing completed. Result saved to: {output_fits}")

def process_radmc_tausurf_image(input_fits, output_fits,overwrite=False):
    hdulist = fits.open(input_fits)
    data = extract_dimensions(hdulist[0].data)
    header = hdulist[0].header
    wcs = WCS(header)

    # Update the header to reflect the new units and beam information
    header['BUNIT'] = '[cm]'

    # Save the smoothed data to a new FITS file
    fits.writeto(output_fits, data, header, overwrite=overwrite)

    print(f" Result saved to: {output_fits}")
    

#function which returns the two largest dimensions of an array
def extract_dimensions(array):
    if array.ndim <= 2:
        return array
    else:
        dimensions_to_remove = np.where(np.array(array.shape) < 2)[0]
        modified_array = np.squeeze(array, axis=tuple(dimensions_to_remove))
        return modified_array
    

def read_image():
    im = image.readImage()
    #data = np.squeeze(im.image[:, ::-1, 0].T)
    data = im.image

