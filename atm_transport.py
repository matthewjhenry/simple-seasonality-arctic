import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import xarray as xr

def my_integration(f_,sigma_full_):

    # define sigma on half levels
    num_lev_    = sigma_full_.shape[0]
    sigma_half_ = np.zeros(num_lev_+1)

    #  lowest half level is ground
    sigma_half_[num_lev_] = 1
    
    sigma_half_[0]=0

    for j in range(num_lev_-1,-1,-1):
        sigma_half_[j] = 2*sigma_full_[j] - sigma_half_[j+1]
    
    sigma_half_[0]=0 
    d_sigma_  = np.diff(sigma_half_)
    d_temp_   = np.transpose(d_sigma_).reshape(f_.shape[0],1)
    d_sigma2_ = np.matlib.repmat(d_temp_,1,f_.shape[1])
    
    fin=np.sum(f_*d_sigma2_,0)
    
    return fin

def press_int_k(f, plev_full):
    plev_full_0 = np.concatenate([[0], plev_full],axis=0)
    d_plev = np.diff(plev_full_0)/100
    d_temp_ = np.transpose(d_plev).reshape(f.shape[0],1)
    d_plev2 = np.matlib.repmat(d_temp_, 1, f.shape[1])
    F      = np.sum(f*d_plev2,0)
    return F

def energy_flux(xar_obj):
    
    # set parameters
    del_sol         = 1.0
    start_xtropics  = 25.0
    deg             = np.pi/180
    l_cond          = 2500000
    radius          = 6371000
    p0              = 100000
    grav            = 9.8000
    cp              = 1.00464e+03
    
    lati, sigi = np.meshgrid(xar_obj.lat,xar_obj.pfull/1000)
    temp_flux = xar_obj.vcomp_temp.mean(('months','lon'))/np.cos(lati*deg)
    z_flux = xar_obj.vcomp_height.mean(('months','lon'))/np.cos(lati*deg)
    sh_flux = xar_obj.sphum_v.mean(('months','lon'))/np.cos(lati*deg)
    
    # variables to be plotted
    energy_flux =  cp*temp_flux + grav*z_flux + l_cond*sh_flux
    dry_energy_flux =  cp*temp_flux+grav*z_flux
    hum_energy_flux =  l_cond*sh_flux
        
    # now integrate over sigma
    energy_flux_1d = sigma_int.my_integration(energy_flux, sig)
    dry_energy_flux_1d = sigma_int.my_integration(dry_energy_flux, sig)
    hum_energy_flux_1d = sigma_int.my_integration(hum_energy_flux, sig)
    # z_energy_flux_1d = sigma_int.my_integration(z_energy_flux, sig)
    
    # convert to power in PW
    rescale = 2.0*pi*radius*p0/grav*np.cos(lat*deg)/1e15
    energy_flux_1d = energy_flux_1d*rescale
    dry_energy_flux_1d = dry_energy_flux_1d*rescale
    hum_energy_flux_1d = hum_energy_flux_1d*rescale
    # z_energy_flux_1d = z_energy_flux_1d*rescale
    
    return energy_flux_1d #,dry_energy_flux_1d,hum_energy_flux_1d #,z_energy_flux_1d


