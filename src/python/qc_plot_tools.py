

def plot_dealiasing( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt
  

    filter_name='Dealiasing'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['texture_code'] ).astype(float) )
    plt.title('Pixels eliminated by texture texture filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels where aliasing was corrected')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_rhofilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='RhoFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)
                            
    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by RHO Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['rho_smooth'][:,:,elev] ),vmin=0,vmax=1.1 )
    plt.title('Smooth RHO HV')
    plt.colorbar()
                                                                                                    
    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return



def plot_echotopfilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='EchoTopFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Echo Top Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['echo_top'][:,:,elev] ),vmin=0,vmax=15000 )
    plt.title('Echo Top')
    plt.colorbar()
													
    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_echodepthfilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='EchoDepthFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Echo Depth Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['echo_top'][:,:,elev] ),vmin=0,vmax=15000 )
    plt.title('Echo Depth')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_refspecklefilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='RefSpeckleFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Ref. Speckle Filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['speckle_ref'][:,:,elev] ),vmin=0,vmax=1 )
    plt.title('Speckle Index')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return


def plot_dopplerspecklefilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerSpeckleFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, qc_output['speckle_v'][:,:,elev] , vmin=0,vmax=1)
    plt.title('Speckle Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Doppler Speckle Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return



def plot_dopplertexturefilter( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='DopplerTextureFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_v=np.ma.masked_array( qc_output['v'] , qc_output['v'] == qc_output['undef_v'] )
    tmp_cv=np.ma.masked_array( qc_output['cv'] , qc_output['cv'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_cv[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3,tmp_v[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, qc_output['texture_v'][:,:,elev] , vmin=0 , vmax=50 )
    plt.title('Texture Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels eliminated by Texture Filter')
    plt.colorbar()

    if show  :

        plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return


def plot_reflectivitytexturefilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='ReflectivityTextureFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['texture_ref'][:,:,elev] ) )
    plt.title('Texture Index')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Ref. Texture Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_attenuationfilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='AttenuationFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['attenuation'][:,:,elev] ) )
    plt.title('Attenuation')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Attenuation Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return

def plot_blockingfilter( qc_output , options , elev=0 , figname='out.png',vmin=-10,vmax=70,cmap='pyart_NWSRef',show=True)  :

    import numpy as np
    import matplotlib.pyplot as plt


    filter_name='BlockingFilter'

    plt.figure(figsize=(8, 8))
    plt.subplot(2,2,1)

    tmp_ref=np.ma.masked_array( qc_output['ref'] , qc_output['ref'] == qc_output['undef_v'] )
    tmp_cref=np.ma.masked_array( qc_output['cref'] , qc_output['cref'] == qc_output['undef_v'] )

    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_cref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, tmp_ref[:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Reflectivity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['blocking'][:,:,elev] ) )
    plt.title('Blocking')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev]/1e3,qc_output['y'][:,:,elev]/1e3, ( qc_output['qcref'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels Eliminated by Blocking Filter')
    plt.colorbar()

    if show  :
       plt.show()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    return




