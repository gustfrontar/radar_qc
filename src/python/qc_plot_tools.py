

def plot_dealiasing( qc_output , options , elev=0 , figname='out.png',vmin=-30,vmax=30,cmap='pyart_NWSVel')  :

    import numpy as np
    import matplotlib.pyplot as plt
  

    filter_name='Dealiasing'

    plt.figure()
    plt.subplot(2,2,1)

    plt.pcolor(qc_output['x'][:,:,elev],qc_output['y'][:,:,elev],qc_output['cv'][:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Corrected Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.pcolor(qc_output['x'][:,:,elev],qc_output['y'][:,:,elev],qc_output['v'][:,:,elev],vmin=vmin,vmax=vmax,cmap=cmap)
    plt.title('Original Doppler Velocity')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.pcolor(qc_output['x'][:,:,elev],qc_output['y'][:,:,elev], ( qc_output['qcv'][:,:,elev]==options[filter_name]['texture_code'] ).astype(float) )
    plt.title('Pixels eliminated by texture texture filter')
    plt.colorbar()

    plt.subplot(2,2,4)
    plt.pcolor(qc_output['x'][:,:,elev],qc_output['y'][:,:,elev], ( qc_output['qcv'][:,:,elev]==options[filter_name]['code'] ).astype(float) )
    plt.title('Pixels where aliasing was corrected')
    plt.colorbar()

    plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)


    return
