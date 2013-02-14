import os
import pylab as pl
import glob
from scitools import easyviz as evz

class DOG_plots(object):
    def __init__( self, t, x, y, casename, savepath, mvi_name):
        self.casename = casename
        self.savepath = savepath
        self.mvi_name = mvi_name
        self.t = t
        self.x = x
        self.y = y
        
    def plot_allpanels( self, stim, ret_r, ff_r, fb_r,
                        im_kwargs=[{},{},{},{}] ):
        fig = pl.figure()
        
        pl.subplot(221)
        extent = im_kwargs[0].pop('extent', (self.x[0], self.x[-1],
                                             self.y[0], self.y[-1]))
        if stim.max()==0:
            stim[0,0]=0.
        pl.imshow(stim, aspect='equal', origin='lower', extent=extent,
                  **im_kwargs[0])
        pl.title('Stimulus')
        
        pl.subplot(222)
        extent = im_kwargs[1].pop('extent', (self.x[0], self.x[-1],
                                             self.y[0], self.y[-1]))
        pl.imshow(ret_r, aspect='equal', origin='lower', extent=extent,
                  **im_kwargs[1])
        pl.title('Retinal response')
        
        pl.subplot(223)
        extent = im_kwargs[2].pop('extent', (self.x[0], self.x[-1],
                                             self.y[0], self.y[-1]))
        pl.imshow(ff_r, aspect='equal', origin='lower', extent=extent,
                  **im_kwargs[2])
        pl.title('LGN response, no feedback')
        pl.xlabel('x-pos (deg)')
        pl.ylabel('y-pos (deg)')
        
        pl.subplot(224)
        extent = im_kwargs[3].pop('extent', (self.x[0], self.x[-1],
                                             self.y[0], self.y[-1]))
        pl.imshow(fb_r, aspect='equal', origin='lower', extent=extent,
                  **im_kwargs[3])
        pl.title('LGN response, with feedback')
        
        return fig

    def plot_responses_allpanels( self, stims, retina_r, lgn_r_ff, lgn_r_fb,
                                  im_kwargs=[{},{},{},{}]):
    
        for f in glob.glob(os.path.join(self.savepath,
                                        self.casename + '_*.png')):
            os.remove(f)
        rcdefault = {'figure.figsize' : pl.rcParams['figure.figsize']}
        rcdict = {'figure.figsize':[14., 14.]}
        pl.rcParams.update(rcdict)
        for itime, tpoint in enumerate(self.t):
            if mod(itime, 50):
                print 'itime = %d' % itime
            savename = self.casename + '_%04d.png' % itime
            stim = stims[:, :, itime]
            ret_r = retina_r[:,:,itime]
            ff_r = lgn_r_ff[:,:,itime]
            fb_r = lgn_r_fb[:,:,itime]
            
            fig = self.plot_allpanels(stim, ret_r, ff_r, fb_r, 
                                      im_kwargs)
            pl.savefig(os.path.join(self.savepath, savename), 
                       format='png')
            pl.close(fig)
        pl.rcParams.update(rcdefault)
