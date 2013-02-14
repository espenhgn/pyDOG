#!/usr/bin/env python
'''

'''
import pylab as pl
import numpy.fft as fft

class thalamoCorticalLoop( object ):
    def __init__( self, t, x, y, A=[1., 0.85], a=[0.25, 0.85] ):
        '''
        
        Initializing generalized thalamo-cortical loop. Full
        functionality is only obtained in the subclasses, like DOG,
        full_eDOG, etc.
        
        Parameters
        ----------
        
        t : array
            1D Time vector
        x : np.array
            1D vector for x-axis sampling points
        y : np.array
            1D vector for y-axis sampling points

        Keyword arguments
        -----------------
        
        A : sequence (default A = [1., 0.85])
            Amplitudes for DOG receptive field for center and surround, 
            respectively
        a : sequence (default a = [0.25, 0.85])
            Width parameter for DOG receptive field for center and surround,
            respectively

        Usage
        -----
            Look at subclasses for example usage            
        '''
        # Set parameteres as attributes
        self.name = 'pyDOG Toolbox'
        self.t = t
        self.A = A
        self.a = a
        self.x = x
        self.y = y
        # Find sampling rates and sampling freqs and 
        self.nu_xs = 1./(x[1]-x[0])
        self.nu_ys = 1./(y[1]-y[0])
        self.fs = 1./(t[1]-t[0])
        self.f = fft.fftfreq(pl.asarray(t).size, t[1]-t[0])
        # fftshift spatial frequency,
        self.nu_x = fft.fftfreq(pl.asarray(x).size, x[1]-x[0])
        self.nu_y = fft.fftfreq(pl.asarray(y).size, y[1]-y[0])
        # Make meshgrids, may come in handy
        self._xx, self._yy = pl.meshgrid(self.x, self.y)
        self._nu_xx, self._nu_yy = pl.meshgrid(self.nu_x, self.nu_y)
        # r is needed for all circular rfs
        self.r = pl.sqrt(self._xx**2 + self._yy**2)
        self.k = 2 * pl.pi * pl.sqrt(self._nu_xx**2 + self._nu_yy**2)

    def impulseResponse( self, domain='fourier', **kwargs):
        '''This function returns the impulse response for the TC loop'''
        G_fourier = self._iR( **kwargs )
        if domain=='fourier':
            ret = G_fourier
        elif domain=='spacetime':
            G = fft.ifftn(G_fourier)
            ret = G.real
 #           ret = self._shift_space(G)

        return ret

    def _iR(self, *args, **kwargs):
        pass
    
    def stimulusResponse( self, stimulus, domain='spacetime', **kwargs):
        '''
        
        '''
        stim_f = fft.fftn(stimulus)
        _ = kwargs.pop(domain, '')
        G = self.impulseResponse(domain='fourier', **kwargs)
        R_fourier = G * stim_f
        if domain=='fourier':
            R = R_fourier
        elif domain=='spacetime':
            R_tmp = fft.ifftn(R_fourier)
            R = R_tmp.real
#            R = self._shift_space(R_tmp)

        return R        

    def ff_kernel( self, domain='f' ):
        pass
    
    def ff_kernel_numeric( self ):
        '''

        '''
        h_t = self.ff_kernel( domain='t' )
        h = fft.fft(h_t)

        return h
    
    def fb_kernel( self ):
        pass

    def f_DOG( self, domain='k'):
        ''' '''
        A = self.A
        a = self.a
        r = self.r
        f_r = A[0] * _singleGaussian(r, a[0]) - A[1] * _singleGaussian(r, a[1])
#        r_sq = self._xx**2 + self._yy**2
#        f_r = A[0]*pl.exp( -r_sq/a[0]**2 )/( pl.pi*a[0]**2 ) \
#            - A[1]*pl.exp( -r_sq/a[1]**2 )/( pl.pi*a[1]**2 )        
        
        if domain=='r':
            f = f_r
        elif domain=='k':
            k = self.k
            f = A[0] * _singleGaussian_k(k, a[0]) \
                - A[1] * _singleGaussian_k(k, a[1])
#            k_sq = 4 * pl.pi**2 * (self._nu_xx**2 + self._nu_yy**2)
#            f = A[0]*pl.exp( -k_sq*a[0]**2/4 ) - A[1]*pl.exp( -k_sq*a[1]**2/4 )
        else:
            msg = 'Domain %s is not recognized. Choose *domain = "r"* for '\
                ' space domain or *domain = "k"* for fourier domain'
            raise Exception, msg
        
        return f

    def f_DOG_numeric( self ):
        '''

        '''
        f_r = self.f_DOG(domain='r')
        f = fft.fft2(f_r)
        
        return f

    def _shift_space( self, R ):
        R_out = pl.zeros(R.shape)
        for itime, _ in enumerate(self.t):
            R_out[:,:,itime].real = pl.ifftshift(R[:,:,itime].real)

        return R_out
        

class DOG( thalamoCorticalLoop ):
    def __init__( self, tau, Delta, *args, **kwargs ):
        self.tau = tau
        self.Delta = Delta
        thalamoCorticalLoop.__init__( self, *args, **kwargs )
        
    def _iR( self, f_mode='analytic', h_mode='analytic' ):

        if f_mode=='analytic':
            f = self.f_DOG( )
        elif f_mode=='numeric':
            f = self.f_DOG_numeric( )
        if h_mode=='analytic':
            h = self.ff_kernel( )
        elif h_mode=='numeric':
            h = self.ff_kernel_numeric( )
        G_fourier = f[:, :, pl.newaxis] * h[pl.newaxis, pl.newaxis, :]

        return G_fourier

    def ff_kernel( self, domain = 'f'):
        '''
        
        '''
        h_t = _singleExp(self.t, self.tau, self.Delta)
        
        if domain=='f':
            h = _singleExp_f(self.f, self.tau, self.Delta)

        elif domain=='t':
            h = h_t
        else:
            msg = 'Domain %s is not recognized. Choose *domain = "t"* for '\
                ' time domain or *domain = "f"* for frequency domain'
            raise Exception, msg

        return h

class full_eDOG( thalamoCorticalLoop ):
    def __init__( self, gs_params, rg_params, rig_params, fb_params,
                  *args):
        '''
        
        Initializing full eDOG loop, as expressed in Einevoll and
        Plesser (2012):

            
        
        Parameters
        ----------
        
        t : array
            1D Time vector
        x : np.array
            1D vector for x-axis sampling points
        y : np.array
            1D vector for y-axis sampling points

        Keyword arguments
        -----------------
        
        A : sequence (default A = [1., 0.85])
            Amplitudes for DOG receptive field for center and surround, 
            respectively
        a : sequence (default a = [0.25, 0.85])
            Width parameter for DOG receptive field for center and surround,
            respectively

        Usage
        -----
            Look at subclasses for example usage            
        '''
        # Ganglion parameters
        A_gs = gs_params[0]
        a_gs = gs_params[1]
        self.tau_gs = gs_params[2][0]
        self.tau_igs = gs_params[2][1]
        self.D_gs = gs_params[3][0]
        self.D_igs = gs_params[3][1]
        # Exitatory ff parameters
        self.A_rg = rg_params[0]
        self.a_rg = rg_params[1]
        self.tau_rg = rg_params[2]
        self.D_rg = rg_params[3]
        # Inhibitory ff params
        self.A_rig = rig_params[0]
        self.a_rig = rig_params[1]
        self.tau_rig = rig_params[2]
        self.D_rig = rig_params[3]
        # Feedback params
        self.A_fb = fb_params[0]
        self.a_fb = fb_params[1]
        self.tau_fb = fb_params[2]
        self.D_fb = fb_params[3]
        
        thalamoCorticalLoop.__init__( self, *args, A = A_gs, a = a_gs)

    def impulseResponseRetina( self, domain='fourier'):
        '''
        This has a DOG for spatial profile and a biphasic temporal kernel
        '''
        # Use analytical version of DOG
        f = self.f_DOG( )
        # Numerical transform of temporal kernel
        h_t = self.h_gs( )
        h = fft.fft(h_t)
        
        G_fourier = f[:, :, pl.newaxis] * h[pl.newaxis, pl.newaxis, :]

        if domain=='fourier':
            G = G_fourier
        elif domain=='spacetime':
            G_tmp = fft.ifftn(G_fourier)
#            G = self._shift_space(G_tmp)

        return G

    def _iR( self, cell_response='lgn' ):
        '''
        
        '''
        if cell_response=='lgn':
            G_fourier = (self.K_rg() + self.K_rig())/(1 - self.K_fb()) \
                * self.impulseResponseRetina()
        elif cell_response=='retina':
            G_fourier = self.impulseResponseRetina()
        else:
            msg =  'only implemented "lgn" and "retina" as cell'\
                +' responses'
            raise Exception, msg

        return G_fourier
        
    def h_gs( self ):
        '''
        This is the feedforward kernel describing processing in retina
        '''
        h1 = _singleAlpha(self.t, self.tau_gs, self.D_gs)
        h2 = _singleAlpha(self.t, self.tau_igs, self.D_igs)

        return h1-h2

    def K_rg( self ):
        K_rg = self.A_rg * self.K_gauss_exp(self.a_rg, self.tau_rg, self.D_rg)
        return K_rg
    
    def K_rig( self ):
        K_rig = self.A_rig * self.K_gauss_exp(self.a_rig, self.tau_rig,
                                              self.D_rig)
        return K_rig

    def K_fb( self ):
        K_fb = self.A_fb * self.K_gauss_exp(self.a_fb, self.tau_fb, self.D_fb)
        return K_fb

    def K_gauss_exp( self, a, tau, D):        
        f = _singleGaussian_k(self.k, a)
        h_t = _singleExp(self.t, tau, D)
        h = fft.fft(h_t)

        K = f[:, :, pl.newaxis] * h[pl.newaxis, pl.newaxis, :]

        return K

class eDOG( thalamoCorticalLoop ):
    def __init__( self, C, c, *args,  **kwargs ):
        self.C = C
        self.c = c
        thalamoCorticalLoop.__init__( self, *args, **kwargs )

    def impulseResponseRetina( self, domain='fourier' ):
        f = self.f_DOG()
        h_t = self.ff_kernel()
        h = fft.fft(h_t)

        G_fourier = f[:, :, pl.newaxis] * h[pl.newaxis, pl.newaxis, :]

        if domain=='fourier':
            G = G_fourier
        elif domain=='spacetime':
            G = fft.ifftn(G_fourier)

        return G    

    def _iR( self, cell_response='lgn' ):
        '''
        
        '''
        if cell_response=='lgn':
            G_fourier = self.impulseResponseRetina()/(1-self.K_fb())
        elif cell_response=='retina':
            G_fourier = self.impulseResponseRetina()
        else:
            msg =  'only implemented "lgn" and "retina" as cell'\
                +' responses'
            raise Exception, msg

        return G_fourier

    def K_fb( self ):
        f = _singleGaussian_k(self.k, self.c)
        h_fb_t = self.fb_kernel()
        h_fb = fft.fft(h_fb_t)
        
        K_fb = self.C * f[:,:,pl.newaxis] * h_fb[pl.newaxis, pl.newaxis, :]
        
        return K_fb
        
    def ff_kernel( self ):
        h1 = _singleAlpha(self.t, self.tau_gs, self.D_gs)
        h2 = _singleAlpha(self.t, self.tau_igs, self.D_igs)

        return h1-h2

    def fb_kernel( self ):
        h_t = _singleExp(self.t, self.tau_fb, self.D_fb)
        return h_t
        

class eDOG_FastLoop( eDOG ):
    def __init__( self, *args, **kwargs ):
        eDOG.__init__( self, *args, **kwargs )


def _singleGaussian(r, a):
    f = 1/(pl.pi*a**2)*pl.exp(-r**2/a**2)
    return f

def _singleGaussian_k(k, a):
    f = pl.exp(-k**2*a**2/4)
    return f

def _singleAlpha(t, tau, Delta):
    h = (t-Delta)/tau**2*pl.exp(-(t-Delta)/tau)
    h[t<Delta]=0

    return h

def _singleExp(t, tau, Delta):
    # Create the convolution kernel
    h = 1./tau*pl.exp(-(t-Delta)/tau)
    h[t<Delta]=0

    return h

def _diffExp(t, tau_1, Delta_1, tau_2, Delta_2):
    # Create difference of exponentials
    h_1 = _singleExp(t, tau_1, Delta_1)
    h_2 = _singleExp(t, tau_2, Delta_2)

    h = h_1 - h_2

    return h

def _singleExp_f(f, tau, Delta):
    # Create convolution kernel in frequency space
    h_f = pl.exp(-2*pl.pi*f*Delta*1j)/(2*pl.pi*f*tau*1j + 1)

    return h_f
    
if __name__ == '__main__':
    print 'Do something meaningful here'
    
