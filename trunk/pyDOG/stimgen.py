'''
This file contains methods for creating stimulus movies for tcloop classes
'''

import pylab as pl

class motherStim( object ):
    def __init__( self, t, x, y, contrast=1 ):
        self.t = t
        self.x = x
        self.y = y
        self.contrast = float(contrast)
        self._xx, self._yy = pl.meshgrid(self.x, self.y)
        self.make_frames()

    def _frame( self, tpoint):
        s = pl.zeros(self._xx.shape)
        return s
        
    def make_frames( self ):
        self.frames = pl.zeros(list(self._xx.shape)+[self.t.size])
        for tidx, tpoint in enumerate(self.t):
            self.frames[:,:,tidx] = self._frame(tpoint)


class slidingBar( motherStim ):
    def __init__( self, theta, speed, width, *args, **kwargs ):
        
        self.theta = float(theta)
        self.speed = float(speed)
        self.width = float(width)
#        y_width = self.width/pl.sin(self.theta)
        motherStim.__init__( self, *args, **kwargs)

    def _frame( self, tpoint):
        s = pl.zeros(self._xx.shape)
        x0 = self.x[0] + tpoint * self.speed
        x1 = x0 + self.width
        
        s[:, self.x[self.x<=x1]>=x0] = 1

        return s        

class flashingSpots( motherStim ):
    def __init__( self, posx, posy, d, t0, duration, *args, **kwargs ):
        self.posx = posx
        self.posy = posy
        self.d = [float(entry) for entry in d]
        self.t0 = t0
        self.duration = duration
        
        motherStim.__init__( self, *args, **kwargs)

    def make_frames( self ):
        s = pl.zeros(list(self._xx.shape)+[self.t.size])
        for ispot in xrange(len(self.t0)):
            posx = self.posx[ispot]
            posy = self.posy[ispot]
            d = self.d[ispot]
            t0 = self.t0[ispot]
            duration = self.duration[ispot]
            s[pl.sqrt((self._xx-posx)**2 + (self._yy-posy)**2)<d/2.,
              :] = 1
            s[:,:,self.t<t0]=0
            s[:,:,self.t>t0+duration]=0
        
        self.frames = s
        

class patchGrating( motherStim ):
    def __init__(self):
        pass

class driftingGrating( motherStim ):
    def __init__(self):
        pass
