#!/usr/bin/env python
from math import *
from numpy import *
from pylab import *

class Flow():
    def __init__(self, omega = 240):
        # Parameters here are written in the non dimensionless form
        self.Tc = 1
        self.a_s = 1
        self.alpha = 1
        self.k0 = 1
        self.n_vort = 2
        self.L = 6*(self.n_vort-1)
        self.u0 = 19.9066 
        self.uE = 2.
        self.delta = 0. 
        self.y0 = 0.5
        self.omega = omega
        self.D = 0.041472
        
    def flow_more_new(self, x, y, t, init_E, Nx, Ny):
        sigma = 1./(2*self.k0)

        R = sqrt(x**2 + y**2)
        f = 1. - exp(-self.a_s*(sqrt(x**2+y**2)-1)**2)
        s = 1. - exp((-(x-1)**2)/self.alpha**2-y**2)


        dys = (-s+1)*2*y
        dfdy = (1-f)*(2*self.a_s*y*(R-1)/R)
        dfdx=(1-f)*(2*self.a_s*x*(R-1)/R)
        dxs=(-s+1)*(2*(x-1)/(self.alpha**2))

        x_p = [1 + self.L*abs(mod(t/self.Tc, 1))]
        y_p = [self.y0 + self.delta*abs(mod(t/self.Tc, 1))]
        h_p = [abs(sin(pi*t/self.Tc))]
        sigma_p = [sigma +  self.delta*abs(mod(t/self.Tc, 1))]
        k_p = [1./(2*sigma_p[0])]
        g_p = [exp(-k_p[0]*((x-x_p[0])**2 + self.alpha*(y-y_p[0])**2))]
        dgdy_p = [g_p[0]*(-k_p[0]*self.alpha*2*(y-y_p[0]))]
        dgdx_p = [g_p[0]*(-k_p[0]*2*(x-x_p[0]))]

        for i in range(1, self.n_vort):
            x_p.append(1 + self.L*abs(mod((t - (self.Tc*i)/self.n_vort)/self.Tc, 1)))
            y_p.append((-1)**(i)*(self.y0 + self.delta*abs(mod((t - (self.Tc*i)/self.n_vort)/self.Tc, 1))))
            sigma_p.append(sigma + self.delta*abs(mod((t - (self.Tc*i)/self.n_vort)/self.Tc, 1)))
            k_p.append(1./(2*sigma_p[i]))
            h_p.append(abs(sin(pi*(t - (self.Tc*i)/self.n_vort)/self.Tc)))
            g_p.append(exp(-k_p[i]*((x-x_p[i])**2 + self.alpha*(y-y_p[i])**2)))
            dgdy_p.append(g_p[i]*(-k_p[i]*self.alpha*2*(y-y_p[i])))
            dgdx_p.append(g_p[i]*(-k_p[i]*2*(x-x_p[i])))        


        g = 0
        dgdy = 0
        dgdx = 0

        for i in range(0, self.n_vort):
            g += ((-1)**(i+1))*self.omega*h_p[i]*g_p[i]
            dgdy += ((-1)**(i+1))*self.omega*h_p[i]*dgdy_p[i]
            dgdx += ((-1)**(i+1))*self.omega*h_p[i]*dgdx_p[i]

        g += self.u0*s*y
        dgdy += self.u0*(dys)*y + self.u0*s
        dgdx += self.u0*(dxs)*y


        g = reshape(g, (Nx, Ny))
        dgdx = reshape(dgdx, (Nx, Ny))
        aux_x = reshape(x, (Nx, Ny))
        g[:, init_E:] = g[:, init_E:] + self.uE*(aux_x[:, init_E:] - 1.0)
        dgdx[:, init_E:] = dgdx[:, init_E:] + self.uE

        g = g.flatten()
        dgdx = dgdx.flatten()

        vx = (dfdy*g + dgdy*f)
        vy = -(dfdx*g + dgdx*f)

        return vx, vy
    
    def advect(self, x_in, t, Nx, Ny, init_E):                

        size = int(len(x_in)/2)
        x = x_in[:size]
        y = x_in[size:]

        vx, vy = self.flow_more_new(x, y, t, init_E, Nx, Ny)

        return concatenate((vx.flatten(), vy.flatten()), axis=0)


