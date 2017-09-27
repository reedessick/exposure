__description__ = "holds things relating to antenna patterns"
__author__ = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import numpy as np
import healpy as hp

#-------------------------------------------------

def antenna_patterns(theta, phi, psi, nx, ny, freqs=None, dt=0.0, dr=None):
    """
    computes the antenna patterns for detector arms oriented along nx and ny (cartesian vectors). 
        if freqs, it computes time-shift phases in the frequency domain using dt. 
        if dr and freq, it will compute dt for itself (save time with cos(theta), etc.

    Antenna patterns are computed accoring to Eqn. B7 from Anderson, et all PhysRevD 63(04) 2003
    """
    n_pix = len(theta)

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)

    Xx = sin_phi*cos_psi - sin_psi*cos_phi*cos_theta
    Xy = -cos_phi*cos_psi - sin_psi*sin_phi*cos_theta
    Xz = sin_psi*sin_theta

    Yx = -sin_phi*sin_psi - cos_psi*cos_phi*cos_theta
    Yy = cos_phi*sin_psi - cos_psi*sin_phi*cos_theta
    Yz = sin_theta*cos_psi

    X = (Xx, Xy, Xz)
    Y = (Yx, Yy, Yz)

    ### iterate over x,y,z to compute F+ and Fx
    Fp = np.zeros((n_pix,),float)
    Fx = np.zeros((n_pix,),float)

    for i in xrange(3):
        nx_i = nx[i]
        ny_i = ny[i]
        Xi = X[i]
        Yi = Y[i]
        for j in xrange(3):
            Xj = X[j]
            Yj = Y[j]
            Dij = 0.5*(nx_i*nx[j] - ny_i*ny[j])
            Fp += (Xi*Xj - Yi*Yj)*Dij
            Fx += (Xi*Yj + Yi*Xj)*Dij

    ### apply time-shits
    if freqs != None:
        freqs = np.array(freqs)
        n_freqs = len(freqs)
        if dr != None:
            dx, dy, dz = dr
            dt = dx*sin_theta*cos_phi + dy*sin_theta*sin_phi + dz*cos_theta
            phs = 2*np.pi*np.outer(dt,freqs)
            phs = np.cos(phs) - 1j*np.sin(phs)
        else:
            phs = np.ones((n_pix,n_freqs),float)

        ones_freqs = np.ones((n_freqs),float)
        Fp = np.outer(Fp, ones_freqs) * phs
        Fx = np.outer(Fx, ones_freqs) * phs

    if n_pix == 1:
        return Fp[0], Fx[0]
    else:
        return Fp, Fx
