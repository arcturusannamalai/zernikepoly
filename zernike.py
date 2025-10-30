# (C) 2025, Muthu Annamalai, <muthu@panmo.cloud>
# This file is released under Apache 2.0 license
# If you find this program useful in your research work please cite the article:
# M. Annamalai, "Techniques to Improve Computation of Zernike Polynomials,"
# Frontiers in Optics - Laser Science, Denver, CO (2025).

import numpy as np
from math import factorial
from numpy import pi
from scipy.linalg import lstsq
from functools import cache
import time

def radial(n, m, rho):
    """
    Compute the Zernike radial polynomial R_n^m(ρ) for each value in ρ.
    Parameters
    ----------
    n : int
        Radial order (n ≥ 0)
    m : int
        Azimuthal order (|m| ≤ n, n−m even)
    rho : ndarray
        Array of radial coordinates (same shape)
    Returns
    -------
    R : ndarray
        Radial polynomial evaluated at each ρ
    """
    if (n - m) % 2:
        return np.zeros_like(rho)
    R = np.zeros_like(rho, dtype=float)
    k_max = (n - m) // 2
    for k in range(k_max + 1):
        num = (-1)**k * factorial(n - k)
        den = factorial(k)
        inv = 1.0 /(factorial((n + m)//2 - k) * factorial((n - m)//2 - k))
        R += (num/den) * inv * rho**(n - 2*k)
    return R

def zernike(n, m, rho, theta, fieldcache = None):
    """
    Compute the full Zernike polynomial Z_n^m(ρ,θ).
    """
    if fieldcache:
        if (n,m) in fieldcache:
            return fieldcache[(n,m)]
        if ((n-1,m-1) in fieldcache) and ((n-1,m+1) in fieldcache) and ((n-2,m) in fieldcache):
            # use the 3-term recurrence
            t1 = fieldcache[(n-1,m-1)]
            t2 = fieldcache[(n-1,m+1)]
            t3 = fieldcache[(n-2,m)]
            field = rho*(t1+t2) - t3
            fieldcache[(n, m)] = field
            return field

    R = radial(n, abs(m), rho)
    field = R * np.sin( (m >= 0)*(pi/2.0 - abs(m) * theta) + (m < 0)*(abs(m)*theta) )
    if fieldcache:
        fieldcache[(n,m)] = field
    return field

@cache
def generate_basis(n_max, resolution):
    """
    Generate all Zernike basis modes up to order n_max on a square grid.
    Parameters
    ----------
    n_max : int
        Maximum radial order.
    resolution : int
        Number of pixels along each axis (grid is resolution×resolution).
    Returns
    -------
    basis : list of ndarray
        List of 2D arrays (resolution×resolution), each a Zernike mode.
    """
    fieldcache = {}
    # build normalized coordinate grid from -1 to 1
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, -y)  # flip Y so positive up
    rho = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)
    mask = rho <= 1.0

    if n_max > 10:
        basis = generate_basis(max(0,n_max-2),resolution)
        range_n = range(n_max-2,n_max + 1)
    else:
        basis = []
        range_n = range(n_max + 1)
    for n in range_n:
        for m in range(-n, n+1):
            if (n - m) % 2 == 0:
                Z = zernike(n, m, rho, theta, fieldcache)
                Z[~mask] = 0.0
                basis.append(Z)
    return basis

def compute_moments(image, basis):
    """
    Compute Zernike moments (least‐squares fit) of an image.
    Parameters
    ----------
    image : ndarray
        2D array (resolution×resolution)
    basis : list of ndarray
        List of Zernike modes (same shape as image)
    Returns
    -------
    coeffs : ndarray
        Vector of Zernike coefficients.
    """
    N = image.size
    K = len(basis)
    A = np.zeros((N, K))
    b = image.ravel()
    for k, Z in enumerate(basis):
        A[:, k] = Z.ravel()
    coeffs, *_ = lstsq(A, b)
    return coeffs

def recreate_image(coeffs, basis):
    """
    Reconstruct an image from Zernike coefficients and basis.
    """
    img = np.zeros_like(basis[0])
    for c, Z in zip(coeffs, basis):
        img += c * Z
    return img

def circular_mask(resolution):
    """
    Return a circular mask (1 inside unit circle, 0 outside).
    """
    x = np.linspace(-1, 1, resolution)
    y = np.linspace(-1, 1, resolution)
    X, Y = np.meshgrid(x, -y)
    rho = np.sqrt(X**2 + Y**2)
    return (rho <= 1.0).astype(float)

if __name__ == "__main__":
    # Example usage
    N = 256        # grid size
    for n_max in range(10,17,2): #range(10,50,2):     # maximum radial order

        # 1) Generate Zernike basis
        s = time.time()
        basis = generate_basis(n_max, N)
        print(f"Elapsed for {n_max} = {time.time()-s}")
        # 2) Create a test image (unit‐circle mask)
        img = circular_mask(N)

        # 3) Compute Zernike moments
        coeffs = compute_moments(img, basis)
        print(f"Computed {len(coeffs)} coefficients.")

        # 4) Reconstruct image from moments
        img_rec = recreate_image(coeffs, basis)

        # 5) (Optional) Compute reconstruction error
        err = np.linalg.norm((img - img_rec).ravel()) / np.linalg.norm(img.ravel())
        print(f"@Order {n_max} Relative reconstruction error: {err:.3e}")
