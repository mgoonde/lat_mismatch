# from: https://www.embeddedrelated.com/showarticle/1620.php
import math

def cfrac_approx(x,kmax=None,nmax=None,include_next=False,verbose=False):
    '''
    x =approx= h/k

    input:
    kmax : max value of k
    nmax : max number of iterations
    include_next : ?

    returns:
    coeffs : integers of continued fraction (dont really care)
    brow[1] : h
    brow[2] : k
    '''
    if kmax is None and nmax is None:
        raise ValueError('Need to specify either max denominator kmax or max number of coefficients nmax')
    if kmax is None:
        kmax = float('inf')
    if nmax is None:
        nmax = float('inf')
    nmin=3
    arow = [1,0,1]
    brow = [x,1,0]
    coeffs = []
    n = 0
    while n < nmax:
        q = math.floor(x)
        r = x - q
        if r == 0:
            break
        n += 1
        x = 1/r
        xrow = [x, arow[1]+q*brow[1], arow[2]+q*brow[2]]
        if n>nmin:
           #if (xrow[2] > kmax and not include_next) or brow[2] > kmax:
           if (xrow[2] > kmax and not include_next):
              break
        arow = brow
        brow = xrow
        coeffs.append(q)
        if verbose:
            print(brow, q)
    if len(coeffs) > 1 and coeffs[-1] == 1:
        coeffs = coeffs[:-1]
        coeffs[-1] += 1
    return coeffs, brow[1], brow[2]


def frac_h(x,kmax):
    '''
    return h
    '''
    cc, h, k = cfrac_approx( x, kmax=kmax, nmax=None )
    return h

def frac_k(x,kmax):
    '''
    return k
    '''
    cc, h, k = cfrac_approx( x, kmax=kmax, nmax=None )
    return k
