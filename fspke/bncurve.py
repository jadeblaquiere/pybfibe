# Copyright (c) 2017, Joseph deBlaquiere <jadeblaquiere@yahoo.com>
# All rights reserved
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of ciphrtxt nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from Crypto.Random import random
from fspke.rabinmiller import isPrime
from ecpy.point import Point, Generator

import math

def _pfunc(x):
    return (36 * (x **4)) + (36 * (x ** 3)) + (24 * (x ** 2)) + (6 * x) + 1

def _nfunc(x):
    return (36 * (x **4)) + (36 * (x ** 3)) + (18 * (x ** 2)) + (6 * x) + 1

def _trfunc(x):
    return (6 * (x ** 2)) + 1

def _legendre(a, p):
    x = pow(a, (p-1)//2, p)
    if x == 0 or x == 1:
        return x
    if x == p-1:
        return -1
    assert False

def _mod_sqrt(a, p):
    assert p % 4 == 3
    return pow(a, (p+1)//4, p)

def _modinv(a, m):
    lastr, r, x, lastx = a, m, 0, 1
    while r:
        lastr, (q, r) = r, divmod(lastr, r)
        x, lastx = lastx - q*x, x
    return lastx % m


# foundation for CHK : https://eprint.iacr.org/2003/083.pdf

class BNCurve(object):
    """Implementation of Barreto-Naehrig pairing-friendly curves 
       BN : https://eprint.iacr.org/2010/429.pdf
       BN implementation : https://eprint.iacr.org/2010/354.pdf
       nbit = bitsize target
       u = parameter setting field and order
       b = curve constant (curve is y**2 = x**3 + ax + b, a=0 in BN)
       G = generator point (as tuple) (must also provide B, u)
    """
    def __init__(self, nbit, b=None, u=None, G=None, friendly=True):
        # first, find umax, umin such that p < 2**nbit, p > 2**(nbit - 1)
        # p = 36u**4 + 36u**3 + 24u**2 + 6u + 1
        # n = 36u**4 + 36u**3 + 18u**2 + 6u + 1
        limit = pow(2,nbit)
        if u is None:
            # derive u randomly based on nbit
            # coarse upper, lower guesses
            upr = int((limit/36) ** 0.25)
            lwr = 0
            mid = (upr + lwr) >> 1
            # midpoint algorithm
            while ((mid != upr) or (mid != lwr)):
                mid = (upr + lwr) >> 1
                pu = _pfunc(upr)
                pl = _pfunc(lwr)
                pm = _pfunc(mid)
                # print("limit goal = %d" % (limit))
                # print("lower p(%d) = %d" % (lwr, pl))
                # print("mid p(%d) = %d" % (mid, pm))
                # print("upper p(%d) = %d" % (upr, pu))
                # print("")
                # pm should be odd, limit even. They should never meet
                assert (pm != limit)
                if pm < limit:
                    if lwr == mid:
                        break
                    lwr = mid
                else:
                    upr = mid
            umax = mid
            # repeat for limit = 2 ** (nbit-1)
            limit = pow(2,(nbit - 1))
            upr = int((limit/36) ** 0.25)
            lwr = 0
            mid = (upr + lwr) >> 1
            while ((mid != upr) or (mid != lwr)):
                mid = (upr + lwr) >> 1
                pu = _pfunc(upr)
                pl = _pfunc(lwr)
                pm = _pfunc(mid)
                #pm should be odd, limit should be even. They should never meet
                assert (pm != limit)
                if pm < limit:
                    if lwr == mid:
                        break
                    lwr = mid
                else:
                    upr = mid
            umin = mid
        else:
            # u is a parameter
            umax = u
            umin = umax
        print("limit goal = %X" % (pow(2,nbit)))
        print("umax = %X, pmax = %X" % (umax, _pfunc(umax)))
        print("umin = %X, pmin = %X" % (umin, _pfunc(umin)))
        print("(%d potential u values, approx 2 ** %d)" % ((umax - umin)+1, int(math.log((umax - umin)+1, 2)+0.5)))
        # choose u at random until valid solution, follow algorithm 1 from:
        # https://www.cryptojedi.org/papers/pfcpo.pdf
        self.u = 0
        if u is None:
            while True:
                urand = random.randint(umin, umax)
                #urand = 6518589491078791937
                p = _pfunc(urand)
                #assert p == 65000549695646603732796438742359905742825358107623003571877145026864184071783
                if isPrime(p) != True:
                    continue
                n = _nfunc(urand)
                #assert n == 65000549695646603732796438742359905742570406053903786389881062969044166799969
                if isPrime(n) != True:
                    continue
                if (p % 4) == 3:
                    if ((p * p) % 16) == 9:
                        self.u = urand
                        self.p = p
                        self.n = n
                        break
        else:
            p = _pfunc(umax)
            n = _nfunc(umax)
            self.u = umax
            self.p = p
            self.n = n
        if b is not None:
            friendly = False
        # print("u = %X (%d)" % (self.u, self.u))
        # print("p = %X (%d)" % (self.p, self.p))
        # print("n = %X (%d)" % (self.n, self.n))
        if friendly:
            # print("searching for b")
            while True:
                #select friendly parameters per Pereira et al
                #we happen to choose d as pure imaginary such that Xi is not
                c = random.randint(0,self.p)
                di = random.randint(0,self.p)
                b = pow(c, 4, self.p) - pow(di, 6, self.p)
                if b != (b % self.p):
                    continue
                if _legendre(b + 1, self.p) != 1:
                    continue
                y = _mod_sqrt(b + 1, self.p)
                # print("trying b = %X, y = %X, y**2 = %X" % (b, y, ((y * y) % p)))
                p2 = self.p * self.p
                # 
                curve1 = { 'p' : self.p, 
                          'a' : 0,
                          'b' : b,
                          'n' : self.n,
                          'bits' : nbit }
                Generator.set_curve(curve1)
                P = Point(pow(di, 2, self.p), pow(c, 2, self.p))
                assert P.is_valid()
                Pnm1 = (self.n - 1) * P
                Ptst = Pnm1 + P
                if Ptst.is_infinite != True:
                    continue
                self.G = Generator(pow(di, 2, self.p), pow(c, 2, self.p))
                self.b = b
                Xi = pow(c, 2, p2) + pow(di, 3, p2)
                bprime = (b * _modinv(Xi, p2))
                h = (2 * self.p) - self.n
                curve2 = { 'p' : self.p * self.p,
                           'a' : 0,
                           'b' : bprime,
                           'n' : self.n,
                           'bits' : 2 * nbit}
                break
        else:
            if G is None:
                if b is None:
                    b = 0
                else:
                    b = b - 1
                while True:
                    b = b + 1
                    if _legendre(b + 1, self.p) != 1:
                        # print("b = %d but %d is not quadratic residue" % (b, b+1))
                        continue
                    y = _mod_sqrt(b + 1, self.p)
                    # print("trying b = %X, y = %X, y**2 = %X" % (b, y, ((y * y) % p)))
                    curve = { 'p' : self.p, 
                              'a' : 0,
                              'b' : b,
                              'n' : self.n,
                              'bits' : nbit }
                    Generator.set_curve(curve)
                    P = Point(1, y)
                    assert P.is_valid()
                    Pnm1 = (self.n - 1) * P
                    Ptst = Pnm1 + P
                    Pa = P.affine()
                    if Ptst.is_infinite == True:
                        self.G = Generator(Pa[0], Pa[1])
                        self.b = b
                        break
            else:
                assert b is not None
                curve = { 'p' : self.p, 
                          'a' : 0,
                          'b' : b,
                          'n' : self.n,
                          'bits' : nbit }
                Generator.set_curve(curve)
                P = Point(G[0], G[1])
                # print("P = %s" % (P.uncompressed_format()))
                assert P.is_valid()
                Pnm1 = (self.n - 1) * P
                Ptst = Pnm1 + P
                Pa = P.affine()
                assert Ptst.is_infinite == True
                self.G = Generator(Pa[0], Pa[1])
                self.b = b
        #self.p = _pfunc(self.u)
        #self.n = _nfunc(self.u)
        print("u = %X (%d)" % (self.u, self.u))
        print("p = %X (%d)" % (self.p, self.p))
        print("n = %X (%d)" % (self.n, self.n))
        print("b = %X (%d)" % (self.b, self.b))
        print("P (compressed) = %s" % (self.G.compress()))
        print("P (uncompressed) = %s" % (self.G.uncompressed_format()))
        
            
