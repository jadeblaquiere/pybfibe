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

from fspke.bncurve import *

# Generate Friendly 256-bit BN curve, random u (and p,n) 
bn1 = BNCurve(256)
# Generate Friendly 256-bit BN curve, fixed u (and p,n)
bn1 = BNCurve(256, u=6518589491078791937)
# BN2,254 Curve benchmarked by Pereira
bn2_254 = BNCurve(254,2,-(pow(2,62) + pow(2,55) + 1))
bn2_254 = BNCurve(254,2,-(pow(2,62) + pow(2,55) + 1), (-1, 1))
# ISO BN2,256 curve
bn2_256 = BNCurve(256, 3, 6518589491078791937, (1, -2))
assert bn2_256.u == 6518589491078791937
assert bn2_256.p == 65000549695646603732796438742359905742825358107623003571877145026864184071783
assert bn2_256.n == 65000549695646603732796438742359905742570406053903786389881062969044166799969
assert bn2_256.G.affine()[0] == (1 % bn2_256.p)
assert bn2_256.G.affine()[1] == (-2 % bn2_256.p)
bn512 = BNCurve(512)
