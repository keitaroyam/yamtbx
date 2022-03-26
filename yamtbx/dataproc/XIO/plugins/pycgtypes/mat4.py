from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
######################################################################
# mat4 - Matrix class (4x4 matrix)
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import types, math, copy
from .vec3 import vec3 as _vec3
from .vec4 import vec4 as _vec4
from .mat3 import mat3 as _mat3


# [  0   1   2   3 ]
# [  4   5   6   7 ]
# [  8   9  10  11 ]
# [ 12  13  14  15 ]


# mat4
class mat4(object):
    """Matrix class (4x4).

    This class represents a 4x4 matrix that can be used to store
    affine transformations.
    """

    def __init__(self, *args):
        "Constructor"

        # No arguments
        if len(args)==0:
            self.mlist = 16*[0.0]

        # 1 argument (list, scalar or mat4)
        elif len(args)==1:
            T = type(args[0])
            if T==float or T==int or T==int:
                self.mlist = [args[0],0.0,0.0,0.0,
                              0.0,args[0],0.0,0.0,
                              0.0,0.0,args[0],0.0,
                              0.0,0.0,0.0,args[0]]
            # mat4
            elif isinstance(args[0], mat4):
                self.mlist = copy.copy(args[0].mlist)
            # String
            elif T==bytes:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                self.mlist=[float(x) for x in s]
            else:
                self.mlist = list(args[0])

        # 4 arguments (sequences)
        elif len(args)==4:
            a,b,c,d=args
            self.mlist = [a[0], b[0], c[0], d[0],
                          a[1], b[1], c[1], d[1],
                          a[2], b[2], c[2], d[2],
                          a[3], b[3], c[3], d[3]]

        # 16 arguments
        elif len(args)==16:
            self.mlist = list(args)

        else:
            raise TypeError("mat4() arg can't be converted to mat4")

        # Check if there are really 16 elements in the list
        if len(self.mlist)!=16:
            raise TypeError("mat4(): Wrong number of matrix elements ("+repr(len(self.mlist))+" instead of 16)")

    def __repr__(self):
        return 'mat4('+repr(self.mlist)[1:-1]+')'

    def __str__(self):
        fmt="%9.4f"
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
        return ('['+fmt%m11+', '+fmt%m12+', '+fmt%m13+', '+fmt%m14+']\n'+
                '['+fmt%m21+', '+fmt%m22+', '+fmt%m23+', '+fmt%m24+']\n'+
                '['+fmt%m31+', '+fmt%m32+', '+fmt%m33+', '+fmt%m34+']\n'+
                '['+fmt%m41+', '+fmt%m42+', '+fmt%m43+', '+fmt%m44+']')

    def __eq__(self, other):
        """== operator"""
        if isinstance(other, mat4):
            return self.mlist==other.mlist
        else:
            return 0

    def __ne__(self, other):
        """!= operator"""
        if isinstance(other, mat4):
            return self.mlist!=other.mlist
        else:
            return 1


    def __add__(self, other):
        """Matrix addition.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M+M
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        """
        if isinstance(other, mat4):
            return mat4(list(map(lambda x,y: x+y, self.mlist, other.mlist)))
        else:
            raise TypeError("unsupported operand type for +")

    def __sub__(self, other):
        """Matrix subtraction.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M-M
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        """
        if isinstance(other, mat4):
            return mat4(list(map(lambda x,y: x-y, self.mlist, other.mlist)))
        else:
            raise TypeError("unsupported operand type for -")

    def __mul__(self, other):
        """Multiplication.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M*2.0
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        >>> print 2.0*M
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        >>> print M*M
        [  90.0000,  100.0000,  110.0000,  120.0000]
        [ 202.0000,  228.0000,  254.0000,  280.0000]
        [ 314.0000,  356.0000,  398.0000,  440.0000]
        [ 426.0000,  484.0000,  542.0000,  600.0000]
        >>> print M*_vec3(1,2,3)
        (0.1765, 0.4510, 0.7255)
        >>> print _vec3(1,2,3)*M
        (0.7083, 0.8056, 0.9028)
        """
        T = type(other)
        # mat4*scalar
        if T==float or T==int or T==int:
            return mat4(list(map(lambda x,other=other: x*other, self.mlist)))
        # mat4*vec3
        if isinstance(other, _vec3):
            m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
            w = float(m41*other.x + m42*other.y + m43*other.z + m44)
            return _vec3(m11*other.x + m12*other.y + m13*other.z + m14, 
                         m21*other.x + m22*other.y + m23*other.z + m24, 
                         m31*other.x + m32*other.y + m33*other.z + m34)/w
        # mat4*vec4
        if isinstance(other, _vec4):
            m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
            return _vec4(m11*other.x + m12*other.y + m13*other.z + m14*other.w, 
                         m21*other.x + m22*other.y + m23*other.z + m24*other.w, 
                         m31*other.x + m32*other.y + m33*other.z + m34*other.w,
                         m41*other.x + m42*other.y + m43*other.z + m44*other.w)
        # mat4*mat4
        if isinstance(other, mat4):
            m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
            n11,n12,n13,n14,n21,n22,n23,n24,n31,n32,n33,n34,n41,n42,n43,n44 = other.mlist
            return mat4( m11*n11+m12*n21+m13*n31+m14*n41,
                         m11*n12+m12*n22+m13*n32+m14*n42,
                         m11*n13+m12*n23+m13*n33+m14*n43,
                         m11*n14+m12*n24+m13*n34+m14*n44,

                         m21*n11+m22*n21+m23*n31+m24*n41,
                         m21*n12+m22*n22+m23*n32+m24*n42,
                         m21*n13+m22*n23+m23*n33+m24*n43,
                         m21*n14+m22*n24+m23*n34+m24*n44,

                         m31*n11+m32*n21+m33*n31+m34*n41,
                         m31*n12+m32*n22+m33*n32+m34*n42,
                         m31*n13+m32*n23+m33*n33+m34*n43,
                         m31*n14+m32*n24+m33*n34+m34*n44,

                         m41*n11+m42*n21+m43*n31+m44*n41,
                         m41*n12+m42*n22+m43*n32+m44*n42,
                         m41*n13+m42*n23+m43*n33+m44*n43,
                         m41*n14+m42*n24+m43*n34+m44*n44)
        # unsupported
        else:
            raise TypeError("unsupported operand type for *")

    def __rmul__(self, other):
        T = type(other)
        # scalar*mat4
        if T==float or T==int or T==int:
            return mat4(list(map(lambda x,other=other: other*x, self.mlist)))
        # vec4*mat4
        if isinstance(other, _vec4):
            m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
            return _vec4(other.x*m11 + other.y*m21 + other.z*m31 + other.w*m41, 
                         other.x*m12 + other.y*m22 + other.z*m32 + other.w*m42,
                         other.x*m13 + other.y*m23 + other.z*m33 + other.w*m43,
                         other.x*m14 + other.y*m24 + other.z*m34 + other.w*m44)
        # vec3*mat4
        if isinstance(other, _vec3):
            m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
            w = float(other.x*m14 + other.y*m24 + other.z*m34 + m44)
            return _vec3(other.x*m11 + other.y*m21 + other.z*m31 + m41, 
                         other.x*m12 + other.y*m22 + other.z*m32 + m42,
                         other.x*m13 + other.y*m23 + other.z*m33 + m43)/w
        # mat4*mat4
        if isinstance(other, mat4):
            return self.__mul__(other)
        # unsupported
        else:
            raise TypeError("unsupported operand type for *")

    def __div__(self, other):
        """Division
        
        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M/2.0
        [   0.5000,    1.0000,    1.5000,    2.0000]
        [   2.5000,    3.0000,    3.5000,    4.0000]
        [   4.5000,    5.0000,    5.5000,    6.0000]
        [   6.5000,    7.0000,    7.5000,    8.0000]
        """
        T = type(other)
        # mat4/scalar
        if T==float or T==int or T==int:
            return mat4(map(lambda x,other=other: x/other, self.mlist))
        # unsupported
        else:
            raise TypeError("unsupported operand type for /")

    def __mod__(self, other):
        """Modulo.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M%5.0
        [   1.0000,    2.0000,    3.0000,    4.0000]
        [   0.0000,    1.0000,    2.0000,    3.0000]
        [   4.0000,    0.0000,    1.0000,    2.0000]
        [   3.0000,    4.0000,    0.0000,    1.0000]
        """
        T = type(other)
        # mat4%scalar
        if T==float or T==int or T==int:
            return mat4(list(map(lambda x,other=other: x%other, self.mlist)))
        # unsupported
        else:
            raise TypeError("unsupported operand type for %")

    def __neg__(self):
        """Negation.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print -M
        [  -1.0000,   -2.0000,   -3.0000,   -4.0000]
        [  -5.0000,   -6.0000,   -7.0000,   -8.0000]
        [  -9.0000,  -10.0000,  -11.0000,  -12.0000]
        [ -13.0000,  -14.0000,  -15.0000,  -16.0000]
        """
        return mat4([-x for x in self.mlist])

    def __pos__(self):
        """
        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print +M
        [   1.0000,    2.0000,    3.0000,    4.0000]
        [   5.0000,    6.0000,    7.0000,    8.0000]
        [   9.0000,   10.0000,   11.0000,   12.0000]
        [  13.0000,   14.0000,   15.0000,   16.0000]
        """
        return mat4([+x for x in self.mlist])


    def __len__(self):
        return 4

    def __getitem__(self, key):
        if type(key)==int:
            if key<0 or key>3:
                raise IndexError("index out of range")
            m=self.mlist
            if   key==0: return [m[0],m[4],m[8],m[12]]
            elif key==1: return [m[1],m[5],m[9],m[13]]
            elif key==2: return [m[2],m[6],m[10],m[14]]
            elif key==3: return [m[3],m[7],m[11],m[15]]
        elif type(key)==tuple:
            i,j=key
            if i<0 or i>3 or j<0 or j>3:
                raise IndexError("index out of range")
            return self.mlist[i*4+j]
        else:
            raise TypeError("index must be integer or 2-tuple")

    def __setitem__(self, key, value):
        if type(key)==int:
            if key<0 or key>3:
                raise IndexError("index out of range")
            m=self.mlist
            if   key==0: m[0],m[4],m[8],m[12]=value
            elif key==1: m[1],m[5],m[9],m[13]=value
            elif key==2: m[2],m[6],m[10],m[14]=value
            elif key==3: m[3],m[7],m[11],m[15]=value
        elif type(key)==tuple:
            i,j=key
            if i<0 or i>3 or j<0 or j>3:
                raise IndexError("index out of range")
            self.mlist[i*4+j] = value
        else:
            raise TypeError("index must be integer or 2-tuple")

    def getRow(self, idx):
        """Return row (as vec4)."""
        m=self.mlist
        if   idx==0: return _vec4(m[0], m[1], m[2], m[3])
        elif idx==1: return _vec4(m[4], m[5], m[6], m[7])
        elif idx==2: return _vec4(m[8], m[9], m[10], m[11])
        elif idx==3: return _vec4(m[12], m[13], m[14], m[15])
        else:
            raise IndexError("index out of range")

    def setRow(self, idx, value):
        """Set row."""
        m=self.mlist
        if   idx==0: m[0],m[1],m[2],m[3] = value
        elif idx==1: m[4],m[5],m[6],m[7] = value
        elif idx==2: m[8],m[9],m[10],m[11] = value
        elif idx==3: m[12],m[13],m[14],m[15] = value
        else:
            raise IndexError("index out of range")

    def getColumn(self, idx):
        """Return column (as vec4)."""
        m=self.mlist
        if   idx==0: return _vec4(m[0], m[4], m[8], m[12])
        elif idx==1: return _vec4(m[1], m[5], m[9], m[13])
        elif idx==2: return _vec4(m[2], m[6], m[10], m[14])
        elif idx==3: return _vec4(m[3], m[7], m[11], m[15])
        else:
            raise IndexError("index out of range")

    def setColumn(self, idx, value):
        """Set column."""
        m=self.mlist
        if   idx==0: m[0],m[4],m[8],m[12] = value
        elif idx==1: m[1],m[5],m[9],m[13] = value
        elif idx==2: m[2],m[6],m[10],m[14] = value
        elif idx==3: m[3],m[7],m[11],m[15] = value
        else:
            raise IndexError("index out of range")

    def toList(self, rowmajor=0):
        """Return a list containing the matrix elements.

        By default the list is in column-major order (which can directly be
        used in OpenGL or RenderMan). If you set the optional argument
        rowmajor to 1, you'll get the list in row-major order.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M.toList()
        [1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]
        >>> print M.toList(rowmajor=1)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        """
        if rowmajor:
            return copy.copy(self.mlist)
        else:
            return self.transpose().mlist
            

    def identity(self):
        """Return identity matrix.

        >>> print mat4().identity()
        [   1.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    1.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    1.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    1.0000]
        """
        return mat4(1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0)

    def transpose(self):
        """Transpose matrix.
        
        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M.transpose()
        [   1.0000,    5.0000,    9.0000,   13.0000]
        [   2.0000,    6.0000,   10.0000,   14.0000]
        [   3.0000,    7.0000,   11.0000,   15.0000]
        [   4.0000,    8.0000,   12.0000,   16.0000]
        """
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
        return mat4(m11,m21,m31,m41,
                    m12,m22,m32,m42,
                    m13,m23,m33,m43,
                    m14,m24,m34,m44)

    def determinant(self):
        """Return determinant.
        
        >>> M=mat4(2.0,0,0,0, 0,2.0,0,0, 0,0,2.0,0, 0,0,0,2.0)
        >>> print M.determinant()
        16.0
        """
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist

        return m11*m22*m33*m44 \
               -m11*m22*m34*m43 \
               +m11*m23*m34*m42 \
               -m11*m23*m32*m44 \
               +m11*m24*m32*m43 \
               -m11*m24*m33*m42 \
               -m12*m23*m34*m41 \
               +m12*m23*m31*m44 \
               -m12*m24*m31*m43 \
               +m12*m24*m33*m41 \
               -m12*m21*m33*m44 \
               +m12*m21*m34*m43 \
               +m13*m24*m31*m42 \
               -m13*m24*m32*m41 \
               +m13*m21*m32*m44 \
               -m13*m21*m34*m42 \
               +m13*m22*m34*m41 \
               -m13*m22*m31*m44 \
               -m14*m21*m32*m43 \
               +m14*m21*m33*m42 \
               -m14*m22*m33*m41 \
               +m14*m22*m31*m43 \
               -m14*m23*m31*m42 \
               +m14*m23*m32*m41
    

    def _submat(self, i,j):
        M=_mat3()
        for k in range(3):
            for l in range(3):
                t=(k,l)
                if k>=i:
                    t=(k+1,t[1])
                if l>=j:
                    t=(t[0],l+1)
                M[k,l] = self[t]
        return M
        
    def inverse(self):
        """Return inverse matrix.

        >>> M=mat4(0,-2.0,0,0, 2.0,0,0,0, 0,0,2,0, 0,0,0,2)
        >>> print M.inverse()
        [   0.0000,    0.5000,    0.0000,    0.0000]
        [  -0.5000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.5000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.5000]
        """
        
        Mi=mat4()
        d=self.determinant()
        for i in range(4):
            for j in range(4):
                sign=1-((i+j)%2)*2
                m3=self._submat(i,j)
                Mi[j,i]=sign*m3.determinant()/d
        return Mi

    def translation(self, t):
        """Return translation matrix."""
        return mat4(1.0, 0.0, 0.0, t.x,
                    0.0, 1.0, 0.0, t.y,
                    0.0, 0.0, 1.0, t.z,
                    0.0, 0.0, 0.0, 1.0)

    def scaling(self, s):
        """Return scaling matrix."""
        return mat4(s.x, 0.0, 0.0, 0.0,
                    0.0, s.y, 0.0, 0.0,
                    0.0, 0.0, s.z, 0.0,
                    0.0, 0.0, 0.0, 1.0)

    def rotation(self, angle, axis):
        """Return rotation matrix.

        angle must be given in radians. axis should be of type vec3.
        """

        sqr_a = axis.x*axis.x
        sqr_b = axis.y*axis.y
        sqr_c = axis.z*axis.z
        len2  = sqr_a+sqr_b+sqr_c

        k2    = math.cos(angle)
        k1    = (1.0-k2)/len2
        k3    = math.sin(angle)/math.sqrt(len2)
        k1ab  = k1*axis.x*axis.y
        k1ac  = k1*axis.x*axis.z
        k1bc  = k1*axis.y*axis.z
        k3a   = k3*axis.x
        k3b   = k3*axis.y
        k3c   = k3*axis.z

        return mat4( k1*sqr_a+k2, k1ab-k3c, k1ac+k3b, 0.0,
                     k1ab+k3c, k1*sqr_b+k2, k1bc-k3a, 0.0,
                     k1ac-k3b, k1bc+k3a, k1*sqr_c+k2, 0.0,
                     0.0, 0.0, 0.0, 1.0)

    def translate(self, t):
        """Concatenate a translation."""
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
        self.mlist[3]  = m11*t.x + m12*t.y + m13*t.z + m14
        self.mlist[7]  = m21*t.x + m22*t.y + m23*t.z + m24
        self.mlist[11] = m31*t.x + m32*t.y + m33*t.z + m34
        self.mlist[15] = m41*t.x + m42*t.y + m43*t.z + m44
        return self

    def scale(self, s):
        """Concatenate a scaling."""
        self.mlist[0]  *= s.x
        self.mlist[1]  *= s.y
        self.mlist[2]  *= s.z
        self.mlist[4]  *= s.x
        self.mlist[5]  *= s.y
        self.mlist[6]  *= s.z
        self.mlist[8]  *= s.x
        self.mlist[9]  *= s.y
        self.mlist[10] *= s.z
        self.mlist[12] *= s.x
        self.mlist[13] *= s.y
        self.mlist[14] *= s.z
        return self

    def rotate(self, angle, axis):
        """Concatenate a rotation.

        angle must be given in radians. axis should be of type vec3.
        """
        R=self.rotation(angle, axis)
        self.mlist = (self*R).mlist
        return self


    def frustum(self, left, right, bottom, top, near, far):
        """equivalent to the OpenGL command glFrustum()"""
        
        return mat4( (2.0*near)/(right-left), 0.0, float(right+left)/(right-left), 0.0,
                     0.0, (2.0*near)/(top-bottom), float(top+bottom)/(top-bottom), 0.0,
                     0.0, 0.0, -float(far+near)/(far-near), -(2.0*far*near)/(far-near),
                     0.0, 0.0, -1.0, 0.0)
    
    def perspective(self, fovy, aspect, near, far):
        """von Mesa ubernommen (glu.c)"""

        top    = near * math.tan(fovy * math.pi / 360.0)
        bottom = -top
        left   = bottom * aspect
        right  = top * aspect

        return self.frustum(left, right, bottom, top, near, far)

    def lookAt(self, pos, target, up=_vec3(0,0,1)):
        """Look from pos to target.

        The resulting transformation moves the origin to pos and
        rotates so that The z-axis points to target. The y-axis is
        as close as possible to the up vector.
        """
        dir = (target - pos).normalize()
        up  = up.normalize()
        up -= (up * dir) * dir
        try:
            up  = up.normalize()
        except:
            # We're looking along the up direction, so choose
            # an arbitrary direction that is perpendicular to dir
            # as new up.
            up = dir.ortho()

        right = up.cross(dir).normalize()

        self.mlist=[right.x, up.x, dir.x, pos.x,
                    right.y, up.y, dir.y, pos.y,
                    right.z, up.z, dir.z, pos.z,
                    0.0, 0.0, 0.0, 1.0]
        return self

    def ortho(self):
        """Return a matrix with orthogonal base vectors.

        Makes the x-, y- and z-axis orthogonal.
        The fourth column and row remain untouched.
        """

        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist

        x = _vec3(m11, m21, m31)
        y = _vec3(m12, m22, m32)
        z = _vec3(m13, m23, m33)

        xl = x.length()
        xl*=xl
        y = y - ((x*y)/xl)*x
        z = z - ((x*z)/xl)*x

        yl = y.length()
        yl*=yl
        z = z - ((y*z)/yl)*y

        return mat4( x.x, y.x, z.x, m14,
                     x.y, y.y, z.y, m24,
                     x.z, y.z, z.z, m34,
                     m41, m42, m43, m44)

    def decompose(self):
        """Decomposes the matrix into a translation, rotation and scaling part.

        Returns a tuple (translation, rotation, scaling). The 
        translation and scaling parts are given as vec3's, the rotation
        is still given as a mat4.
        """
        dummy = self.ortho()
        dummy.setRow(3,_vec4(0.0, 0.0, 0.0, 1.0))

        x = dummy.getColumn(0)
        y = dummy.getColumn(1)
        z = dummy.getColumn(2)
        xl = x.length()
        yl = y.length()
        zl = z.length()
        scale = _vec3(xl,yl,zl)
        
        x/=xl
        y/=yl
        z/=zl
        dummy.setColumn(0,x)
        dummy.setColumn(1,y)
        dummy.setColumn(2,z)
        if dummy.determinant()<0.0:
            dummy.setColumn(0,-x)
            scale.x=-scale.x

        return (_vec3(self.mlist[3], self.mlist[7], self.mlist[11]),
                dummy,
                scale)

    def getMat3(self):
        """Convert to mat3 by discarding 4th row and column.
        """
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = self.mlist
        return _mat3(m11,m12,m13,
                     m21,m22,m23,
                     m31,m32,m33)
        
######################################################################

def _test():
    import doctest, mat4
    failed, total = doctest.testmod(mat4)
    print("%d/%d failed" % (failed, total))

if __name__=="__main__":

    _test()


