from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
######################################################################
# mat3 - Matrix class (3x3 matrix)
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import types, math, copy
from .vec3 import vec3 as _vec3

# [  0   1   2 ]
# [  3   4   5 ]
# [  6   7   8 ]

# mat3
class mat3(object):
    """Matrix class (3x3).

    This class represents a 3x3 matrix that can be used to store
    linear transformations.
    """

    def __init__(self, *args):
        """Constructor.

        There are several possibilities how to initialize a matrix,
        depending on the number of arguments you provide to the constructor.

        - 0 arguments: Every component is zero.
        - 1 number argument: The diagonal is initialized to that number,
          all the other elements are zero.
        - 1 sequence argument: The elements are initialized with the numbers
          in the sequence (the sequence must contain 9 numbers).
        - 1 mat3 argument: The matrix is copied.
        - 3 sequence arguments: The columns are initialized with the
          respective sequence (each sequence must contain 3 numbers).
        - 9 number arguments: The matrix is initialized with those values
          (row-major order).
        """

        # No arguments
        if len(args)==0:
            self.mlist = 9*[0.0]

        # 1 argument (list, scalar or mat3)
        elif len(args)==1:
            T = type(args[0])
            # Scalar
            if T==float or T==int or T==int:
                self.mlist = [args[0],0.0,0.0,
                              0.0,args[0],0.0,
                              0.0,0.0,args[0]]
            # mat3
            elif isinstance(args[0], mat3):
                self.mlist = copy.copy(args[0].mlist)
            # String
            elif T==bytes:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                self.mlist=[float(x) for x in s]
            else:
                self.mlist = list(args[0])

        # 3 arguments (sequences)
        elif len(args)==3:
            a,b,c=args
            self.mlist = [a[0], b[0], c[0],
                          a[1], b[1], c[1],
                          a[2], b[2], c[2]]

        # 9 arguments
        elif len(args)==9:
            self.mlist = list(args)

        else:
            raise TypeError("mat3() arg can't be converted to mat3")

        # Check if there are really 9 elements in the list
        if len(self.mlist)!=9:
            raise TypeError("mat4(): Wrong number of matrix elements ("+repr(len(self.mlist))+" instead of 9)")


    def __repr__(self):
        return 'mat3('+repr(self.mlist)[1:-1]+')'

    def __str__(self):
        fmt="%9.4f"
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
        return ('['+fmt%m11+', '+fmt%m12+', '+fmt%m13+']\n'+
                '['+fmt%m21+', '+fmt%m22+', '+fmt%m23+']\n'+
                '['+fmt%m31+', '+fmt%m32+', '+fmt%m33+']')

    def __eq__(self, other):
        """== operator"""
        if isinstance(other, mat3):
            return self.mlist==other.mlist
        else:
            return 0

    def __ne__(self, other):
        """!= operator"""
        if isinstance(other, mat3):
            return self.mlist!=other.mlist
        else:
            return 1


    def __add__(self, other):
        if isinstance(other, mat3):
            return mat3(list(map(lambda x,y: x+y, self.mlist, other.mlist)))
        else:
            raise TypeError("unsupported operand type for +")

    def __sub__(self, other):
        if isinstance(other, mat3):
            return mat3(list(map(lambda x,y: x-y, self.mlist, other.mlist)))
        else:
            raise TypeError("unsupported operand type for -")

    def __mul__(self, other):
        T = type(other)
        # mat3*scalar
        if T==float or T==int or T==int:
            return mat3(list(map(lambda x,other=other: x*other, self.mlist)))
        # mat3*vec3
        if isinstance(other, _vec3):
            m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
            return _vec3(m11*other.x + m12*other.y + m13*other.z, 
                         m21*other.x + m22*other.y + m23*other.z, 
                         m31*other.x + m32*other.y + m33*other.z)            
        # mat3*mat3
        if isinstance(other, mat3):
            m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
            n11,n12,n13,n21,n22,n23,n31,n32,n33 = other.mlist
            return mat3( m11*n11+m12*n21+m13*n31,
                         m11*n12+m12*n22+m13*n32,
                         m11*n13+m12*n23+m13*n33,

                         m21*n11+m22*n21+m23*n31,
                         m21*n12+m22*n22+m23*n32,
                         m21*n13+m22*n23+m23*n33,

                         m31*n11+m32*n21+m33*n31,
                         m31*n12+m32*n22+m33*n32,
                         m31*n13+m32*n23+m33*n33)
        # unsupported
        else:
            raise TypeError("unsupported operand type for *")

    def __rmul__(self, other):
        T = type(other)
        # scalar*mat3
        if T==float or T==int or T==int:
            return mat3(list(map(lambda x,other=other: other*x, self.mlist)))
        # vec3*mat3
        if isinstance(other, _vec3):
            m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
            return _vec3(other.x*m11 + other.y*m21 + other.z*m31, 
                         other.x*m12 + other.y*m22 + other.z*m32, 
                         other.x*m13 + other.y*m23 + other.z*m33)
        # mat3*mat3
        if isinstance(other, mat3):
            return self.__mul__(other)
        # unsupported
        else:
            raise TypeError("unsupported operand type for *")


    def __div__(self, other):
        T = type(other)
        # mat3/scalar
        if T==float or T==int or T==int:
            return mat3(list(map(lambda x,other=other: x/other, self.mlist)))
        # unsupported
        else:
            raise TypeError("unsupported operand type for /")

    def __mod__(self, other):
        T = type(other)
        # mat3%scalar
        if T==float or T==int or T==int:
            return mat3(list(map(lambda x,other=other: x%other, self.mlist)))
        # unsupported
        else:
            raise TypeError("unsupported operand type for %")


    def __neg__(self):
        return mat3([-x for x in self.mlist])

    def __pos__(self):
        return mat3([+x for x in self.mlist])


    def __len__(self):
        return 3

    def __getitem__(self, key):
        """Return a column or an individual element."""
        if   key==0: return _vec3(self.mlist[0],self.mlist[3],self.mlist[6])
        elif key==1: return _vec3(self.mlist[1],self.mlist[4],self.mlist[7])
        elif key==2: return _vec3(self.mlist[2],self.mlist[5],self.mlist[8])
        elif type(key)==tuple:
            i,j=key
            if i<0 or i>2 or j<0 or j>2:
                raise IndexError("index out of range")
            return self.mlist[i*3+j]
        else:
            raise IndexError("index out of range")

    def __setitem__(self, key, value):
        """Set a column or an individual element."""
        if   key==0: self.mlist[0],self.mlist[3],self.mlist[6] = value
        elif key==1: self.mlist[1],self.mlist[4],self.mlist[7] = value
        elif key==2: self.mlist[2],self.mlist[5],self.mlist[8] = value
        elif type(key)==tuple:
            i,j=key
            if i<0 or i>2 or j<0 or j>2:
                raise IndexError("index out of range")
            self.mlist[i*3+j] = value
        else:
            raise TypeError("index must be integer or 2-tuple")

    def getRow(self, idx):
        """Return a row (as vec3)."""
        if   idx==0: return _vec3(self.mlist[0], self.mlist[1], self.mlist[2])
        elif idx==1: return _vec3(self.mlist[3], self.mlist[4], self.mlist[5])
        elif idx==2: return _vec3(self.mlist[6], self.mlist[7], self.mlist[8])
        else:
            raise IndexError("index out of range")

    def setRow(self, idx, value):
        """Set a row (as vec3)."""
        if   idx==0: self.mlist[0], self.mlist[1], self.mlist[2] = value
        elif idx==1: self.mlist[3], self.mlist[4], self.mlist[5] = value
        elif idx==2: self.mlist[6], self.mlist[7], self.mlist[8] = value
        else:
            raise IndexError("index out of range")

    def getColumn(self, idx):
        """Return a column (as vec3)."""
        if   idx==0: return _vec3(self.mlist[0], self.mlist[3], self.mlist[6])
        elif idx==1: return _vec3(self.mlist[1], self.mlist[4], self.mlist[7])
        elif idx==2: return _vec3(self.mlist[2], self.mlist[5], self.mlist[8])
        else:
            raise IndexError("index out of range")

    def setColumn(self, idx, value):
        """Set a column."""
        if   idx==0: self.mlist[0], self.mlist[3], self.mlist[6] = value
        elif idx==1: self.mlist[1], self.mlist[4], self.mlist[7] = value
        elif idx==2: self.mlist[2], self.mlist[5], self.mlist[8] = value
        else:
            raise IndexError("index out of range")

    def toList(self, rowmajor=0):
        """Create a list containing the matrix elements.

        By default the list is in column-major order. If you set the
        optional argument rowmajor to 1, you'll get the list in row-major
        order.
        """
        if rowmajor:
            return copy.copy(self.mlist)
        else:
            return self.transpose().mlist

    def identity(self):
        """Return the identity matrix."""
        return mat3(1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0)

    def transpose(self):
        """Return the transposed matrix."""
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
        return mat3(m11,m21,m31,
                    m12,m22,m32,
                    m13,m23,m33)

    def determinant(self):
        """Return determinant."""
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
        return m11*m22*m33+ \
               m12*m23*m31+ \
               m13*m21*m32- \
               m31*m22*m13- \
               m32*m23*m11- \
               m33*m21*m12

    def inverse(self):
        """Return inverse matrix."""
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
        d = 1.0/self.determinant()
        return mat3( m22*m33-m23*m32, m32*m13-m12*m33, m12*m23-m22*m13,
                     m23*m31-m21*m33, m11*m33-m31*m13, m21*m13-m11*m23,
                     m21*m32-m31*m22, m31*m12-m11*m32, m11*m22-m12*m21 )*d

    def scaling(self, s):
        """Return a scale transformation."""
        return mat3(s.x, 0.0, 0.0,
                    0.0, s.y, 0.0,
                    0.0, 0.0, s.z)

    def rotation(self, angle, axis):
        """Return a rotation matrix."""
        
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

        return mat3( k1*sqr_a+k2, k1ab-k3c, k1ac+k3b,
                     k1ab+k3c, k1*sqr_b+k2, k1bc-k3a,
                     k1ac-k3b, k1bc+k3a, k1*sqr_c+k2)

    def scale(self, s):
        self.mlist[0] *= s.x
        self.mlist[1] *= s.y
        self.mlist[2] *= s.z
        self.mlist[3] *= s.x
        self.mlist[4] *= s.y
        self.mlist[5] *= s.z
        self.mlist[6] *= s.x
        self.mlist[7] *= s.y
        self.mlist[8] *= s.z
        return self

    def rotate(self, angle, axis):
        R=self.rotation(angle, axis)
        self.mlist = (self*R).mlist
        return self

    def ortho(self):
        """Return a matrix with orthogonal base vectors.
        """

        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist

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

        return mat3( x.x, y.x, z.x,
                     x.y, y.y, z.y,
                     x.z, y.z, z.z)

    def decompose(self):
        """Decomposes the matrix into a rotation and scaling part.

        Returns a tuple (rotation, scaling). The scaling part is given
        as a vec3, the rotation is still a mat3.
        """
        dummy = self.ortho()

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

        return (dummy, scale)


######################################################################

if __name__=="__main__":

    vec3 = _vec3
    a=vec3(1,2,3)

    M = mat3("2,4,5,6")

    a=mat3(M)
    a[0,0]=17
    print(M)
    print(a)
