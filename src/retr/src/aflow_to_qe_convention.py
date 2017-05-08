# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

import functools

import numpy as np
from numpy import linalg as la


__all__ = ('aflow_prim2conv', 'aflow_conv2prim', 'qe_prim2conv', 'qe_conv2prim',
           'conv_aflow2qe', 'conv_qe2aflow', 'prim_aflow2qe', 'prim_qe2aflow',
           'InvalidIbravError')


class InvalidIBravError(ValueError):
    """Supplied Bravais lattice index (ibrav) is not valid.
    
    Attributes:
        value: The invalid Bravais lattice index that caused this error.
    """
    
    def __init__(self, value):
        self.value = value
        super(ValueError, self).__init__('%s is not a valid Bravais lattice' % value)


def validate_ibrav(func):
    """Decorates a basis transformation function to validate the Bravais lattice
    index in the range [0, 14].
    
    Arguments:
        func: The decorated function of two arguments, cell and ibrav.
    
    Returns:
        Decorated function that raises `InvalidIBravError` when bad values are
        encountered.
    """
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        
        ibrav = args[1]
        
        try:
            assert ibrav in range(15)
        
        except AssertionError:
            raise InvalidIBravError(ibrav)
        
        return func(*args, **kwargs)
    
    return decorator


def matmul(a, b, *others):
    """Calculates the accumulated matrix or dot product of the arguments."""
    args = (a, b) + others
    return functools.reduce(np.dot, args)


class IBrav(object):
    """ibrav namespace"""
    FREE = 0
    CUB = 1
    FCC = 2
    BCC = 3
    HEX = 4
    RHL = 5
    TET = 6
    BCT = 7
    ORC = 8
    ORCC = 9
    ORCF = 10
    ORCI = 11
    MCL = 12
    MCLC = 13
    TRI = 14


def aflow_primitive_to_aflow_conventional_transform(ibrav):
    """
    Builds a transformation matrix for converting primtive cell vectors to
    conventional cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # face-centered lattices
    if ibrav in (IBrav.FCC, IBrav.ORCF):
        transform = np.matrix([[-1.,  1.,  1.],
                               [ 1., -1.,  1.],
                               [ 1.,  1., -1.]])
        
    # body-centered lattices
    elif ibrav in (IBrav.BCC, IBrav.BCT, IBrav.ORCI):
        transform = np.matrix([[ 0.,  1.,  1.],
                               [ 1.,  0.,  1.],
                               [ 1.,  1.,  0.]])
        
    # orthorhombic base-centered lattice
    elif ibrav == IBrav.ORCC:
        transform = np.matrix([[ 1.,  1.,  0.],
                               [-1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # monoclinic base-centered lattice
    elif ibrav == IBrav.MCLC:
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # lattices without centering share the same cell vectors between
    # their primitive cells and conventional cells, so no change is necessary
    else:
        transform = np.matrix(np.identity(3))
    
    return transform
    
    
def aflow_conventional_to_aflow_primitive_transform(ibrav):
    """
    Builds a transformation matrix for converting conventional cell vectors to
    primitive cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # primitive to conventional is the inverse of the transform we want to do
    inverse_transform = aflow_primitive_to_aflow_conventional_transform(ibrav)
    
    # invert the inverse transform
    return la.inv(inverse_transform)


def qe_primitive_to_qe_conventional_transform(ibrav):
    """
    Builds a transformation matrix for converting primitive cell vectors to
    conventional cell vectors using Quantum Espresso convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    if ibrav == IBrav.FCC:
        transform = np.matrix([[-1.,  1., -1.],
                               [-1.,  1.,  1.],
                               [ 1.,  1., -1.]])
        
    # body-centered lattices
    elif ibrav in (IBrav.BCC, IBrav.ORCI):
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 0.,  1., -1.],
                               [ 1.,  0.,  1.]])
        
    # can't stuff BCT in with the other body-centered lattices
    # like in AFLOW because QE uniquely defines its lattice vectors
    elif ibrav == IBrav.BCT:
        transform = np.matrix([[ 1.,  0., -1.],
                               [-1.,  1.,  0.],
                               [ 0.,  1.,  1.]])
    
    # orthorhombic base-centered lattice
    elif ibrav == IBrav.ORCC:
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # orthorhombic face-centered lattice
    elif ibrav == IBrav.ORCF:
        transform = np.matrix([[ 1.,  1., -1.],
                               [-1.,  1.,  1.],
                               [ 1., -1.,  1.]])
    
    # monoclinic base-centered lattice
    elif ibrav == IBrav.MCLC:
        transform = np.matrix([[ 1.,  0.,  1.],
                               [ 0.,  1.,  0.],
                               [-1.,  0.,  1.]])
        
    # lattices without centering share the same cell vectors between
    # their primitive cells and conventional cells, so no change is necessary
    else:
        transform = np.matrix(np.identity(3))
    
    return transform


def qe_conventional_to_qe_primitive_transform(ibrav):
    """
    Builds a transformation matrix for converting conventional cell vectors to
    primitive cell vectors using Quantum Espresso convention.
    """
    # primitive to conventional is the inverse of the transform we want to do
    inverse_transform = qe_primitive_to_qe_conventional_transform(ibrav)
    
    # invert the inverse transform
    return la.inv(inverse_transform)


def aflow_conventional_to_qe_conventional_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting conventional cell vectors
    using AFLOW convention to conventional cell vectors using Quantum Espresso
    convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    if ibrav == IBrav.HEX:
        transform = np.matrix([[ 1.,  1.,  0.],
                               [-1.,  0.,  0.],
                               [ 0.,  0.,  1.]])
    
    elif ibrav == IBrav.RHL:
        raise NotImplementedError
        
        assert angle is not None, "characteristic angle must be supplied for RHL"
        
        cos = np.cos(angle)
        cos_half = np.cos(angle / 2)
        sin_half = np.sin(angle / 2)
        
        aflow = np.matrix([
            [      cos_half, -sin_half,                                0],
            [      cos_half,  sin_half,                                0],
            [cos / cos_half,         0, np.sqrt(1 - (cos / cos_half)**2)]
        ])
        
        tx = np.sqrt((1 - cos) / 2)
        ty = np.sqrt((1 - cos) / 6)
        tz = np.sqrt((1 + 2*cos) / 3)
        
        qe = np.matrix([
            [ tx,  -ty, tz],
            [  0, 2*ty, tz],
            [-tx,  -ty, tz]
        ])
        
        transform = matmul(qe, la.inv(aflow))
    
    elif ibrav in (IBrav.MCL, IBrav.MCLC):
        raise NotImplementedError
    
    else:
        transform = np.identity(3)
    
    return transform


def qe_conventional_to_aflow_conventional_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting conventional cell vectors
    using Quantum Espresso convention to conventional cell vectors using AFLOW
    convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    inverse_transform = aflow_conventional_to_qe_conventional_transform(ibrav, angle)
    return la.inv(inverse_transform)


def aflow_primitive_to_qe_primitive_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting primitive cell vectors using
    AFLOW convention to primitive cell vectors using Quantum Espresso convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # successive transforms AFLOW prim -> AFLOW conv -> QE conv -> QE prim
    transforms = (aflow_primitive_to_aflow_conventional_transform(ibrav),
                  aflow_conventional_to_qe_conventional_transform(ibrav, angle),
                  qe_conventional_to_qe_primitive_transform(ibrav))
    
    return matmul(*reversed(transforms))


def qe_primitive_to_aflow_primitive_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting primitive cell vectors using
    Quantum Espresso convention to primitive cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # successive transforms QE prim -> QE conv -> AFLOW conv -> AFLOW prim
    transforms = (qe_primitive_to_qe_conventional_transform(ibrav),
                  qe_conventional_to_aflow_conventional_transform(ibrav, angle),
                  aflow_conventional_to_aflow_primitive_transform(ibrav))
    
    return matmul(*reversed(transforms))


@validate_ibrav
def aflow_primitive_to_aflow_conventional(cell, ibrav):
    """
    Converts primitive cell vectors to conventional cell vectors using AFLOW
    convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW conventional lattice vectors.
    """
    transform = aflow_primitive_to_aflow_conventional_transform(ibrav)
    return matmul(transform, cell)
    

@validate_ibrav
def aflow_conventional_to_aflow_primitive(cell, ibrav):
    """
    Converts conventional cell vectors to primitive cell vectors using AFLOW
    convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW primitive lattice vectors.
    """
    transform = aflow_conventional_to_aflow_primitive_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def qe_primitive_to_qe_conventional(cell, ibrav):
    """
    Converts primitive cell vectors to conventional cell vectors using Quantum
    Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_primitive_to_qe_conventional_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def qe_conventional_to_qe_primitive(cell, ibrav):
    """
    Converts conventional cell vectors to primitive cell vectors using Quantum
    Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_conventional_to_qe_primitive_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def aflow_conventional_to_qe_conventional(cell, ibrav, angle=None):
    """
    Converts conventional cell vectors using AFLOW convention to conventional
    cell vectors using Quantum Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso conventional lattice
        vectors.
    """
    transform = aflow_conventional_to_qe_conventional_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def qe_conventional_to_aflow_conventional(cell, ibrav, angle=None):
    """
    Converts conventional cell vectors using Quantum Espresso convention to
    conventional cell vectors using AFLOW convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW conventional lattice vectors.
    """
    transform = qe_conventional_to_aflow_conventional_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def aflow_primitive_to_qe_primitive(cell, ibrav, angle=None):
    """
    Converts primitive cell vectors using AFLOW convention to primitive cell
    vectors using Quantum Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = aflow_primitive_to_qe_primitive_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def qe_primitive_to_aflow_primitive(cell, ibrav, angle=None):
    """
    Converts primitive cell vectors using Quantum Espresso convention to
    primitive cell vectors using AFLOW convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_primitive_to_aflow_primitive_transform(ibrav, angle)
    return matmul(transform, cell)


# aliases
aflow_prim2conv = aflow_primitive_to_aflow_conventional
aflow_conv2prim = aflow_conventional_to_aflow_primitive
conv_aflow2qe = aflow_conventional_to_qe_conventional
conv_qe2aflow = qe_conventional_to_aflow_conventional
qe_prim2conv = qe_primitive_to_qe_conventional
qe_conv2prim = qe_conventional_to_qe_primitive
prim_aflow2qe = aflow2qe = aflow_primitive_to_qe_primitive
prim_qe2aflow = qe2aflow = qe_primitive_to_aflow_primitive

