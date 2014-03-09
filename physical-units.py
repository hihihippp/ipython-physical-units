# -*- coding: utf-8 -*-
""" IPython extension for physical quantity input """

#-----------------------------------------------------------------------------
#  Copyright (C) 2013  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------
# Original author: Georg Brandl <georg@python.org>.
#                  https://bitbucket.org/birkenfeld/ipython-physics

import re
import sys
from math import pi
import numpy as np

from IPython.core.inputtransformer import StatelessInputTransformer
from IPython.core.inputtransformer import CoroutineInputTransformer
from IPython.display import display, Math, Latex, HTML

class UnitError(ValueError):
    pass

# Adapted from ScientificPython:
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# with contributions from Greg Ward
# last revision: 2007-5-25


class NumberDict(dict):
    """Dictionary storing numerical values.

    An instance of this class acts like an array of number with generalized
    (non-integer) indices. A value of zero is assumed for undefined
    entries. NumberDict instances support addition, and subtraction with other
    NumberDict instances, and multiplication and division by scalars.
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return 0

    def __coerce__(self, other):
        if isinstance(other, dict):
            other = NumberDict(other)
        return self, other

    def __add__(self, other):
        sum_dict = NumberDict()
        for key in self.keys():
            sum_dict[key] = self[key]
        for key in other.keys():
            sum_dict[key] = sum_dict[key] + other[key]
        return sum_dict

    def __sub__(self, other):
        sum_dict = NumberDict()
        for key in self.keys():
            sum_dict[key] = self[key]
        for key in other.keys():
            sum_dict[key] = sum_dict[key] - other[key]
        return sum_dict

    def __mul__(self, other):
        new = NumberDict()
        for key in self.keys():
            new[key] = other*self[key]
        return new
    __rmul__ = __mul__

    def __div__(self, other):
        new = NumberDict()
        for key in self.keys():
            new[key] = self[key]/other
        return new


# Type checks


def isPhysicalUnit(x):
    return hasattr(x, 'factor') and hasattr(x, 'powers')


def isPhysicalQuantity(x):
    return hasattr(x, 'value') and hasattr(x, 'unit')
   

class PhysicalUnit(object):
    """Physical unit.

    A physical unit is defined by a name (possibly composite), a scaling factor,
    and the exponentials of each of the SI base units that enter into it. Units
    can be multiplied, divided, and raised to integer powers.
    """

    def __init__(self, names, factor, powers, offset=0,url='',comment=''):
        """
        @param names: a dictionary mapping each name component to its
                      associated integer power (e.g. C{{'m': 1, 's': -1}})
                      for M{m/s}). As a shorthand, a string may be passed
                      which is assigned an implicit power 1.
        @param factor: a scaling factor
        @param powers: the integer powers for each of the nine base units
        @param offset: an additive offset to the base unit (used only for
                       temperatures)
        """

        self.prefixed = False
        self.prefixunit = self
        self.baseunit = self
        self.comment = comment
        self.url = url        
        if isinstance(names, basestring):
            self.names = NumberDict()
            self.names[names] = 1
        else:
            self.names = names
        self.factor = factor
        self.offset = offset
        self.powers = powers

    def set_name(self, name):
        self.names = NumberDict()
        self.names[name] = 1

    @property
    def name(self):
        num = ''
        denom = ''
        for unit in self.names.keys():
            power = self.names[unit]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return num + denom

    @property
    def is_dimensionless(self):
        return not reduce(lambda a, b: a or b, self.powers)

    @property
    def is_angle(self):
        return self.powers[7] == 1 and \
            reduce(lambda a, b: a + b, self.powers) == 1

    def __str__(self):
        return self.name

    def __repr__(self):
        unit = self.name.replace('**', '^').replace('µ', 'mu').replace('°', 'deg')
        return '<PhysicalUnit ' + self.name + '>'

    def _repr_latex_(self):   
        unit = self.name.replace('**', '^').replace('µ', 'mu').replace('°', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')
        if self.prefixed == False:
            info = '(<a href="' + self.url + '" target="_blank">'+ self.comment + '</a>)'
        else:
            baseunit = _unit_table[self.prefixunit]
#            print  baseunit.comment, baseunit.url
            info = r'$ = %s \cdot %s$ (' % (self.factor, self.prefixunit) +\
                    '<a href="' + baseunit.url + '" target="_blank">'+ baseunit.comment + '</a>)'            
        s = r'$%s$ %s' % (unit, info)
        return s

    @property
    def latex(self):
        return Latex(self._repr_latex_())

    def __cmp__(self, other):
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        return cmp(self.factor, other.factor)

    def __mul__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot multiply units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names + other.names,
                                self.factor * other.factor,
                                map(lambda a, b: a+b, self.powers, other.powers))
        else:
            return PhysicalUnit(self.names + {str(other): 1},
                                self.factor * other, self.powers,
                                self.offset * other)

    __rmul__ = __mul__

    def __div__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names - other.names,
                                self.factor / other.factor,
                                map(lambda a, b: a-b, self.powers, other.powers))
        else:
            return PhysicalUnit(self.names+{str(other): -1},
                                self.factor/other, self.powers)

    def __rdiv__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(other.names - self.names,
                                other.factor/self.factor,
                                map(lambda a, b: a-b, other.powers, self.powers))
        else:
            return PhysicalUnit({str(other): 1} - self.names,
                                other / self.factor,
                                map(lambda x: -x, self.powers))

    def __pow__(self, other):
        if self.offset != 0:
            raise UnitError('Cannot exponentiate units with non-zero offset')
        if isinstance(other, int):
            return PhysicalUnit(other*self.names, pow(self.factor, other),
                                map(lambda x, p=other: x*p, self.powers))
        if isinstance(other, float):
            inv_exp = 1./other
            rounded = int(np.floor(inv_exp + 0.5))
            if abs(inv_exp-rounded) < 1.e-10:
                if reduce(lambda a, b: a and b,
                          map(lambda x, e=rounded: x%e == 0, self.powers)):
                    f = pow(self.factor, other)
                    p = map(lambda x, p=rounded: x/p, self.powers)
                    if reduce(lambda a, b: a and b,
                              map(lambda x, e=rounded: x%e == 0,
                                  self.names.values())):
                        names = self.names/rounded
                    else:
                        names = NumberDict()
                        if f != 1.:
                            names[str(f)] = 1
                        for i in range(len(p)):
                            names[_base_names[i]] = p[i]
                    return PhysicalUnit(names, f, p)
                else:
                    raise UnitError('Illegal exponent %f' % other)
        raise UnitError('Only integer and inverse integer exponents allowed')

    def conversion_factor_to(self, other):
        """Return conversion factor to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        if self.offset != other.offset and self.factor != other.factor:
            raise UnitError(('Unit conversion (%s to %s) cannot be expressed ' +
                            'as a simple multiplicative factor') %
                            (self.name, other.name))
        return self.factor/other.factor

    def conversion_tuple_to(self, other):
        """Return conversion factor and offset to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')

        # let (s1,d1) be the conversion tuple from 'self' to base units
        #   (ie. (x+d1)*s1 converts a value x from 'self' to base units,
        #   and (x/s1)-d1 converts x from base to 'self' units)
        # and (s2,d2) be the conversion tuple from 'other' to base units
        # then we want to compute the conversion tuple (S,D) from
        #   'self' to 'other' such that (x+D)*S converts x from 'self'
        #   units to 'other' units
        # the formula to convert x from 'self' to 'other' units via the
        #   base units is (by definition of the conversion tuples):
        #     ( ((x+d1)*s1) / s2 ) - d2
        #   = ( (x+d1) * s1/s2) - d2
        #   = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
        #   = ( (x+d1) - (d1*s2/s1) ) * s1/s2
        #   = (x + d1 - d2*s2/s1) * s1/s2
        # thus, D = d1 - d2*s2/s1 and S = s1/s2
        factor = self.factor / other.factor
        offset = self.offset - (other.offset * other.factor / self.factor)
        return (factor, offset)

    def list(self):
        """ List all defined units """
        str = "<table>"
        str += "<tr><th>Name</th><th>Base Unit</th><th>Quantity</th></tr>"
        for name in _unit_table:
            unit = _unit_table[name]
            if isinstance(unit,PhysicalUnit):
                if unit.prefixed == False:
                    baseunit = '$ %s $' % unit.baseunit 
                    baseunit = baseunit.strip().replace('^', '**').replace('µ', 'mu').replace('°', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')
                    str+= "<tr><td>" + unit.name + '</td><td>' + baseunit +\
                          '</td><td><a href="' + unit.url+'" target="_blank">'+ unit.comment +\
                          "</a></td></tr>"
        str += "</table>"
        return HTML(str)


# Helper functions
def _findUnit(unit):
    if isinstance(unit, basestring):
        name = unit.strip().replace('^', '**').replace('µ', 'mu').replace('°', 'deg') #.replace('*', r' \cdot ').replace(' pi', r' \pi ')
        try:
            unit = eval(name, _unit_table)
        except NameError:
            raise UnitError('Invalid or unknown unit in %s' % name)
        for cruft in ['__builtins__', '__args__']:
            try:
                del _unit_table[cruft]
            except:
                pass
    if not isPhysicalUnit(unit):
        raise UnitError(str(unit) + ' is not a unit')
    return unit


def _convertValue(value, src_unit, target_unit):
    (factor, offset) = src_unit.conversion_tuple_to(target_unit)
    return (value + offset) * factor


class PhysicalQuantity(object):
    """Physical quantity with units.

    PhysicalQuantity instances allow addition, subtraction, multiplication, and
    division with each other as well as multiplication, division, and
    exponentiation with numbers.  Addition and subtraction check that the units
    of the two operands are compatible and return the result in the units of the
    first operand. A limited set of mathematical functions (from numpy) is
    applicable as well.
    """

    __array_priority__ = 1000 # make sure numpy arrays do not get iterated
    
    def __init__(self, value, unit=None,  **kwargs):
        """There are two constructor calling patterns:

        1. PhysicalQuantity(value, unit), where value is any number and unit is
           a string defining the unit

        2. PhysicalQuantity(value_with_unit), where value_with_unit is a string
           that contains both the value and the unit, i.e. '1.5 m/s'. This form
           is provided for more convenient interactive use.
        """
        for key, val in kwargs.iteritems():
            pass
#            if key is 'islog':
#                islog = val    # convert to log at initialization
        
        if unit is not None:
            self.value = value
            self.unit = _findUnit(unit)
        else:
            raise UnitError('No number found in %r' % value)

    def __getattr__(self,attr):
        try:
            attrunit = _unit_table[attr]
        except:
            raise AttributeError
        if self.unit.prefixed == True:
            base = self.unit.comment
        else:
            base = self.unit.name
        if isPhysicalUnit(attrunit):
            return self.to(attrunit.name)
        raise AttributeError

    def __str__(self):
        unit = self.unit.name.replace('**', '^')
        return '%s %s' % ( self.value, unit)

    def __complex__(self):
        return self.base.value

    def __float__(self):
            return self.base.value

    @property
    def Q(self):
        return self.__float__()

    def __repr__(self):
        return self.__str__()

    def _repr_latex_(self):
        unit = self.unit.name.replace('**', '^').replace('µ', 'mu').replace('°', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')        
        s = r'%s $%s$' % (self.value, unit)
        return s

    @property
    def latex(self):
        return Latex(self._repr_latex_())

    def _sum(self, other, sign1, sign2):
        if not isPhysicalQuantity(other):
            raise UnitError('Incompatible types')
        new_value = sign1 * self.value + \
            sign2 * other.value * other.unit.conversion_factor_to(self.unit)
        return self.__class__(new_value, self.unit)

    def __add__(self, other):
        return self._sum(other, 1, 1)

    __radd__ = __add__

    def __sub__(self, other):
        return self._sum(other, 1, -1)

    def __rsub__(self, other):
        return self._sum(other, -1, 1)

    def __cmp__(self, other):
        diff = self._sum(other, 1, -1)
        return cmp(diff.value, 0)

    def __mul__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value * other, self.unit)
        value = self.value * other.value
        unit = self.unit * other.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    __rmul__ = __mul__

    def __div__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value / other, self.unit)
        value = self.value / other.value
        unit = self.unit / other.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    def __rdiv__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(other / self.value, pow(self.unit, -1))
        value = other.value / self.value
        unit = other.unit / self.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        if isPhysicalQuantity(other):
            raise UnitError('Exponents must be dimensionless')
        return self.__class__(pow(self.value, other), pow(self.unit, other))

    def __rpow__(self, other):
        raise UnitError('Exponents must be dimensionless')

    def __abs__(self):
        return self.__class__(abs(self.value), self.unit)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value, self.unit)

    def __nonzero__(self):
        return self.value != 0

    def __format__(self, *args, **kw):
        return "{1:{0}} {2}".format(args[0], self.value, self.unit)

    def convert(self, unit):
        """Change the unit and adjust the value such that the combination is
        equivalent to the original one. The new unit must be compatible with the
        previous unit of the object.
        """
        unit = _findUnit(unit)
        self.value = _convertValue(self.value, self.unit, unit)
        self.unit = unit

    def _round(self, x):
        if np.greater(x, 0.):
            return np.floor(x)
        else:
            return np.ceil(x)

    def to(self, *units):
        """Express the quantity in different units. If one unit is specified, a
        new PhysicalQuantity object is returned that expresses the quantity in
        that unit. If several units are specified, the return value is a tuple
        of PhysicalObject instances with with one element per unit such that the
        sum of all quantities in the tuple equals the the original quantity and
        all the values except for the last one are integers. This is used to
        convert to irregular unit systems like hour/minute/second.
        """
        units = map(_findUnit, units)
        if len(units) == 1:
            unit = units[0]
            value = _convertValue(self.value, self.unit, unit)
            return self.__class__(value, unit)
        else:
            units.sort()
            result = []
            value = self.value
            unit = self.unit
            for i in range(len(units)-1, -1, -1):
                value = value*unit.conversion_factor_to(units[i])
                if i == 0:
                    rounded = value
                else:
                    rounded = self._round(value)
                result.append(self.__class__(rounded, units[i]))
                value = value - rounded
                unit = units[i]
            return tuple(result)

    @staticmethod
    def any_to(qty, unit):
        if not isPhysicalQuantity(qty):
            qty = PhysicalQuantity(qty, 'rad')
        return qty.to(unit)

    @property
    def base(self):
        """Returns the same quantity converted to base units."""
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in xrange(len(_base_names)): 
            unit = _base_names[i]
            power = self.unit.powers[i]
            if power < 0:
                denom += '/' + unit
                if power < -1:
                    denom += '**' + str(-power)
            elif power > 0:
                num += '*' + unit
                if power > 1:
                    num += '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return self.__class__(new_value, num + denom)

    # implementations of special functions, used by numpy ufuncs
    def sqrt(self):
        return pow(self, 0.5)

    def sin(self):
        if self.unit.is_angle:
            return unp.sin(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of sin must be an angle')

    def cos(self):
        if self.unit.is_angle:
            return unp.cos(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of cos must be an angle')

    def tan(self):
        if self.unit.is_angle:
            return unp.tan(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of tan must be an angle')


_base_names = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [
    ('m',   PhysicalUnit('m',   1.,    [1, 0, 0, 0, 0, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Metre', comment='Metre')),
    ('g',   PhysicalUnit('g',   0.001, [0, 1, 0, 0, 0, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Kilogram', comment='Kilogram')),
    ('s',   PhysicalUnit('s',   1.,    [0, 0, 1, 0, 0, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Second', comment='Second')),
    ('A',   PhysicalUnit('A',   1.,    [0, 0, 0, 1, 0, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Ampere', comment='Ampere')),
    ('K',   PhysicalUnit('K',   1.,    [0, 0, 0, 0, 1, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Kelvin', comment='Kelvin')),
    ('mol', PhysicalUnit('mol', 1.,    [0, 0, 0, 0, 0, 1, 0, 0, 0],url='https://en.wikipedia.org/wiki/Mole_(unit)', comment='Mol')),
    ('cd',  PhysicalUnit('cd',  1.,    [0, 0, 0, 0, 0, 0, 1, 0, 0],url='https://en.wikipedia.org/wiki/Candela', comment='Candela')),
    ('rad', PhysicalUnit('rad', 1.,    [0, 0, 0, 0, 0, 0, 0, 1, 0],url='https://en.wikipedia.org/wiki/Radian', comment='Radian')),
    ('sr',  PhysicalUnit('sr',  1.,    [0, 0, 0, 0, 0, 0, 0, 0, 1],url='https://en.wikipedia.org/wiki/Steradian', comment='Streradian')), ]


_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]


def _addUnit(name, unit, comment='',prefixed=False, prefixunit=None, url=''):
    if name in _unit_table:
        raise KeyError('Unit ' + name + ' already defined')
    if isinstance(unit, str):
        newunit = eval(unit, _unit_table)
        for cruft in ['__builtins__', '__args__']:
            try:
                del _unit_table[cruft]
            except:
                pass
    else:
        newunit = unit
    newunit.set_name(name)
    newunit.comment = comment
    newunit.prefixunit = prefixunit
    newunit.baseunit = unit
    newunit.prefixed = prefixed
    newunit.url = url
    _unit_table[name] = newunit


def _addPrefixed(unit):
    _prefixed_names = []
    for prefix in _prefixes:
        name = prefix[0] + unit
        _addUnit(name, prefix[1]*_unit_table[unit],prefixed=True,prefixunit=unit)
        _prefixed_names.append(name)

_unit_table['kg'] = PhysicalUnit('kg',   1., [0, 1, 0, 0, 0, 0, 0, 0, 0])
_addUnit('Hz', '1/s', 'Hertz', url='https://en.wikipedia.org/wiki/Hertz')
_addUnit('N', 'm*kg/s**2', 'Newton', url='https://en.wikipedia.org/wiki/Newton_(unit)')
_addUnit('Pa', 'N/m**2', 'Pascal', url='https://en.wikipedia.org/wiki/Pascal_(unit)')
_addUnit('J', 'N*m', 'Joule', url='https://en.wikipedia.org/wiki/Joule')
_addUnit('W', 'J/s', 'Watt', url='https://en.wikipedia.org/wiki/Watt')
_addUnit('C', 's*A', 'Coulomb', url='https://en.wikipedia.org/wiki/Coulomb')
_addUnit('V', 'W/A', 'Volt', url='https://en.wikipedia.org/wiki/Volt')
_addUnit('F', 'C/V', 'Farad', url='https://en.wikipedia.org/wiki/Farad')
_addUnit('Ohm', 'V/A', 'Ohm', url='https://en.wikipedia.org/wiki/Ohm_(unit)')
_addUnit('S', 'A/V', 'Siemens', url='https://en.wikipedia.org/wiki/Siemens_(unit)')
_addUnit('Wb', 'V*s', 'Weber', url='https://en.wikipedia.org/wiki/Weber_(unit)')
_addUnit('T', 'Wb/m**2', 'Tesla', url='https://en.wikipedia.org/wiki/Tesla_(unit)')
_addUnit('H', 'Wb/A', 'Henry', url='https://en.wikipedia.org/wiki/Henry_(unit)')
_addUnit('lm', 'cd*sr', 'Lumen', url='https://en.wikipedia.org/wiki/Lumen_(unit)')
_addUnit('lx', 'lm/m**2', 'Lux', url='https://en.wikipedia.org/wiki/Lux')
del _unit_table['kg']

# add scaling prefixes
_prefixes = [
    ('Y',  1.e24), ('Z',  1.e21), ('E',  1.e18), ('P',  1.e15), ('T',  1.e12),
    ('G',  1.e9),  ('M',  1.e6),  ('k',  1.e3),  ('h',  1.e2),  ('da', 1.e1),
    ('d',  1.e-1), ('c',  1.e-2), ('m',  1.e-3), ('mu', 1.e-6), ('n',  1.e-9),
    ('p',  1.e-12), ('f',  1.e-15), ('a',  1.e-18), ('z',  1.e-21),
    ('y',  1.e-24),
]

for unit in _unit_table.keys():
    _addPrefixed(unit)

_unit_table['pi'] = np.pi
_unit_table['kg'] = PhysicalUnit('kg',   1., [0, 1, 0, 0, 0, 0, 0, 0, 0],url='https://en.wikipedia.org/wiki/Kilogram')

# Angle units
_addUnit('deg', 'pi*rad/180', 'degrees')
_addUnit('arcmin', 'pi*rad/180/60', 'minutes of arc')
_addUnit('arcsec', 'pi*rad/180/3600', 'seconds of arc')
_unit_table['cycles'] = 2*np.pi

_addUnit('min', '60*s', 'Minute', url='https://en.wikipedia.org/wiki/Hour')
_addUnit('h', '60*60*s', 'Hour', url='https://en.wikipedia.org/wiki/Hour')

name = r'([_a-zA-Z]\w*)'
number = r'(-?[\d0-9.eE-]+)'
unit = r'([a-zA-Z1°µ][a-zA-Z0-9°µ/*^-]*)'
quantity = number + r'(?:\s+\+\/-\s+' + number + ')?' + r'\s+' + unit

inline_unit_re = re.compile(r'\((%s)\)' % quantity)

nice_assign_re = re.compile(r'^%s\s*=\s*(%s)$' % (name, quantity))
quantity_re = re.compile(quantity)
subst_re = re.compile(r'\?' + name)

# sort units after length for Regex matching
_li = sorted(_unit_table.iterkeys(), key=len)
_li.sort(key=len)
_unit_list = '('
for unit in _li[::-1]:
    _unit_list += unit + '|'
_unit_list = _unit_list[0:-1] + ')'

# regex for finding units and quoted strings
number = r'(\s+|=|^)(-?[\d0-9-]+' +r'-?[\d0-9.eE-]*)'
stringmatch = r'(["\'])(?:(?=(\\?))\2.)*?\1'
match = stringmatch + '|' + number + r'(\s*)' + _unit_list + r'(?:\W+|$)'
line_match = re.compile(match)

# regex to match unit after it has been found using line_match
number = r'(-?[\d0-9-]+' +r'-?[\d0-9.eE-]*)'

match = number + r'(.\s|\s*)' + _unit_list
unit_match = re.compile(match)


def replace_inline(match):
    """Replace an inline unit expression, e.g. ``(1 m)``, by valid Python code
    using a Quantity call.
    """
    return '(Quantity(\'' + match.group(1) + '\'))'


def replace_assign(match):
    """Replace a pretty assignment, e.g. ``B = 1 T``, by valid Python code using
    a Quantity call.
    """
    return '%s = Quantity(\'%s\')' % (match.group(1), match.group(2))


def replace_inline(ml):
    """Replace an inline unit expression by valid Python code
    """
    if ml.group()[0][0] in '"\'':
        return ml.group()

    def replace_unit(mo):
        try:
            return '(' + mo.group(1) + '*' + 'Quantity( 1,\'' + mo.group(3) + '\'))'
        except KeyError:
            return mo.group()
    return unit_match.sub(replace_unit, ml.group())


@StatelessInputTransformer.wrap
def _transform(line):
    line = line_match.sub(replace_inline, line)
    return line

__transformer = _transform()


def load_ipython_extension(ip):
    global __transformer
    ip.input_transformer_manager.logical_line_transforms.insert(0, __transformer)

    # set up simplified quantity input
    ip.user_ns['Quantity'] = PhysicalQuantity
       
    # active true float division
    exec ip.compile('from __future__ import division', '<input>', 'single') \
        in ip.user_ns


def unload_ipython_extension(ip):
    global __transformer
    if type(__transformer) is StatelessInputTransformer:
        ip.input_transformer_manager.logical_line_transforms.remove(__transformer)
        ip.user_ns.pop('Quantity')
