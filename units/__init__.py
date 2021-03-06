# -*- coding: utf-8 -*-
"""
Physical Units for Python
"""
from __future__ import absolute_import
from math import pi
import numpy as np

from .Unit import *

class UnitError(ValueError):
    pass


def isPhysicalQuantity(x):
    return hasattr(x, 'value') and hasattr(x, 'unit')

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
        
        if unit is not None:
            self.value = value
            self.unit = _findUnit(unit)
        else:
            raise UnitError('No number found in %r' % value)

    def __getattr__(self,attr):
        """ Convert to different scaling in the same unit.
            If a '_' is appended, drop unit after rescaling and return value only.
        """
        dropunit = (attr[-1] == '_')
        attr = attr.strip('_')
        try:
            attrunit = _unit_table[attr]
        except:
            raise AttributeError
        if self.unit.prefixed == True:
            base = self.unit.comment
        else:
            base = self.unit.name
        if isPhysicalUnit(attrunit):
            if dropunit == True :
                return self.to(attrunit.name).value
            else:
                return self.to(attrunit.name)
        raise AttributeError

    def __len__(self):
        if isinstance(self.value, np.ndarray) or isinstance(self.value, list):
            return len(self.value)
        raise TypeError

    def __str__(self):
        unit = self.unit.name.replace('**', '^')
        return '%s %s' % ( self.value, unit)

    def __complex__(self):
        return self.base.value

    def __float__(self):
            return self.base.value

    def __repr__(self):
        return self.__str__()

    def _repr_latex_(self):
        unit = self.unit.name.replace('**', '^').replace('�', 'mu').replace('�', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')        
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
            return np.sin(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of sin must be an angle')

    def cos(self):
        if self.unit.is_angle:
            return np.cos(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of cos must be an angle')

    def tan(self):
        if self.unit.is_angle:
            return np.tan(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of tan must be an angle')

        
# Helper functions
def _findUnit(unit):
    if isinstance(unit, basestring):
        name = unit.strip().replace('^', '**').replace('�', 'mu').replace('�', 'deg') #.replace('*', r' \cdot ').replace(' pi', r' \pi ')
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



global _base_names
global _base_units

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


Q=PhysicalQuantity
U=PhysicalUnit


# numpy linspace wrapper for units
def Qlinspace(start, stop, num = 50,  endpoint=True, retstep=False):
    if not isinstance(start,PhysicalQuantity) and not isinstance(stop,PhysicalQuantity):
        return np.linspace(start, stop, num,  endpoint, retstep)

    if isinstance(start,PhysicalQuantity) and isinstance(stop,PhysicalQuantity):
        if start.base.unit != stop.base.unit:
            print start.base.unit, stop.base.unit
            raise UnitError("Cannot match units")
    
    if isinstance(start,PhysicalQuantity):
        start_value = start.value
        unit = start.unit
    else:
        start_value = start

    if isinstance(stop,PhysicalQuantity):
        stop_value = stop.value
        unit = stop.unit
    else:
        stop_value = stop

    array = np.linspace(start_value, stop_value, num,  endpoint, retstep)

    if retstep:
        return PhysicalQuantity(array[0], unit), PhysicalQuantity(array[1], unit)
    else:
        return array * PhysicalQuantity(1, unit)

