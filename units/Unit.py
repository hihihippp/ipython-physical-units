# -*- coding: utf-8 -*-
import numpy as np
from NDict import *

def isPhysicalUnit(x):
    return hasattr(x, 'factor') and hasattr(x, 'powers')

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
        unit = self.name.replace('**', '^').replace('�', 'mu').replace('�', 'deg')
        return '<PhysicalUnit ' + self.name + '>'

    def _repr_latex_(self):   
        unit = self.name.replace('**', '^').replace('�', 'mu').replace('�', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')
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

    def html_list(self):
        """ List all defined units """
        str = "<table>"
        str += "<tr><th>Name</th><th>Base Unit</th><th>Quantity</th></tr>"
        for name in _unit_table:
            unit = _unit_table[name]
            if isinstance(unit,PhysicalUnit):
                if unit.prefixed == False:
                    baseunit = '$ %s $' % unit.baseunit 
                    baseunit = baseunit.strip().replace('^', '**').replace('�', 'mu').replace('�', 'deg').replace('*', r' \cdot ').replace(' pi', r' \pi ')
                    str+= "<tr><td>" + unit.name + '</td><td>' + baseunit +\
                          '</td><td><a href="' + unit.url+'" target="_blank">'+ unit.comment +\
                          "</a></td></tr>"
        str += "</table>"
        return HTML(str)
        
    def list(self):
        """ List all defined units """
        str=[]
        for name in _unit_table:
            unit = _unit_table[name]
            if isinstance(unit,PhysicalUnit):
                str.append(unit.name)
        return str

