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

#from units import *
import units 


name = r'([_a-zA-Z]\w*)'
number = r'(-?[\d0-9.eE-]+)'
unit = r'([a-zA-Z1°µ][a-zA-Z0-9°µ/*^-]*)'
quantity = number + r'(?:\s+\+\/-\s+' + number + ')?' + r'\s+' + unit

inline_unit_re = re.compile(r'\((%s)\)' % quantity)

nice_assign_re = re.compile(r'^%s\s*=\s*(%s)$' % (name, quantity))
quantity_re = re.compile(quantity)
subst_re = re.compile(r'\?' + name)

# sort units after length for Regex matching
_li = sorted(units._unit_table.iterkeys(), key=len)
_li.sort(key=len)
_unit_list = '('
for unit in _li[::-1]:
    _unit_list += unit + '|'
_unit_list = _unit_list[0:-1] + ')'

# regex for finding units and quoted strings
#number = r'(\s+|=|^)(-?[\d0-9-]+' +r'-?[\d0-9.eE-]*)'
number = r'(-?[\d0-9-]+' +r'-?[\d0-9.eE-]*)'
stringmatch = r'(["\'])(?:(?=(\\?))\2.)*?\1'
match = stringmatch + '|' + number + r'(\s*)' +_unit_list + r'(?:\W+|$)'
line_match = re.compile(match)

# regex to match unit after it has been found using line_match
number = r'(-?[\d0-9-]+' +r'-?[\d0-9.eE-]*)'

match = number + r'(.\s|\s*)' + _unit_list
unit_match = re.compile(match)

def replace_inline(ml):
    """Replace an inline unit expression by valid Python code
    """
    if ml.group()[0][0] in '"\'':
        return ml.group()

    def replace_unit(mo):
        try:
            return '(' + mo.group(1) + '*' + 'Q( 1,\'' + mo.group(3) + '\'))'
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
    ip.user_ns['Q'] = units.PhysicalQuantity
       
    # active true float division
    exec ip.compile('from __future__ import division', '<input>', 'single') \
        in ip.user_ns


def unload_ipython_extension(ip):
    global __transformer
    if type(__transformer) is StatelessInputTransformer:
        ip.input_transformer_manager.logical_line_transforms.remove(__transformer)
        ip.user_ns.pop('Q')
