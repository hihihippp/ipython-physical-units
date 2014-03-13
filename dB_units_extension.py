""" IPython extension for dB calculations """

#-----------------------------------------------------------------------------
#  Copyright (C) 2013  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

import re
import numpy as np
import IPython
ip = IPython.get_ipython()

from IPython.core.inputtransformer import StatelessInputTransformer, CoroutineInputTransformer

import units.dB_units as dB

# regex: number + space + dB-unit
# valid: 0dBm, 0 dBm, 0. dBm
# invalid: 0.dBm
number = r'(-?[\d0-9-.]+)'
_unit_list = '('
for x in dB._dB_units.keys()[0:-1]:
    _unit_list += x + '|'
_unit_list += dB._dB_units.keys()[-1] + ')'
match = number + r'(\s*)' + _unit_list

# regex for finding units and quoted strings
number = r'-?[\d0-9.]+'  # problem: 1dBm-0dB
stringmatch = r'(["\'])(?:(?=(\\?))\2.)*?\1'
match = stringmatch + '|' + number + r'(\s*)' + _unit_list
line_match = re.compile(match)

# regex to match unit after it has been found using line_match
number = r'(-?[\d0-9-]+)'
match = number + r'(.\s|\s*)' + _unit_list
unit_match = re.compile(match)


def replace_inline(ml):
    """Replace an inline unit expression by valid Python code
    """
    if ml.group()[0][0] in '"\'':
        return ml.group()

    def replace_unit(mo):
        try:
            return "dBUnit(" + mo.group(1) + ", '" + mo.group(3) + "', islog=True)"
        except KeyError:
            return mo.group()

    return unit_match.sub(replace_unit, ml.group())


@StatelessInputTransformer.wrap
def dBUnit_transform(line):
    line = line_match.sub(replace_inline, line)
    return line

__transformer = dBUnit_transform()


def load_ipython_extension(ip):
    global __transformer
    ip.input_transformer_manager.logical_line_transforms.insert(0, __transformer)
    ip.user_ns['dBUnit'] = dB.dBUnit


def unload_ipython_extension(ip):
    global __transformer
    if type(__transformer) is StatelessInputTransformer:
        ip.input_transformer_manager.logical_line_transforms.remove(__transformer)
        ip.user_ns.pop('dBUnit')
