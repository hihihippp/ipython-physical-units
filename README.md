physical-units
--------------

Python module and IPython notebook extension to add physical units to Python, IPython and the IPython notebook.

For an example notebook see here:
http://nbviewer.ipython.org/github/juhasch/ipython-physical-units/blob/master/notebooks/units-example.ipynb

This extension is based on code from Georg Brandl here: https://bitbucket.org/birkenfeld/ipython-physics

Features
--------

- Provides standard units like m, s, g, A, K , mol, cd, rad, sr and derived units
- All units can be prefixed with a scaling factor: nm, mm, cm, km
- Optimized to work with the IPython notebook. Provides _repr_latex function to pretty print units
- Basic support of Numpy ndarrays
- Ignores units inside strings and quotation marks
- Complex number handling

TODO
----

- Cleanup
- Provide a more exhaustive list of physical constants
