ipython-physical-units
----------------------

IPython notebook extension to add physical units to the IPython notebook.
This extension can also be used without the notebook.

For an example notebook see here:
[https://github.com/juhasch/ipython-physical-units/raw/master/notebooks/140126-units-example.ipynb]


This extension is based on code from Georg Brandl here:
[https://bitbucket.org/birkenfeld/ipython-physics]

Additional features
-------------------

- _repr_ function to display units as TeX
- ignores units inside strings / quotation marks
- works with ndarrays

TODO:
- I removed some of the units included in the original implementation, add them back.
- Provide a more exhaustive list of physical constants
- Improve interaction with the uncertainties module 
- Improve numpy handling. Automatically convert back to built-in type if units cancel out
- Complex number handling ?