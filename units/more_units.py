
# SI derived units; these automatically get prefixes
_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])

_addUnit('Hz', '1/s', 'Hertz')
_addUnit('N', 'm*kg/s**2', 'Newton')
_addUnit('Pa', 'N/m**2', 'Pascal')
_addUnit('J', 'N*m', 'Joule')
_addUnit('W', 'J/s', 'Watt')
_addUnit('C', 's*A', 'Coulomb')
_addUnit('V', 'W/A', 'Volt')
_addUnit('F', 'C/V', 'Farad')
_addUnit('ohm', 'V/A', 'Ohm')
_addUnit('S', 'A/V', 'Siemens')
_addUnit('Wb', 'V*s', 'Weber')
_addUnit('T', 'Wb/m**2', 'Tesla')
_addUnit('H', 'Wb/A', 'Henry')
_addUnit('lm', 'cd*sr', 'Lumen')
_addUnit('lx', 'lm/m**2', 'Lux')
_addUnit('Bq', '1/s', 'Becquerel')
_addUnit('Gy', 'J/kg', 'Gray')
_addUnit('Sv', 'J/kg', 'Sievert')
_addUnit('kat', 'mol/s', 'Katal')

_addUnit('abA', '10*A', 'Abampere')

del _unit_table['kg']

for unit in _unit_table.keys():
    _addPrefixed(unit)


# Fundamental constants, as far as needed to define other units
_unit_table['pi'] = np.pi
_addUnit('c0', '299792458.*m/s', 'speed of light')
_addUnit('mu0', '4.e-7*pi*N/A**2', 'permeability of vacuum')
_addUnit('eps0', '1/mu0/c0**2', 'permittivity of vacuum')
_addUnit('hplanck', '6.62606957e-34*J*s', 'Planck constant')
_addUnit('hbar', 'hplanck/(2*pi)', 'Planck constant / 2pi')
_addUnit('e0', '1.602176565e-19*C', 'elementary charge')
_addUnit('me', '9.10938291e-31*kg', 'electron mass')
_addUnit('kb', '1.3806488e-23*J/K', 'Boltzmann constant')

# Time units
_addUnit('min', '60*s', 'minute')
_addUnit('h', '60*min', 'hour')
_addUnit('d', '24*h', 'day')
_addUnit('wk', '7*d', 'week')
_addUnit('yr', '365.25*d', 'year')
_addPrefixed('yr')
_addUnit('fortnight', '1209600*s', '14 days')

# Length units
_addUnit('inch', '2.54*cm', 'inch')
_addUnit('ft', '12*inch', 'foot')
_addUnit('yd', '3*ft', 'yard')
_addUnit('mi', '5280.*ft', '(British) mile')
_addUnit('nmi', '1852.*m', 'Nautical mile')
_addUnit('Ang', '1.e-10*m', 'Angstrom')
_addUnit('AA', '1.e-10*m', 'Angstrom')
_addUnit('lyr', 'c0*yr', 'light year')
_addUnit('Bohr', '4*pi*eps0*hbar**2/me/e0**2', 'Bohr radius')
_addUnit('furlong', '201.168*m', 'furlongs')
_addUnit('au', '149597870691*m', 'astronomical unit')

# Area units
_addUnit('ha', '10000*m**2', 'hectare')
_addUnit('acres', 'mi**2/640', 'acre')
_addUnit('b', '1.e-28*m', 'barn')

# Volume units
_addUnit('l', 'dm**3', 'liter')
_addUnit('dl', '0.1*l', 'deci liter')
_addUnit('cl', '0.01*l', 'centi liter')
_addUnit('ml', '0.001*l', 'milli liter')
_addUnit('mul', '0.000001*l', 'micro liter')
_addUnit('tsp', '4.92892159375*ml', 'teaspoon')
_addUnit('tbsp', '3*tsp', 'tablespoon')
_addUnit('floz', '2*tbsp', 'fluid ounce')
_addUnit('cup', '8*floz', 'cup')
_addUnit('pt', '16*floz', 'pint')
_addUnit('qt', '2*pt', 'quart')
_addUnit('galUS', '4*qt', 'US gallon')
_addUnit('galUK', '4.54609*l', 'British gallon')

# Mass units
_addUnit('t', '1000*kg', 'Metric ton')
_addUnit('amu', '1.660538921e-27*kg', 'atomic mass units')
_addUnit('Da', '1*amu', 'Dalton')
_addUnit('oz', '28.349523125*g', 'ounce')
_addUnit('lb', '16*oz', 'pound')
_addUnit('ton', '2000*lb', 'US ton')

# Force units
_addUnit('dyn', '1.e-5*N', 'dyne (cgs unit)')

# Energy units
_addUnit('erg', '1.e-7*J', 'erg (cgs unit)')
_addUnit('eV', 'e0*V', 'electron volt')
_addUnit('Hartree', 'me*e0**4/16/pi**2/eps0**2/hbar**2', 'Wavenumbers/inverse cm')
_addUnit('Ken', 'kb*K', 'Kelvin as energy unit')
_addUnit('cal', '4.184*J', 'thermochemical calorie')
_addUnit('kcal', '1000*cal', 'thermochemical kilocalorie')
_addUnit('cali', '4.1868*J', 'international calorie')
_addUnit('kcali', '1000*cali', 'international kilocalorie')
_addUnit('Btu', '1055.05585262*J', 'British thermal unit')

_addPrefixed('eV')

# Electromagnetic units
_addUnit('G', '1e-4*T', 'Gauss')
_addUnit('Oe', '79.5774715*A/m', 'Oersted')

_addPrefixed('G')
_addPrefixed('Oe')

# Power units
_addUnit('hp', '745.7*W', 'horsepower')

# Pressure units
_addUnit('bar', '1.e5*Pa', 'bar (cgs unit)')
_addUnit('mbar', '1.e2*Pa', 'millibar')
_addUnit('kbar', '1.e8*Pa', 'kilobar')
_addUnit('atm', '101325.*Pa', 'standard atmosphere')
_addUnit('torr', 'atm/760', 'torr = mm of mercury')
_addUnit('psi', '6894.75729317*Pa', 'pounds per square inch')

# Angle units
_addUnit('deg', 'pi*rad/180', 'degrees')
_addUnit('arcmin', 'pi*rad/180/60', 'minutes of arc')
_addUnit('arcsec', 'pi*rad/180/3600', 'seconds of arc')
_unit_table['cycles'] = 2*np.pi

# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
kelvin = _findUnit('K')
_addUnit('degR', '(5./9.)*K', 'degrees Rankine')
_addUnit('degC', PhysicalUnit(None, 1.0, kelvin.powers, 273.15),
         'degrees Celcius')
_addUnit('degF', PhysicalUnit(None, 5./9., kelvin.powers, 459.67),
         'degree Fahrenheit')
del kelvin

# Radiation-related units
_addUnit('Ci', '3.7e10*Bq', 'Curie')
_addUnit('rem', '0.01*Sv', 'Rem')

_addPrefixed('Ci')
_addPrefixed('rem')

# Astronomical units
_addUnit('Msol', '1.98892e30*kg', 'solar mass')
_addUnit('Lsol', '3.839e26*W', 'solar luminosity')
_addUnit('pc', '3.08568025e16*m')
_addPrefixed('pc')
