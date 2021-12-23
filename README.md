# Black-Body Radiation UNDER CONSTRUCTION

In this project I attempt to reproduce the Planck Black Body spectrum using a computer simulation.

## Introduction

An object in thermodynmaic equilibrium with its environment emits electromagnetic 
radiation; the spectrum of this raiation is defined by the temperature of the object.
This radiation is called Black Body Radiation.

The higher the object's temperature, the higher the frequency of the Black-Body radiation.

| TEMPERATURE  | DESCRIPTION                  |  FREQ-RANGE-THz  |  WAVELENGTH-RANGE-NM |
| -----------  | -----------                  |  --------------  |  ------------------- |
| 2.7          | Cosmic Microwave Background  |   .1  -    .7    |  500000  -  3400000  |
| 300          | Room Temperature             |   10  -    80    |    4000  -    30000  |
| 5800         | Surface of Sun               |  200  -  1400    |     210  -     1500  |

The wavelength of visible light ranges from 380 nm (violet) to 700 nm (red).

![plot-planck-5800.png](/assets/plot_planck_5800.png)

## Brief History

```
1687  Newton's Laws of Motion
1689  Leibniz, vis viva, Kinetic Energy
1801  Young Double Slit Experiment
1834  Ideal Gas Law
1843  Equipartition of kinetic energy first proposed
1861  Maxwell's Equations
1900  Planck's Law , empirically derived, uses quantized energy
1900  Rayleigh Jeans law, based on thermodynamics
1905  Einstein conceived of quanta of light, Photoelectric Effect
1918  Max Planck awarded the Nobel Prize in Physics for his discovery of the enegy quanta
1921  Albert Einstein awarded the Nobel Prize in Physics for the Photoelectirc Effect
1926  The name 'Photon' is used to describe Einstein's quanta of light
```

In the late 1600s, Newton and Leibniz discovered the laws of motion and conservation
of kinetic energy. These discoveries provided the ability to predict the result of an
elastic collision of 2 objects, given that you know the objects masses and initial velocities.

In 1801 Thomas Young discovered that light is a wave, using a technique similar to the modern
double slit experiment. Prior to this, light was believed to be made of tiny particles called
corpuscles. In 1818 the wave theory of light replaced the corpuscle theory when an experiment
observed the Poisson Spot, which could only be explained if light is a wave. A century later,
Einstein discovered that light is both a particle and a wave.

In 1834 the Ideal Gas Law was conceived. This law assumes that a gas is made up of many tiny
particles that collide elastically. Newton/Leibniz could predict the outcome of individual
collisions. The Ideal Gas Law instead predicts the average behaviour of the entire volume 
of gas. The modelling of large number of particles is known as Statistical Mechanics.

In 1843 Equipartition of kinetic energy first proposed. In statistical mechanics this 
relates system temperature to the average kinetic energy for each degree of freedom. For example, 
consider many pucks on an Air Hockey Table, where the sides of the table vibrate to impart
energy to the pucks. In a short time, the pucks will reach an equilibrium state. The pucks
have 2 degrees of freedom, they can move in either the X or Y direction. The average energy
of the pucks due to their motion in the X direction will be equal to the average energy of 
the pucks in the Y direction. For a gas of molecules, the degrees of freedom are: motion in
the X,Y,Z direction; rotation about the X,Y,Z axis; and longitudinal osciallations of the
atoms that make up the moulecules. According to the Equipartition theorem each of these degrees 
of freedom will have the same average energy. The average energy for each degree of freedom 
(aka mode) is the Boltzmann constant times the temperature of the gas (Kb * T). Gasses that 
are comprised of molecules that have more degrees of freedom will have a higher capacity to store 
thermal energy.

In 1861 Maxwell built on the work of Gauss, Faraday and Ampere, and published Maxwell's Equations.
Maxwell's equations show that light waves propogate due to an oscillation between electric and 
mangetic fields. The speed of light propogation in vacuum is related to the Vacuum Permittivity constants.
These constants define the ability of a vacuum to hold an electric and a magnetic field.

In 1900 Maxwell Planck discovered a formula for the black body spectrum. At the time thee were
accurate experimental measurements of the Black Body spectrum. Planck used empirical techniques
to derive the formula, and Planck did not recognize the significance of his accomplishment at
the time. Planck was later awarded the Nobel Prize in Phyiscs for his discovery of the energy
quanta. Planck's Black Body spectrum formula (Planck's Law) was derived by assuming that the 
oscillators (that create the light waves) can have energy values only in minimal increments.
This is the quantization of energy. The Planck formula for the black body spectrum defines a new
physical constant 'h' (the Planck constant), which defines the value of the minimal energy 
increment. Planck's Law was developed without making use of the Rayleigh-Jeans Law.
The Planck Constant 'h', is now widely used in quantum mechanics.

In 1900  The Rayleigh-Jeans law made use of the Equipartion therom to predict the Black Body
spectrum. Their approach starts by assuming that light waves are in thermodynamic equilibrium with
a chamber containing the waves. Rayleigh and Jeans calculated the number of modes of light
waves that could be contained in a chamber, vs the frequency of the light waves. According to the
Equipartion Therom, each mode would have Kb * T energy. The Rayleigh-Jeans Law for black body
spectral radiance says that the Spectral Radiance at frequency 'f' is Kb * T *  the number of
light wave modes at that frequency. This equation yields results that are very close to 
experiment and Planck for low frequencies; however the equation is way off at high frequencies, and
goes to infinity as the frequency increases. This problem was called the Ultraviolet Catastrophe.
The Ultraviolet Catastrophe was later resolved by taking into account the quantization of energy.

In 1905 Albert Einstein worked on the photolectric effect, which led him to conceive of the
quanta of light.  Einstein was later awared the Nobel Prize in Physics for this. The quanta
of light discovered by Einstein in the Photoelectric effect was equivalent to the quanta 
of energy that Planck had proposed. However Einstein showed that it is light that is quantized; 
Planck had earlier thought it was the oscaillator energy that was quantized.

In 1926, the name 'Photon' is first used to describe the Einstein quanta of light.

## Black Body Spectrum

An object in thermodynmaic equilibrium with its environment emits electromagnetic 
radiation; the spectrum of this raiation is defined by the temperature of the object.
This radiation is called Black Body Radiation.

Formulas for the black body spectrum come in two forms, where the spectral radiance is either:
probdist_get_value
- power emitted per unit emitting area, per steradian, per unit wavelength
- power emitted per unit emitting area, per steradian, per unit frequency

I use the formula which returns power per unit wavelength, because this formula generates a 
spectrum that shows the peak intensity for T=5800 (the surface of the sun) to be at center of 
the visible light spectrum.

## Planck's Law

Reference: https://en.wikipedia.org/wiki/Planck%27s_law
```
                       2 h c^2          1
    SpectralRadiance = ------- *  -------------------
                       wvlen^5           h c 
                                  e^ (----------)  - 1
                                      wvlen * KT
```

## Rayleigh-Jeans Law

Reference: https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Jeans_law
```
                        2 c KT
    SpectralRadiance = ---------
                        wvlen^4
```

The bb.c calc_rj() function:
```
                             8 * π       c
    MODE_DENSITY(wvlen) =  --------- * --------
                           wvlen^4      4 * π

    SpectralRadiance = MODE_DENSITY(wvlen) * KT

    The first term in MODE_DENSITY is taken from the mode density derivation provided here
    http://www.reading.ac.uk/physicsnet/units/3/3pha4/Lectures/l1.pdf, 
    and allowing for 2 polarizations.

    The second term in MODE_DENSITY is a units conversion factor, so that the
    calc_rj() function in bb.c returns the correct value.
```

The spectral radiance calculated using the Rayleigh-Jeans law agrees with experiment and Planck 
at low frequencies; however, at high frequencies Rayleigh-Jeans diverges dramatically from Planck. 
RJ goes to infinity as wavelength goes to zero. This was called the Ultraviolet Catastrophe.

![plot-rj-vs-planck.png](/assets/plot_rj_vs_planck.png)

## Solving the Ultraviolet Catastrophe

The RJ Law's Ultraviolet Catastrophe is resolved by incorporating energy quantization.

Recall that the mode density is a function of wavelength (or frequency), and each mode
has a constant energy content equal to KT.

For example, assuming KT=1000 ...
* Low Frequency, hf is much smaller than KT: Rounding down to a multiple of hf makes just a small change. For example, KT=1000, and a random energy value conforming to this KT could be 977.5. Assuming the low frequency hf value is 2; then the 977.5 is rounded down (quantized) to a multiple of 2. It's value changes from 977.5 to 976.
* High Frequency, hf is close to KT: Rounding down to a multiple of hf makes a large change. For example, KT=1000, and a random energy value conforming to this KT could be 977.5. Assuming the high frequency hf value is 500; then the 977.5 is rounded down (quantized) to a multiple of 500. It's value changes from 977.5 to 500.

## My Black Body Calculation

The bb program adds energy quantization to the Rayleigh-Jeans Law.

The spectral radiance for a given wavlength is computed as follows:
- an energy value conforming to the maxwell-boltzmann energy distribution is obtained many times,
- each time this energy value is quantized by rounding down to a multiple of hf
- the average of these quantized energy values is calculated
- this average quantized energy value is then multiplied by the MODE_DENSITY to obtain the Spectral Radiance

The plots of Spectral Radiance vs Log Frequency for the Rayleigh-Jeans law, the Planck Law, and 
my calculation are shown below. My calculation is a close match to Planck, however:
- to get this agreement I had to use a value of h that is .777 of the real value of h
- their is a small diffeence between the plots of Mine vs Placnk

I have been unable to account for the .777 or the small difference between the plots of 
Mine vs Planck. 
- The small difference could be the result of the numerical methods used, however
  I have attempted to increase the accuracy, for example by averaging more samples, without 
  improving the result
- I have no idea why the .777 fudge factor for the value of h is needed.

Notes:
- The 'Maxwell-Boltzmann kinetic energy distribution' is from here:
[tec-science Maxwell–Boltzmann distribution](https://www.tec-science.com/thermodynamics/kinetic-theory-of-gases/maxwell-boltzmann-distribution/)
- The Mode Density derivation is from here: 
[reading.ac.uk Radiation in a Cavity (Rayleigh, Planck (1900))](http://www.reading.ac.uk/physicsnet/units/3/3pha4/Lectures/l1.pdf)

![plot-rj-planck-mine.png](/assets/plot_rj_planck_mine.png)

## APPENDIX References

[Wikipedia Planck's law](https://en.wikipedia.org/wiki/Planck%27s_law)

[Wikipedia Rayleigh–Jeans law](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Jeans_law)

[Wikipedia Equipartition theorem](https://en.wikipedia.org/wiki/Equipartition_theorem)

[Wikipedia Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant)

[Wikipedia Maxwell–Boltzmann distribution](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution)

[tec-science Maxwell–Boltzmann distribution](https://www.tec-science.com/thermodynamics/kinetic-theory-of-gases/maxwell-boltzmann-distribution/)

[reading.ac.uk Radiation in a Cavity (Rayleigh, Planck (1900))](http://www.reading.ac.uk/physicsnet/units/3/3pha4/Lectures/l1.pdf)

[Reading Feynman](https://readingfeynman.org/2014/09/22/plancks-constant-h/)

## APPENDIX - MORE REFERENCES

[Writing-On-Github](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)

[Gnuplot doc](http://gnuplot.sourceforge.net/docs_4.2/gnuplot.html)

UTF-8 References
- https://en.wikipedia.org/wiki/UTF-8
- https://en.wikipedia.org/wiki/Block_Elements
- https://www.compart.com/en/unicode/block
- https://www.compart.com/en/unicode/U+03C0

Using vi, insert  ctrl-v u code
```
  example: ^vu03c0 inserts π

```
