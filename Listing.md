# Listing
This notebook is an interactive way to see what data is available to MCNP. It is a modern version of the Listing of Available ACE Data Tables ([LA-UR-17-20709](https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-17-20709)). It parses the XSDIR on your machine and shows the data that is available to you.

If you are concerned about what this does you are free to examine the code yourself. It is available as an Open Source project on GitHub at <https://github.com/NuclearData/DataListing>.

## Types of Nuclear Data
In the table below are given the different kinds of data, generally they are referred to by the incident particle. Each type is a link to a notebook where the details can be explored. In the table below lists the different types of data and the associated `lib_type`, the one- or two-letter extension of the full ZAID. 


### Ordering of Listing Tables
The order of the elements in the listing tables is important; it is the same order as in the XSDIR file. When MCNP looks for data from the material card, it scans the XSDIR file---starting at the top---until it finds the *first* matching ZAID (or partial ZAID). Thus, the order of the XSDIR file is important and it is reflected in the order of the listing tables.

| Type                                                           | `lib_type`  |
|:-------------------------------------------------------------- | --------    |
| [Continuous-Energy Neutron](ContinuousEnergyNeutron.ipynb)     | `c` or `nc` |
| [Thermal Scattering S(⍺,ß) Neutron](ThermalScattering.ipynb)   | `t`         |
| [Discrete-Energy Neutron](DiscreteEnergy.ipynb)                | `d`         |
| [Coupled Discrete-Energy Neutron-Photon](NeutronPhoton.ipynb)  | `m`         |
| [Photoatomic](Photoatomic.ipynb)                               | `p`         |
| [Photonuclear](Photonuclear.ipynb)                             | `u`         |
| [Dosimetry](Dosimetry.ipynb)                                   | `y`         |
| [Electron](Electron.ipynb)                                     | `e`         |
| [Charged-Particle](ChargedParticle.ipynb)                      | varied      |

## Reaction Types and Tallying Reaction Rates
In may of the notebooks linked above, you'll find a table of reaction types after the listing of the data. The reaction types are identified using the ENDF `MT` numbers. These numbers uniquely define a specific reaction. For example, `MT` 102, $(z,\gamma)$ is the "radiative capture" cross section induced by some incident particle $z$. For more information, see [Appendix B of the ENDF manual](https://www.nndc.bnl.gov/csewg/docs/endf-manual.pdf).

The `MT` numbers can be used in MCNP with the tally multiplier card, `FM`. The `FM` card multiplies any tallied quantity (flux, current) by the specified cross section to produce reaction rates. The `MT` numbers are also used to identify reactions when plotting cross sections. In addition to the ENDF `MT` numbers, MCNP also defines `FM` numbers that are similar to the `MT` numbers. **Note:** `MT` numbers are always positive, while `FM` numbers are always negative. In each of the table of reation types, both the `MT` and `FM` numbers are given. 

For more information about tallying reaction rates, see the MCNP manual. 