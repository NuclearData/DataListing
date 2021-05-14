# Listing
This notebook is an interactive way to see what data is available to MCNP. It is a modern version of the Listing of Available ACE Data Tables ([LA-UR-17-20709](https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-17-20709)). It parses the XSDIR on your machine and shows the data that is available to you.

If you are concerned about what this does you are free to examine the code yourself. It is available as an Open Source project on GitHub at <https://github.com/NuclearData/DataListing>.

## Types of Nuclear Data
In the table below are given the different kinds of data, generally they are referred to by the incident particle. Each type is a link to a notebook where the details can be explored. In the table below lists the different types of data and the associated `lib_type`, the one- or two-letter extension of the full ZAID. 

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
| [Proton](Proton.ipynb)                                         | `h`         |
| [Charged-Particle](ChargedParticle.ipynb)                      | varied      |
