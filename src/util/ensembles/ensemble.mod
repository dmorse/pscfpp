namespace Util
{

/**
* \defgroup Ensemble_Module Statistical Ensembles
*
* This module contains classes that represent statistical mechanical
* ensembles for particular macroscopic variables. The classes 
* EnergyEnsemble, BoundaryEnsemble, and SpeciesEnsemble represent 
* statistical ensembles for fluctuations of energy, boundary volume 
* or shape, and number of molecules of a single species, respectively. 
*
* Each ensemble has a type, which specifies whether it is constrained
* or Boltzmann ensemble, and stores whatever parameters are needed to
* describe the Boltzmann ensemble.  The type of a EnergyEnsemble can
* be "adiabatic" or "isothermal", and a value of temperature is stored
* if it is isothermal. Similarly, BoundaryEnsemble can be rigid or 
* isobaric, and stores a pressure value if isobaric. Species ensemble 
* can be "closed" or "grand", and stores a value for the chemical 
* potential if it is grand.
*
* We envision generalizing the BoundaryEnsemble to allow for constant
* stress ensembles, and may generalize some or all of the ensembles 
* to allow for non-Boltzmann statistical ensembles. 
*
* \ingroup Util_NS_Module
*/

}

