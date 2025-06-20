/*----------------------------------------------------------------------------*
|   wtActuatorFOAM:                                                           |
|                                                                             |
|    Yet another actuator library to simulate wind turbines in OpenFOAM       |
|                                                                             |
| Copyright (C) 2025 Computational Simulation Center (CSC-CONICET)--Argentina |
|                                                                             |
| This library is based on OpenFOAM: The Open Source CFD Toolbox,             |
| Copyright (C) 2011-2016 OpenFOAM Foundation,                                |
| Copyright (C) 2018-2025 OpenCFD Ltd.                                        |
\*----------------------------------------------------------------------------*

License

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library.  If not, see <https://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------

Class
    Foam::actuatorModel_airfoil

Description
    Derived class computing the nodal forces according to the blade element
    method (BEM).
    (needs detailed airfoils information, may be corrected by tip/root
    factors)

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::fv::wtActuatorSource
  - Foam::actuatorModel

SourceFiles
    actuatorModel_airfoil.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModel_airfoil_H
#define actuatorModel_airfoil_H

#include "actuatorModel.H"

#include "bladeModel.H"
#include "profileModelList.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class actuatorModel_airfoil Declaration
\*---------------------------------------------------------------------------*/

// template <class RhoFieldType>
class actuatorModel_airfoil
:
    public actuatorModel
{
protected:
    // Protected Data
        // - Blade data
        bladeModel *blade_;

        //- Profile data
        profileModelList *profiles_;

        void airfoilAnglesCoeffs(scalar &rNode, scalar &phi, scalar &chord, scalar &twist,
                                 scalar &alpha, scalar &beta, scalar &Cl, scalar &Cd
        );

public:

    //- Run-time type information
    TypeName("airfoil");


    // Constructors

        //- Constructor from rotor and dictionary
        actuatorModel_airfoil(const fv::wtActuatorSource& rotor, const dictionary& dict);

        //- No copy construct
        actuatorModel_airfoil(const actuatorModel_airfoil&) = delete;

        //- No copy assignment
        void operator=(const actuatorModel_airfoil&) = delete;


    //- Destructor
    virtual ~actuatorModel_airfoil() = default;


    // Member Functions
        //- Read
        void read(const dictionary& dict);

        virtual void init();

        virtual void applyForce
        (
            const volVectorField &U,
            // const RhoFieldType &rho, // TODO: use when templatized
            const geometricOneField &rho, // TODO: remove when templatized
            vectorField &Usource,
            scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
