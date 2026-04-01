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
    Foam::actuatorModel_genanalytic

Description
    Derived class computing the nodal forces according to the actuator disk
    model from Sørensen (2020) doi: 10.1016/j.renene.2019.09.134
    (may use tip/root correction factors)

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::fv::wtActuatorSource
  - Foam::actuatorModel

SourceFiles
    actuatorModel_genanalytic.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModel_genanalytic_H
#define actuatorModel_genanalytic_H

#include "actuatorModel.H"

#include "bladeModel.H"
#include "profileModelList.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class actuatorModel_genanalytic Declaration
\*---------------------------------------------------------------------------*/

// template <class RhoFieldType>
class actuatorModel_genanalytic
:
    public actuatorModel
{
public:

    //- Run-time type information
    TypeName("genanalytic");


    // Constructors

        //- Constructor from rotor and dictionary
        actuatorModel_genanalytic(const fv::wtActuatorSource& rotor, const dictionary& dict);

        //- No copy construct
        actuatorModel_genanalytic(const actuatorModel_genanalytic&) = delete;

        //- No copy assignment
        void operator=(const actuatorModel_genanalytic&) = delete;


    //- Destructor
    virtual ~actuatorModel_genanalytic() = default;


    // Member Functions
        //- Read
        void read(const dictionary& dict);

        virtual void init();

        scalar induction1(scalar &UdDir);

        scalar induction_1(scalar Uref);

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
