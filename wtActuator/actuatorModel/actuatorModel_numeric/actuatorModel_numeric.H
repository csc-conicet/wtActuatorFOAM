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
    Foam::actuatorModel_numeric

Description
    Derived class computing the nodal forces according to the actuator disk
    model from from van der Laan etal (2015) doi: 10.1002/we.1816
    (needs forcesRadialDistribution table with forces at each radial position
    for each wind speed, tip/root factors do not apply)

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::fv::wtActuatorSource
  - Foam::actuatorModel

SourceFiles
    actuatorModel_numeric.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModel_numeric_H
#define actuatorModel_numeric_H

#include "actuatorModel.H"

#include "bladeModel.H"
#include "profileModelList.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class actuatorModel_numeric Declaration
\*---------------------------------------------------------------------------*/

// template <class RhoFieldType>
class actuatorModel_numeric
:
    public actuatorModel
{
protected:
    // Protected Data
        //- list for r positions in the original table
        List<scalar> posrList_;

        //- Fn values along r list
        List<scalar> fnList_;

        //- Ft values along r list
        List<scalar> ftList_;

        //- list of unique Uinf values
        DynamicList<scalar> UinfOnlyList_;
        label nUinf_orig_; //- Size of UinfOnlyList_

        //- list of unique Uref values
        DynamicList<scalar> UrefOnlyList_;
        label nUref_orig_; //- Size of UrefOnlyList_

        //- list of unique node radial positions
        DynamicList<scalar> rNodeList_;
        label nR_; //- Size of rNodeList_

        //- list of radial positions of actuator nodes
        DynamicList<scalar> rNodePosList_;

        //- list of radial positions of each Uref value in UinfOnlyList_
        DynamicList<scalar> posUrefInUinf_;

        //- Total Thrust and Torque from original forces distribution (used to normalize nodal forces)
        DynamicList<scalar> sumFaero_t_;
        DynamicList<scalar> sumFaero_n_;

public:

    //- Run-time type information
    TypeName("numeric");


    // Constructors

        //- Constructor from rotor and dictionary
        actuatorModel_numeric(const fv::wtActuatorSource& rotor, const dictionary& dict);

        //- No copy construct
        actuatorModel_numeric(const actuatorModel_numeric&) = delete;

        //- No copy assignment
        void operator=(const actuatorModel_numeric&) = delete;


    //- Destructor
    virtual ~actuatorModel_numeric() = default;


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
