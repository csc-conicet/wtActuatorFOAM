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
    Foam::actuatorModel_adaptive

Description
    Derived class computing the nodal forces according to the actuator disk
    model from Navarro Diaz etal (2019) doi: 10.1016/j.jweia.2018.12.018
    (needs forcesRadialDistribution table with forces at each radial position
    for each wind speed and turbine wind speed reference, tip/root factors
    do not apply).

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::fv::wtActuatorSource
  - Foam::actuatorModel

SourceFiles
    actuatorModel_adaptive.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModel_adaptive_H
#define actuatorModel_adaptive_H

#include "actuatorModel.H"
#include "bladeModel.H"
#include "profileModelList.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class actuatorModel_adaptive Declaration
\*---------------------------------------------------------------------------*/

// template <class RhoFieldType>
class actuatorModel_adaptive
:
    public actuatorModel
{
protected:
    // Protected Data
        //- Interpolated lists
        //- UdiList
        List<scalar> UdiList_;

        //- fnList
        List<scalar> fnList_;

        //- ftList
        List<scalar> ftList_;

        //- list for r positions in the original table
        List<scalar> posrList_;

        //- list of Uinf values
        DynamicList<scalar> UinfOnlyList_;
        label nUinf_orig_; //- Size of UinfOnlyList_

        //- list of Uref values
        DynamicList<scalar> UrefOnlyList_;
        label nUref_orig_; //- Size of UrefOnlyList_

        //- list of node unique radial positions
        DynamicList<scalar> rNodeList_;
        label nR_; //- Size of rNodeList_

        //- list of radial positions of actuator nodes
        DynamicList<scalar> rNodePosList_;

        void posInTableUdi(const label &posUref, const scalar &posrNode, const scalar &magUnode, label &pos1, label &pos2);

public:

    //- Run-time type information
    TypeName("adaptive");


    // Constructors

        //- Constructor from rotor and dictionary
        actuatorModel_adaptive(const fv::wtActuatorSource& rotor, const dictionary& dict);

        //- No copy construct
        actuatorModel_adaptive(const actuatorModel_adaptive&) = delete;

        //- No copy assignment
        void operator=(const actuatorModel_adaptive&) = delete;


    //- Destructor
    virtual ~actuatorModel_adaptive() = default;


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
