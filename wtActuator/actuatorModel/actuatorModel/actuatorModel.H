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
    Foam::actuatorModel

Description
    Base class for actuator models.

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::fv::wtActuatorSource

SourceFiles
    actuatorModel.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModel_H
#define actuatorModel_H

#include "wtActuatorSource.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class actuatorModel Declaration
\*---------------------------------------------------------------------------*/

// TODO: templatize this class with RhoFieldType for compressible cases
// template <class RhoFieldType>
class actuatorModel
{
public:
    // Public Enumerations

        //- Options for the tip factor specification
        enum tipFactorName
        {
            TIP_NONE,
            TIP_SHEN,
            TIP_GLAUERT,
            TIP_PRANDTL
        };

        //- Names for tip factor
        static const Enum<tipFactorName> tipFactorName_;

        //- Options for the root factor specification
        enum rootFactorName
        {
            ROOT_NONE,
            ROOT_GLAUERT,
            ROOT_SORENSEN
        };

        //- Names for root factor
        static const Enum<rootFactorName> rootFactorName_;

protected:

    // Protected Data

        //- Reference to the rotor source model
        const fv::wtActuatorSource& rotor_;

        //- Name of model
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;

        //- tipFactor type
        tipFactorName tipFactorSelect;

        //- rootFactor type
        rootFactorName rootFactorSelect;

        //- root factor characteristic distance
        scalar rDist_;

        bool gradInterpolation_ = false;
        int saveLevel_;
        bool saveNodeForces_;

    // Protected Functions
        scalar (Foam::actuatorModel::*tipfactor_)(const scalar rNode, const scalar tsr, const scalar phi);
        scalar tipfactor_Shen(const scalar rNode, const scalar tsr, const scalar phi);
        scalar tipfactor_Glauert(const scalar rNode, const scalar tsr, const scalar phi);
        scalar tipfactor_Prandtl(const scalar rNode, const scalar tsr, const scalar phi);
        scalar tipfactor_none(const scalar rNode, const scalar tsr, const scalar phi);

        scalar (Foam::actuatorModel::*rootfactor_)(const scalar rNode, const scalar phi, const scalar rDist);
        scalar rootfactor_Glauert(const scalar rNode, const scalar phi, const scalar rDist);
        scalar rootfactor_Sorensen(const scalar rNode, const scalar phi, const scalar rDist);
        scalar rootfactor_none(const scalar rNode, const scalar phi, const scalar rDist);

        void ntrVectors(vector &coordNode, vector &uniBladeDir, vector &uniThetaDir, tensor &iTransform);
        void nodeUrel(vector &Unode_ntr, scalar &rNode, vector &Urel, scalar &phi);

        void posInList(const List<scalar> &inputList, const scalar &searchValue, label &pos1, label &pos2);

public:

    //- Run-time type information
        TypeName("actuatorModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            actuatorModel,
            dictionary,
            (
                const fv::wtActuatorSource& rotor,
                const dictionary& dict
            ),
            (rotor, dict)
        );

    // Constructors

        //- Construct from components
        actuatorModel
        (
            const fv::wtActuatorSource& rotor,
            const dictionary& dict,
            const word& name
        );


    // Selectors

        //- Return a reference to the selected trim model
        static autoPtr<actuatorModel> New
        (
            const fv::wtActuatorSource& rotor,
            const dictionary& dict
        );


    //- Destructor
        virtual ~actuatorModel() = default;


    // Member Functions

        //- Read
        virtual void read(const dictionary& dict);

        virtual void init() = 0;

        virtual void applyForce
        (
            const volVectorField &U,
            // const RhoFieldType &rho, // TODO: use when templatized
            const geometricOneField &rho, // TODO: remove when templatized
            vectorField &Usource,
            scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
        ) = 0;

        // TEMPLATE EXAMPLES
        // //- Correct the model
        // virtual void correct
        // (
        //     const vectorField& U,
        //     vectorField& force
        // ) = 0;

        // //- Correct the model for compressible flow
        // virtual void correct
        // (
        //     const volScalarField rho,
        //     const vectorField& U,
        //     vectorField& force
        // ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
