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
    Foam::fv::wtActuatorSource

Group
    grpFvOptionsSources

Description
    Applies cell-based momentum sources on velocity (i.e. \c U)
    within a specified cellSet to approximate the effects of rotor forces.

    Sources applied to:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

Usage
    Refer to the fvOptions file of the cases in the tutorials directory.

See also
  - Foam::actuatorModel
  - Foam::bladeModel
  - Foam::lookupProfile

SourceFiles
    wtActuatorSource.C
    wtActuatorSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef wtActuatorSource_H
#define wtActuatorSource_H

#include "cellSetOption.H"
#include "Enum.H"
#include <fstream>
#include <ostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations

// template <class RhoFieldType>
class actuatorModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class wtActuatorSource Declaration
\*---------------------------------------------------------------------------*/

class wtActuatorSource
:
    public fv::cellSetOption
{
public:

    // // Public Enumerations
    //- Options for the actuator mesh
    enum meshType
    {
        MESH_RINGS
        // MESH_LINES
    };

    //- Names for actuator mesh
    static const Enum<meshType> meshType_;

    //- Options for induction calculation
    enum inductionType
    {
        INDUCTION_SORENSEN,
        INDUCTION_ANALYTIC,
        INDUCTION_TABLE
    };

    //- Names for actuator mesh
    static const Enum<inductionType> inductionType_;

protected:
    // Protected Data
        //- actuatorModel type
        autoPtr<actuatorModel> actuatorModel_;

        //- mesh type
        meshType meshTypeSelect;

        //- induction type
        inductionType inductionTypeSelect;

        //- Disk area
        scalar diskArea_;

        //- Number of turbine blades
        scalar bladesNumber_;

        //- Disk center point sample
        point diskPoint_;

        //- Disk area normal - in flow sense
        vector diskDir_;

        //- disk cell ID
        label diskCellId_;

        // Yaw angle with respect to x-axis (yaw == 360 means self orientation)
        scalar yaw_;

        //- Wind turbine reference velocity
        scalar Uref_;

        //- Power from table interpolation
        scalar Power_;

        //- Power coefficient from table interpolation
        scalar Cp_;

        //- Thrust from table interpolation
        scalar Thrust_;

        //- Thrust coefficient from table interpolation
        scalar Ct_;

        //- Torque from table interpolation
        scalar Torque_;

        //- Rotational velocity [rad/s]
        scalar omega_;

        //- Pitch angle [deg]
        scalar pitch_;

        //- Tip speed ratio
        scalar lambda_;

        //- Reference density [kg/m^3] to compute loads in incompressible solvers
        scalar density_;

        //- cell size in the disc
        scalar cellSize_;
        scalar minSep_;

        //- Parameters
        scalar orientRadiusFrac_;       // Diameter of sphere (wrt maxR_) where flow direction is computed
        scalar selectRadialBoundFrac_;  // Diameter of sphere (wrt maxR_) where actuator cells are selected
        scalar selectAxialBoundFrac_;   // Distance to actuator plane (wrt minSep_) that limits actuator cells selected
        bool unifWeight_;               // Uniform weight on every diskCells_ cell [true, false]
        bool axialWeight_;              // Weighting with distance to actuator plane (gaussian) [true, false]
        scalar sigmaAxialWeightFrac_;   // Width of axialWeight gussian (wrt minSep_)
        bool sphereWeight_;             // Weighting with distance to actuator center (gaussian) [true, false]
        scalar maxRweightFrac_;         // Max radial distance where cells are considered for weighting (wrt maxR_)
        scalar sigmaSphereWeightFrac_;  // Width of sphereWeight Gaussian (wrt maxR_)
        scalar sphereTolDiskFrac_;      // Max radial distance where cells are considered for force smearing wrt maxR_
        vector distribGaussianCoeffs_;  // Width of Gaussians for force smearing wrt minSep_ (n,t,r)

        // Has the AD moved since the last time step (not yet implemented)
        bool movingAD_ = false;

        // First iteration
        bool firstStep_ = true;

        // Actuator radius
        scalar maxR_;

        // Reduced set of cells near to actuator
        DynamicList<label> diskCells_;

        // Total nodes in the AD
        label nodesNumber_;

        // List of vectors with cylindrical coordinates and area of nodes (r, theta, area)
        DynamicList<vector> nodesList_;

        // List of node radial positions
        DynamicList<scalar> rNodeList_;

        // List of vectors with Cartesian  coordinates of nodes (x, y, z)
        DynamicList<vector> nodesPosList_;

        // List of cell ID for each node
        DynamicList<label> nodesCellId_;

        // Weights for Ud averaging and total weighted volume
        std::map<label, float> weightADcells_;
        scalar weightedADvol_;

        // direction of the disc
        vector diskYawed_;
        vector uniDiskDir_;

        // Actuator configuration
        bool UdCorrection_; // Use Ud correction from Shapiro etal (2019) doi: 10.1002/we.2376.
        bool calibration_;  // Use the gradient to extrapolate Ud to node location

        //- Output level
        int saveLevel_;
        bool saveNodeForces_;

        //- Wind turbine Uref curve list
        List<scalar> wtUrefList_;

        //- Wind turbine Ud curve list
        List<scalar> wtUdList_;

        //- Wind turbine Cp curve list
        List<scalar> wtCpList_;

        //- Wind turbine Ct curve list
        List<scalar> wtCtList_;

        //- Wind turbine omega curve list
        List<scalar> wtOmegaList_;

        //- Wind turbine pìtch curve list
        List<scalar> wtPitchList_;

        // Is the topoSet decomposed?
        bool parted_;

        // Has this processor cells of the topoSet?
        bool hasTopoCells_;

    // Protected Member Functions
        //- Check data
        void checkData() const;

        void checkIfParted();

        void contructWTLists();

        void computeMinSep();

        // Nodes list construction
        void constructNodesList_rings();

        // Nodes coordinates list construction
        void constructNodesPosList();

        vector computeNodePosition(
            scalar rNode,
            scalar titaNode_Rad
        );

        void computeDiskOrientation(const vectorField &U);

        void selectADcells();

        void computeUdWeights();

        vector computeUd(const vectorField &U);

        void correctUd(
            scalar &UdDir,
            vector &Ud,
            scalar epsilon
        );

        scalar induction1(scalar UdDir);
        scalar induction2(scalar UdDir);

        // If UdDir is provided interpolation in terms of Ud, else in terms of Uref_
        void interpolateWTCurves(scalar UdDir);

        void computeWTmagnitudes();

        //- Add force to the UEqn
        template <class RhoFieldType>
        void addActuatorForce
        (
            const volVectorField &U,
            const RhoFieldType &rho,
            vectorField &Usource
        );

        //- actuators output file
        std::ofstream *outActuators = NULL;

        //- actuators complementary output file
        std::ofstream *outActuators2 = NULL;

        //- blades output file
        std::ofstream *outNodes = NULL;

public:

    //- Runtime type information
    TypeName("wtActuator");


    // Constructors
        //- Construct from components
        wtActuatorSource
        (
            const word &name,
            const word &modelType,
            const dictionary &dict,
            const fvMesh &mesh
        );

        //- No copy construct
        wtActuatorSource(const wtActuatorSource&) = delete;

        //- No copy assignment
        void operator=(const wtActuatorSource&) = delete;


    //- Destructor
    ~wtActuatorSource();

    // Member Functions
        void distributeActuatorForces(
            vectorField &Usource,
            vector Faero, // in (x, y, z) coordinates
            label node,
            tensor iTransform
        ) const;

        vector getNodeVelocity(
            label node,
            const volVectorField &U,
            const volTensorField *gradU
        ) const;

        void saveNodeForces(
            label node,
            vector Unode,
            scalar Faero_n,
            scalar Faero_t
        ) const;

        // Access
            //- Return rotor radius
            inline scalar maxR() const;

            //- Return rotor area
            inline scalar diskArea() const;

            //- Return rotor thrust
            inline scalar Thrust() const;

            //- Return rotor thrust coeff
            inline scalar Ct() const;

            //- Return rotor Torque
            inline scalar Torque() const;

            //- Return rotor Uref
            inline scalar Uref() const;

            //- Return fvOptions coeff dictionary
            inline dictionary getCoeffsDict() const;

            //- Return the number of actuator nodes
            inline scalar nodesNumber() const;

            //- Return the number of actuator nodes
            inline scalar bladesNumber() const;

            //- Return the
            inline DynamicList<vector> nodesList() const;

            //- Return the
            inline DynamicList<vector> nodesPosList() const;

            //- Return the list por node radial positions
            inline DynamicList<scalar> rNodeList() const;

            //- Return the
            inline vector uniDiskDir() const;

            //- Return the
            inline point diskPoint() const;

            //- Return the rotational speed [rad/s]
            //  Positive anti-clockwise when looking along the lift direction
            inline scalar omega() const;

            //- Return the pitch angle
            inline scalar pitch() const;

            //- Return the pitch angle
            inline scalar lambda() const;

            //- Return the pitch angle
            inline scalar density() const;

        // Source term addition
            //- Add explicit contribution to momentum equation
            virtual void addSup(
                fvMatrix<vector> &eqn,
                const label fieldI
            );

            // TODO: make class actuatorModel templated on RhoFieldType
            // //- Add explicit contribution to compressible momentum equation
            // virtual void addSup(
            //     const volScalarField &rho,
            //     fvMatrix<vector> &eqn,
            //     const label fieldI);

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary &dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wtActuatorSourceI.H"

#endif

// ************************************************************************* //
