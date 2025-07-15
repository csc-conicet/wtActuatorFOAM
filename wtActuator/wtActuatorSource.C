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

\*---------------------------------------------------------------------------*/

#include "wtActuatorSource.H"
#include "actuatorModel.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "fvCFD.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

#include "wtActuatorSourceTemplates.C"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(wtActuatorSource, 0);
        addToRunTimeSelectionTable(option, wtActuatorSource, dictionary);
    }
}

const Foam::Enum<Foam::fv::wtActuatorSource::meshType>
Foam::fv::wtActuatorSource::meshType_
({
    { wtActuatorSource::MESH_RINGS, "rings" },
    // { wtActuatorSource::MESH_LINES, "lines" },
});

const Foam::Enum<Foam::fv::wtActuatorSource::inductionType>
Foam::fv::wtActuatorSource::inductionType_
({
    { wtActuatorSource::INDUCTION_SORENSEN, "sorensen" },
    { wtActuatorSource::INDUCTION_ANALYTIC, "analytic" },
    { wtActuatorSource::INDUCTION_TABLE, "table" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::wtActuatorSource::computeMinSep()
{
    minSep_ = cellSize_; // minSep_ from user input or default value (1000)

    // In case the user didn't input any fixed cellSize_, it is computed looping over the cells
    if (cellSize_ == 1000)
    {
        forAll(cells(), c1) // Loop on all topoSet cells
        {
            List<label> cellNeighb = mesh().cellCells()[cells()[c1]];

            // Loop on neighbour cells
            forAll(cellNeighb, c2)
            {
                vector arrow = mesh().cellCentres()[cells_[c1]] - mesh().cellCentres()[cellNeighb[c2]];

                scalar vdist = mag(arrow & vector(0, 0, 1));
                scalar hdist = mag(arrow - (arrow & vector(0, 0, 1)) * vector(0, 0, 1));

                if ( hdist > vdist && hdist < minSep_)
                {
                    minSep_ = hdist;
                }
            }
        }

        // Select the minimum minSep_ among all processors
        reduce(minSep_ , minOp<scalar>());
    }
    Info << "Minimum horizontal cell centers separation: " << minSep_ << endl;
}

void Foam::fv::wtActuatorSource::constructNodesList_rings()
{
    Info << "* Rings node distribution" << endl;
    scalar nodesCellsRatio = readScalar(coeffs_.lookup("nodesCellsRatio"));
    scalar rThicknessCellsizeRatio = readScalar(coeffs_.lookup("rThicknessCellsizeRatio"));

    // Estimate the number of nodes in the actuator mesh
    scalar cellsInAd = diskArea_ / pow(minSep_, 2);
    label estimatedNodes = nodesCellsRatio * cellsInAd + 0.5;
    Info << "Estimated number of AD nodes: " << estimatedNodes << endl;

    scalar secArea = diskArea_ / estimatedNodes; // Target area corresponding to each node
    scalar rInt = sqrt(secArea / M_PI); // Inner radius of inner ring
    scalar rExt = maxR_; // Outer radius of outer ring

    // Estimate ringThickness according to cellsize
    scalar ringThickness = rThicknessCellsizeRatio * minSep_;

    label numberRings = (rExt - rInt) / ringThickness + 0.5; // Number of rings
    Info << "Number of rings: " << numberRings << endl;

    // Correct ringThickness according to numberRings
    ringThickness = (rExt - rInt) / numberRings;

    nodesNumber_ = 0;       // Nodes in the actuator mesh (to be accumulated)
    label nodesRing ;      // Nodes in each ring
    scalar dTitaRing;       // Angle between nodes in ring [deg]
    scalar titaNode_Rad;    // Angular position of a node [rad] (couter-clockwise from the positive z as seen from the front)
    scalar areaNodeRing;    // Area of each ring nodes

    // Center node first
    vector nodePolarCoord = vector(0, 0, secArea);
    nodesList_.append(nodePolarCoord);
    nodesNumber_++;
    rNodeList_.append(0); // Set up a list with unique node radial positions

    scalar rMedRing = rInt + ringThickness / 2; // Mean radius of first ring (then updated)
    while (rMedRing  < rExt)
    {
        nodesRing = estimatedNodes * (pow(rMedRing + ringThickness / 2, 2) - pow(rMedRing - ringThickness / 2, 2))
                          / (pow(rExt, 2) - pow(rInt, 2)) + 0.5;

        nodesNumber_ += nodesRing;
        dTitaRing = 360.0 / nodesRing; // deg

        rNodeList_.append(rMedRing);

        areaNodeRing = M_PI * (pow(rMedRing + ringThickness / 2, 2) - pow(rMedRing - ringThickness / 2, 2)) / nodesRing;

        // Fill nodesList_
        for (label nodeIterator = 0; nodeIterator < nodesRing; nodeIterator++)
        {
            titaNode_Rad = 2 * M_PI * (dTitaRing * (nodeIterator)) / 360;
            vector nodePolarCoord = vector(rMedRing, titaNode_Rad, areaNodeRing);
            nodesList_.append(nodePolarCoord); // (R, Theta, Area) list
        }
        rMedRing += ringThickness;
    }
    Info << "Final number of AD nodes: " << nodesNumber_ << endl;

    nodesPosList_.reserve(nodesNumber_);
    nodesCellId_.reserve(nodesNumber_);
}

void Foam::fv::wtActuatorSource::checkIfParted()
{
    parted_ = false;
    hasTopoCells_ = cells().size() > 0;
    scalar partsNo = 0;
    if (hasTopoCells_)
    {
        partsNo = 1;
    }
    reduce(partsNo, sumOp<scalar>());

    if (partsNo > 1)
    {
        parted_ = true;
    }
}

void Foam::fv::wtActuatorSource::contructWTLists()
{
    //- Wind turbine curves (Uref, Ud, Cp, Ct, omega, pitch)
    List<List<scalar>> wtCurvesTable_;
    coeffs_.lookup("wtCurvesTable") >> wtCurvesTable_;

    // Construct individual lists for each variable in wtCurvesTable_
    wtUrefList_.setSize(wtCurvesTable_.size());
    wtUdList_.setSize(wtCurvesTable_.size());
    wtCpList_.setSize(wtCurvesTable_.size());
    wtCtList_.setSize(wtCurvesTable_.size());
    wtOmegaList_.setSize(wtCurvesTable_.size());
    wtPitchList_.setSize(wtCurvesTable_.size());

    forAll(wtCurvesTable_, i)
    {
        wtUrefList_[i] = wtCurvesTable_[i][0];
        wtUdList_[i] = wtCurvesTable_[i][1];
        wtCpList_[i] = wtCurvesTable_[i][2];
        wtCtList_[i] = wtCurvesTable_[i][3];
        wtOmegaList_[i] = (wtCurvesTable_[0].size() > 3) ? wtCurvesTable_[i][4] : 0.0;
        wtPitchList_[i] = (wtCurvesTable_[0].size() > 4) ? wtCurvesTable_[i][5] : 0.0;
    }
}

scalar sumArea(const DynamicList<vector> &coordinates)
{
    scalar nodesAreaSum = 0;
    forAll(coordinates, c)
    {
        nodesAreaSum += coordinates[c][2];
    }
    return nodesAreaSum;
}

void Foam::fv::wtActuatorSource::checkData() const
{
    if (magSqr(diskArea_) <= VSMALL)
    {
        FatalErrorInFunction
            << "diskArea is approximately zero"
            << abort(FatalError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorInFunction
            << "disk direction vector is approximately zero"
            << abort(FatalError);
    }

    if (returnReduce(diskCellId_, maxOp<label>()) == -1)
    {
        // This could not be an error if the actuator center belongs to other proc
        FatalErrorInFunction
            << "disk center location " << diskPoint_ << " not found in mesh"
            << abort(FatalError);
    }

    scalar nodesAreaSum = sumArea(nodesList_);
    if (fabs(nodesAreaSum - diskArea_) > (0.01 * diskArea_)) // difference should be less than 1% of disk area
    {
        FatalErrorInFunction
            << "fabs(nodesAreaSum - diskArea_) > (0.01 * diskArea_)"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::wtActuatorSource::wtActuatorSource(
    const word &name,
    const word &modelType,
    const dictionary &dict,
    const fvMesh &mesh)
    : fv::cellSetOption(name, modelType, dict, mesh),
        actuatorModel_(actuatorModel::New(*this, coeffs_)),
        diskArea_(readScalar(coeffs_.lookup("diskArea"))),
        bladesNumber_(coeffs_.lookupOrDefault("blades", 3)),
        diskPoint_(coeffs_.lookup("diskPoint")),
        diskDir_(coeffs_.lookup("diskDir")),
        yaw_(coeffs_.lookupOrDefault("yaw", 360.0)),
        density_(coeffs_.lookupOrDefault("density", 1.225)),
        cellSize_(coeffs_.lookupOrDefault("cellSize", 1000)),
        calibration_(readBool(coeffs_.lookup("calibration"))),
        UdCorrection_(readBool(coeffs_.lookup("UdCorrection"))),
        saveLevel_(coeffs_.lookupOrDefault("saveLevel", 1)),
        saveNodeForces_(coeffs_.lookupOrDefault("saveNodeForces", false)),
        orientRadiusFrac_(coeffs_.lookupOrDefault("orientationRadiusFrac", 1.0)),
        selectRadialBoundFrac_(coeffs_.lookupOrDefault("selectRadialBoundFrac", 1.15)),
        selectAxialBoundFrac_(coeffs_.lookupOrDefault("selectAxialBoundFrac", 3.0)),
        unifWeight_(readBool(coeffs_.lookup("unifWeight"))),
        axialWeight_(readBool(coeffs_.lookup("axialWeight"))),
        sphereWeight_(readBool(coeffs_.lookup("sphereWeight"))),
        sigmaAxialWeightFrac_(coeffs_.lookupOrDefault("sigmaAxialWeightFrac", 1.0)),
        maxRweightFrac_(coeffs_.lookupOrDefault("maxRweight", 1.0)),
        sigmaSphereWeightFrac_(coeffs_.lookupOrDefault("sigmaSphereWeightFrac", 0.5)),
        sphereTolDiskFrac_(coeffs_.lookupOrDefault("sphereTolDiskFrac", 1.15)),
        distribGaussianCoeffs_(coeffs_.lookup("distribGaussianCoeffs")
)
{
    Info << "    - creating wtActuatorSource: " << name_ << endl;

    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    if (yaw_ == 360) // Self orienting disk - compute uniDiskDir_ in each time step
    {
        movingAD_ = true;
    }

    diskCellId_ = mesh.findCell(diskPoint_);
    maxR_ = sqrt(diskArea_ / M_PI); // actuator radius

    contructWTLists();

    if (calibration_)
    {
        Uref_ = readScalar(coeffs_.lookup("Uref"));
        interpolateWTCurves(); // What if no table?
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
        coeffs_.readIfPresent("omega", omega_);
        coeffs_.readIfPresent("pitch", pitch_);
        lambda_ = omega_ * maxR_ / Uref_;

        computeWTmagnitudes();
    }
    else
    {
        Ct_ = coeffs_.lookupOrDefault("Ct", 8.0/9.0); // Ct is needed for the Ud correction and analytical induction
        Cp_ = coeffs_.lookupOrDefault("Cp", 1.0/3.0); // Cp is needed for the analytical induction
    }

    checkIfParted();

    computeMinSep();

    inductionType_.readIfPresent("actuatorInductionModel", coeffs_, inductionTypeSelect);

    meshType_.readIfPresent("actuatorMesh", coeffs_, meshTypeSelect);

    // Select which actuator mesh (e.g. rings, lines)
    switch (meshTypeSelect)
    {
        // case MESH_LINES: // Not yet implemented
        // {
        //     constructNodesList_lines();
        //     break;
        // }
        case MESH_RINGS:
        default:
        {
            constructNodesList_rings();
        }
    }

    actuatorModel_->init();

    checkData();

    if (Pstream::master()) //In Master node
    {
        // outActuators is opened by each AD and the header is overwritten, so
        // we opened once in write mode to have one header and then in append
        // to write results.
        if (saveLevel_)
        {
            outActuators = new std::ofstream("outActuators.csv");
            (*outActuators) << "----- wtActuator output file -----" << std::endl;
            (*outActuators) << "Actuator name, time [s], Uref [m/s], Ud [m/s], Cp, Ct, omega [rad/s], pitch [deg], "
                            << "Power(Uref, Cp) [W], Thrust(Uref, Ct) [N], Torque [Nm]" << std::endl;
            outActuators->close();
            delete outActuators;

            // Reopen outActuators.csv in append mode
            outActuators = new std::ofstream("outActuators.csv", std::ios::app);
        }
        if (saveLevel_ > 1)
        {
            outActuators2 = new std::ofstream("outActuators_extended.csv");
            (*outActuators2) << "----- wtActuator extended output file -----" << std::endl;
            (*outActuators2) << "Actuator name, time [s], Thrust_actuator [N], Torque_actuator [Nm], "
                             << "Thrust_nodes [N], Torque_nodes [Nm]" << std::endl;
            outActuators2->close();
            delete outActuators2;

            // Reopen outActuators_extended.csv in append mode
            outActuators2 = new std::ofstream("outActuators_extended.csv", std::ios::app);
        }

        if (saveNodeForces_)
        {
            fileName rootDir = "outActuatorsForces";
            if (!isDir(rootDir))
            {
                mkDir(rootDir);
            }

            outNodes = new std::ofstream(rootDir + "/" + name_ + "_nodeForces.csv");
            (*outNodes) << "Actuator name, time [s], node#, r [m], theta [rad], area [m^2], x [m], y [m], z [m], ";
            (*outNodes) << "Unode_x [m/s], Unode_y [m/s], Unode_z [m/s], Faero_n [N/m^2], Faero_t [N/m^2]" << std::endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::wtActuatorSource::addSup(
    fvMatrix<vector> &eqn,
    const label fieldI
)
{
    if (hasTopoCells_ || parted_ || firstStep_ || mesh_.time().writeTime())
    {
        vectorField &Usource = eqn.source();
        const volVectorField &U = eqn.psi();

        addActuatorForce(
            U,
            geometricOneField(),
            Usource);

        firstStep_ = false;
        }
}

// TODO: make class actuatorModel templated on RhoFieldType
// void Foam::fv::wtActuatorSource::addSup(
//     const volScalarField &rho,
//     fvMatrix<vector> &eqn,
//     const label fieldI)
// {
//      vectorField &Usource = eqn.source();
//      const volVectorField &U = eqn.psi();
//
//      if (hasTopoCells_ || parted_ || mesh_.time().writeTime())
//      {
//         addActuatorForce(
//             U,
//             rho,
//             Usource);
//      }
// }

bool Foam::fv::wtActuatorSource::read(const dictionary &dict)
{
    return fv::cellSetOption::read(dict);
}

Foam::fv::wtActuatorSource::~wtActuatorSource()
{
    if (Pstream::master())
    {
        if(outActuators != NULL) {
            outActuators->close();
            delete outActuators;
        }
        if(outActuators2 != NULL) {
            outActuators2->close();
            delete outActuators2;
        }
        if(outNodes != NULL) {
            outNodes->close();
            delete outNodes;
        }
    }
}
// ************************************************************************* //
