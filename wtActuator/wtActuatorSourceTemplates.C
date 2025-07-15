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
#include "volFields.H"
#include <functional>
#include <math.h>
#include "fvc.H"
#include "fvCFD.H"
#include <map>

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //
void Foam::fv::wtActuatorSource::computeDiskOrientation(
    const vectorField &U
)
{
    if (yaw_ == 360) // Self orienting disk - compute uniDiskDir_ in each time step
    {
        scalar sphereTolSqr = pow(orientRadiusFrac_ * maxR_, 2);

        // Velocity inside sphere to be accumulated for orientation calculation
        vector center_orientUd = vector(0, 0, 0);
        // Volume of the center cells of the sphere
        scalar center_orientV = 0.0;

        forAll(cells(), c)
        {
            if (magSqr(mesh().cellCentres()[cells()[c]] - diskPoint_) < sphereTolSqr)
            {
                scalar cellV = mesh().V()[cells()[c]];
                center_orientV += cellV;
                // average Ud vector for the center cells of the sphere
                center_orientUd += U[cells()[c]] * cellV;
            }
        }

        if (parted_ || firstStep_ || mesh_.time().writeTime())
        {
            reduce(center_orientV, sumOp<scalar>());
            reduce(center_orientUd, sumOp<vector>());
        }

        center_orientUd /= center_orientV;

        // No orientation in vertical direction
        center_orientUd[2] = 0;

        // for self orientation depending on the velocity on the disc
        uniDiskDir_ = center_orientUd / mag(center_orientUd);
    }
    else // Orientation imposed (yaw fixed or externally defined)
    {
        scalar yawRad = yaw_ * M_PI / 180;

        // Rotate the orginal diskDir with the yaw angle
        uniDiskDir_ = vector(diskDir_[0]*cos(yawRad) - diskDir_[1]*sin(yawRad),
                             diskDir_[0]*sin(yawRad) + diskDir_[1]*cos(yawRad),
                             diskDir_[2]);

        uniDiskDir_ /= mag(uniDiskDir_);
    }
}

void Foam::fv::wtActuatorSource::constructNodesPosList()
{
    // TODO: Improve this function to easy the cellID search when the actuator is moving
    // DynamicList<label> nodesCellId_old;
    // if (movingAD_){
    //     nodesCellId_old = nodesCellId_; //AO: store old values to seed findCell
    // }

    // label findSeed = -1;
    nodesPosList_.clear();
    nodesCellId_.clear();

    forAll(nodesList_, n)
    {
        point cartCoord = computeNodePosition(nodesList_[n][0], nodesList_[n][1]);

        nodesPosList_.append(cartCoord);

        // if (movingAD_) findSeed = nodesCellId_old[n];
        // cellId = mesh().findCell(cartCoord, findSeed); // Add seed to search when moving. The first time there is no seed!!
        label cellId = mesh().findCell(cartCoord); // TODO: cellID may be -1 if cell is in other proc
        nodesCellId_.append(cellId);
    }
}

point Foam::fv::wtActuatorSource::computeNodePosition(scalar rNode, scalar titaNode_rad)
{
    point Bi = vector(0, rNode * sin(titaNode_rad), rNode * cos(titaNode_rad));

    point BiDir = vector(-Bi[1] * uniDiskDir_[1], Bi[1] * uniDiskDir_[0], Bi[2]);

    return BiDir + diskPoint_;
}

void Foam::fv::wtActuatorSource::selectADcells()
{
    scalar sphereTol = selectRadialBoundFrac_ * maxR_;
    scalar planeTol = selectAxialBoundFrac_ * minSep_;
    diskCells_.clear();
    forAll(cells(), c)
    {
        vector arrow = mesh().cellCentres()[cells()[c]] - diskPoint_; // cell position relative to actuator center
        scalar dPlane = mag(uniDiskDir_ & arrow) / mag(uniDiskDir_); // cell distance from actuator plane
        scalar dSphere = mag(arrow); // cell distance from actuator center

        if ((dSphere <= sphereTol) and (dPlane <= planeTol))
        {
            diskCells_.append(cells()[c]);
        }
    }
    // Pout << "Actuator cells: " << diskCells_.size() << endl;
}

void Foam::fv::wtActuatorSource::computeUdWeights()
{
    scalar sigmaAxialWeight = sigmaAxialWeightFrac_ * minSep_;
    scalar maxRweight = maxRweightFrac_ * maxR_;
    scalar sigmaSphereWeight =  sigmaSphereWeightFrac_ * maxR_;

    scalar uWeight, vWeight, cWeight;
    if (unifWeight_) uWeight = 1; else
    {
        uWeight = 0;
        if (!sphereWeight_ and !axialWeight_) Info << "No weighting set!!" << endl;
    }

    weightADcells_.clear();
    weightedADvol_ = 0;
    forAll(diskCells_, c)
    {
        vector arrow = mesh().cellCentres()[diskCells_[c]] - diskPoint_; // cell position relative to actuator center
        scalar dPlane = mag(uniDiskDir_ & arrow) / mag(uniDiskDir_); // cell distance from actuator plane
        scalar dSphere = mag(arrow); // cell distance from actuator center

        if (dSphere <= maxRweight){
            if (!sphereWeight_ and !axialWeight_) vWeight = 0;
            else {
                vWeight = 1;
                if (sphereWeight_){
                    vWeight *= (1 / (sigmaSphereWeight * std::sqrt(M_PI))) * std::exp(-1 * std::pow((dSphere / sigmaSphereWeight), 2));
                }
                if (axialWeight_){
                    vWeight *= (1 / (sigmaAxialWeight * std::sqrt(M_PI))) * std::exp(-1 * std::pow((dPlane / sigmaAxialWeight), 2));
                }
            }
            cWeight = mesh().V()[diskCells_[c]] * (uWeight + vWeight);
        }
        else cWeight = 0;

        weightADcells_[diskCells_[c]] = cWeight;
        weightedADvol_ += cWeight;
    }
    if (parted_ || firstStep_ || mesh_.time().writeTime())
    {
        reduce(weightedADvol_, sumOp<scalar>());
    }
}

vector Foam::fv::wtActuatorSource::computeUd(
    const vectorField &U
)
{
    vector Ud = vector(0, 0, 0);
    forAll(diskCells_, c)
    {
        Ud += U[diskCells_[c]] * weightADcells_[diskCells_[c]] / weightedADvol_;
    }
    if (parted_ || firstStep_ || mesh_.time().writeTime())
    {
        reduce(Ud, sumOp<vector>());
    }

    return Ud;
}

void Foam::fv::wtActuatorSource::correctUd(
    scalar &UdDir,
    vector &Ud,
    scalar epsilon
)
{
    scalar axialInductionFactor = 0.5 - 0.5 * std::sqrt(1 - Ct_);
    scalar Ct_modified = Ct_ / std::pow(1 - axialInductionFactor, 2);
    scalar correctionScalar = std::pow(1 + (Ct_modified * epsilon) / (2 * std::sqrt(2 * M_PI) * maxR_), -1);
    UdDir *= correctionScalar;
    Ud *= correctionScalar;
}

static void interpolate(const List <scalar> &list, scalar value, label &i1, label &i2, scalar &ddx) {
    if(list.size() < 2 || list[0] >= value) {
        i1 = i2 = 0;
        ddx = 0;
        return;
    }

    for(i2 = 0; i2 < list.size(); i2++)
        if(list[i2] >= value) {
            i1 = i2 - 1;
            ddx = (value - list[i1]) / (list[i2] - list[i1]);
            return;
        }

    i1 = i2 = list.size() - 1;
    ddx = 0;
}

void Foam::fv::wtActuatorSource::interpolateWTCurves(scalar UdDir = VGREAT)
{
    label i1, i2;
    scalar ddx;

    if(UdDir == VGREAT) {
        // Interpolates Cp, Ct, omega and pitch in terms of Uref_
        interpolate(wtUrefList_, Uref_, i1, i2, ddx);
    }
    else {
        // Interpolates Uref_, Cp, Ct, omega and pitch in terms of Ud
        interpolate(wtUdList_, UdDir, i1, i2, ddx);
        // Interpolate Uref_
        Uref_ = ddx*(wtUrefList_[i2] - wtUrefList_[i1]) + wtUrefList_[i1];
    }

    // Interpolate Cp, Ct, omega and pitch
    Cp_ = ddx*(wtCpList_[i2] - wtCpList_[i1]) + wtCpList_[i1];
    Ct_ = ddx*(wtCtList_[i2] - wtCtList_[i1]) + wtCtList_[i1];
    omega_ = ddx*(wtOmegaList_[i2] - wtOmegaList_[i1]) + wtOmegaList_[i1];
    pitch_ = ddx*(wtPitchList_[i2] - wtPitchList_[i1]) + wtPitchList_[i1];
    lambda_ = omega_ * maxR_ / Uref_;
}

void Foam::fv::wtActuatorSource::computeWTmagnitudes()
{
    // TODO: Manage incompressible/compressible cases with density_
    Power_ =  diskArea_ * Cp_ * pow(Uref_, 3) / 2;
    Thrust_ = diskArea_ * Ct_ * pow(Uref_, 2) / 2;
    Torque_ = (omega_ < VSMALL) ? 0.0 : Power_ / omega_;
}

vector Foam::fv::wtActuatorSource::getNodeVelocity(
    label node,
    const volVectorField &U,
    const volTensorField *gradU
) const
{
    vector U_node = vector(1e9, 1e9, 1e9);
    label cellId = nodesCellId_[node];

    if (cellId != -1) // if the closer cell is in this procesor
    {
        U_node = U[cellId];

        if (gradU != NULL)
        {
            vector dx = nodesPosList_[node] - mesh().cellCentres()[cellId];
            vector dU = dx & (*gradU)[cellId];

            U_node += dU;
        }
    }

    if (parted_ || firstStep_ || mesh_.time().writeTime())
    {
        reduce(U_node, minOp<vector>());
    }

    if (mag(U_node) > 1e9) // Flag in case it does not find a cell
    {
        Info << U_node << endl;
        U_node = vector(10, 0, 0);
        Info << "OpenFOAM cell not found for actuator node: " << node << endl;
    }
    return U_node;
}

scalar Foam::fv::wtActuatorSource::induction1(scalar UdDir)
{
    // Computed induction according to Sørensen etal (2020) doi: 10.1016/j.renene.2019.09.134
    // Info << "Computes induction according to Sorensen formula, Uref_ in terms of Ud = " << UdDir << endl;
    return 2 * UdDir / (1 + sqrt(1 - Ct_));
}


scalar Foam::fv::wtActuatorSource::induction2(scalar UdDir)
{
    // Info << "Computes induction according to analytic formula, Uref_ in terms of Ud = " << UdDir << endl;
    return UdDir * Ct_ / Cp_;
}

void Foam::fv::wtActuatorSource::saveNodeForces(
    label node,
    vector Unode,
    scalar Faero_n,
    scalar Faero_t
) const
{
        // "Actuator name, time [s], node#, r [m], theta [rad], area, x, y, z,
        // Unode_x, Unode_y, Unode_z, Faero_n, Faero_t"
        (*outNodes) << name() << "," << mesh().time().value() << "," << node << ","
                    << nodesList()[node][0] << "," << nodesList()[node][1] << "," << nodesList()[node][2] << ","
                    << nodesPosList_[node][0] << "," << nodesPosList_[node][1] << "," << nodesPosList_[node][2] << ","
                    << Unode[0] << "," << Unode[1] << "," << Unode[2]
                    << "," << Faero_n << "," << Faero_t << std::endl;
}

void Foam::fv::wtActuatorSource::distributeActuatorForces(
    vectorField &Usource,
    vector Faero, // in (x, y, z) coordinates
    label node,
    tensor iTransform
) const
{
    scalar sphereTolDisk = sphereTolDiskFrac_ * maxR_;
    // Vector of Gaussian for force smearing (n, t, r)
    vector gaussianCoeff = distribGaussianCoeffs_ * minSep_;

    // If sphereTolDisk is 0 just put the force on the cell corresponding to node
    if (sphereTolDisk < VSMALL)
    {
        label cellId = nodesCellId_[node];
        if (cellId != -1)
        {
            Usource[cellId] += Faero;
        }
    }
    else // Distribute force with Gaussian kernel
    {
        std::map<label, float> weightCells;
        DynamicList<label> nodeCells;
        vector coordNode = nodesPosList_[node];
        scalar En = gaussianCoeff[0];
        scalar Et = gaussianCoeff[1];
        scalar Er = gaussianCoeff[2];

        scalar gaussScale = (1 / (En * Et * Er * pow(sqrt(M_PI), 3)));
        scalar nodeV = 0; // Cumulated volume of the node cells weighted for force distribution

        // loop over all the actuator cells to weight for force distribution
        forAll(diskCells_, c)
        {
            point cellCentre = mesh().cellCentres()[diskCells_[c]];
            vector arrow = cellCentre - coordNode;  // Segment from node to cell in Cartesian coordinates
            vector arrow_ntr = iTransform & arrow;  // Segment from node to cell in ntr coordinates

            // calculate the distances in blade coordinate system
            scalar dn = abs(arrow_ntr[0]);
            scalar dt = abs(arrow_ntr[1]);
            scalar dr = abs(arrow_ntr[2]);

            // Calculate weight if cell is inside the sphere (cut excess outside the disc)
            if (mag(cellCentre - diskPoint_) <= sphereTolDisk)
            {
                scalar gaussian = exp(-1 * (pow(dn / En, 2) + pow(dt / Et, 2) + pow(dr / Er, 2)));

                if (gaussian > 0.0001) // Disregard the tails of the Gaussian
                {
                    scalar weight = mesh().V()[diskCells_[c]] * gaussScale * gaussian;

                    weightCells[diskCells_[c]] = weight;
                    nodeV += weight;
                    nodeCells.append(diskCells_[c]);
                }
            }
        }
        if (parted_ || firstStep_ || mesh_.time().writeTime())
        {
            reduce(nodeV, sumOp<scalar>());
        }

        // Loop over cells subset to apply force to Usource
        forAll(nodeCells, c)
        {
            Usource[nodeCells[c]] += weightCells[nodeCells[c]] * Faero / nodeV;
        }

    }
}

template <class RhoFieldType>
void Foam::fv::wtActuatorSource::addActuatorForce(
    const volVectorField &U,
    const RhoFieldType &rho,
    vectorField &Usource
)
{
    // TODO: detect movement

    // Recompute if actuator moved
    if (movingAD_ or firstStep_){
        computeDiskOrientation(U);

        constructNodesPosList();

        // Create list of cells that belong to the disk
        selectADcells();

        computeUdWeights();
    }

    vector Ud = computeUd(U);
    scalar UdDir = mag(Ud[0] * uniDiskDir_[0] + Ud[1] * uniDiskDir_[1]);

    if (UdCorrection_) correctUd(UdDir, Ud, 3*minSep_);

    if (!calibration_)
    {
        switch (inductionTypeSelect)
        {
            case INDUCTION_TABLE:
            {
                interpolateWTCurves(UdDir);
                break;
            }
            case INDUCTION_SORENSEN:
            {
                Uref_ = induction1(UdDir);
                interpolateWTCurves();
                break;
            }
            case INDUCTION_ANALYTIC:
            {
                Uref_ = induction2(UdDir);
                interpolateWTCurves();
                break;
            }
            default:
            {
                interpolateWTCurves(UdDir);
            }
        }

        computeWTmagnitudes();
    }

    scalar Thrust_nodes = 0.0; // Thrust computed from nodes contribution
    scalar Torque_nodes = 0.0; // Torque computed from nodes contribution
    // TODO: accumulate thrust and torque from forces applied in cells
    bool flagWrite = mesh_.time().writeTime();

    actuatorModel_->applyForce(U, rho, Usource, Thrust_nodes, Torque_nodes, flagWrite);

    //----------------------------------------------------------------------------------------
    if ((Pstream::master()) and (flagWrite)) //In Master node
    {
        scalar t = mesh().time().value();
        if (saveLevel_)
        {
            // "Actuator_name , time [s], Uref [m/s], Ud [m/s], Cp, Ct, omega [rad/s], pitch [deg],"
            // "Power(Uref, Cp) [W], Thrust(Uref, Ct) [N], Torque [Nm]" << std::endl;
            (*outActuators) << name() << "," << t << "," << Uref_ << "," << mag(Ud) << "," << Cp_ << "," << Ct_ << ","
                            << omega_ << "," << pitch_ << "," << Power_ * density_ << "," << Thrust_ * density_
                            << "," << Torque_ * density_
                            << std::endl;

        }
        if (saveLevel_ > 1)
        {
            // "Actuator_name, time [s], Thrust_actuator [N], Torque_actuator [Nm], "
            // "Thrust_nodes [N], Torque_nodes [Nm]" << std::endl;
            (*outActuators2) << name() << "," << t << ","
                             << Thrust_ * density_ << "," << Torque_ * density_ << ","
                             << Thrust_nodes * density_ << "," << Torque_nodes * density_
                             << std::endl;
        }
    }

}
