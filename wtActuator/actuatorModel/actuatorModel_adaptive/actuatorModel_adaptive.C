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

#include "actuatorModel_adaptive.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_adaptive, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_adaptive, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_adaptive::actuatorModel_adaptive
(
    const fv::wtActuatorSource& rotor,
    const dictionary& dict
)
:
    actuatorModel(rotor, dict, typeName)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::actuatorModel_adaptive::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_adaptive::init()
{
    Info << "    - initializing actuatorModel: adaptive" << endl;
    Info << "       - defaulting to no tip factor" << endl;
    Info << "       - defaulting to no root factor" << endl;

    rNodeList_ = rotor_.rNodeList();

    //- forcesRadialDistribution
    List<List<scalar>> forcesRadialDistribution_;
    coeffs_.lookup("forcesRadialDistribution") >>  forcesRadialDistribution_;
    label table1size = forcesRadialDistribution_.size();

    //--- Define lists using the table forcesRadialDistribution
    //- Original lists
    List<scalar> Uref2List_orig(table1size);
    List<scalar> UinfList_orig(table1size);
    List<scalar> rList_orig(table1size);
    List<scalar> UdiList_orig(table1size);
    List<scalar> fnList_orig(table1size);
    List<scalar> ftList_orig(table1size);

    forAll(forcesRadialDistribution_, i)
    {
        Uref2List_orig[i] = forcesRadialDistribution_[i][0];
        UinfList_orig[i] = forcesRadialDistribution_[i][1];
        rList_orig[i] = forcesRadialDistribution_[i][2];
        UdiList_orig[i] = forcesRadialDistribution_[i][3];
        fnList_orig[i] = forcesRadialDistribution_[i][4];
        ftList_orig[i] = forcesRadialDistribution_[i][5];
    }

    // Count how many radius sections are in the original table
    label nR_orig = 1;
    DynamicList<scalar> rListOnly_orig;
    rListOnly_orig.append(rList_orig[0]);
    while (rList_orig[nR_orig] > rList_orig[nR_orig - 1])
    {
        rListOnly_orig.append(rList_orig[nR_orig]);
        nR_orig++;
    }
    nR_ = rNodeList_.size(); // Radius sections in the new table

    // Count how many unique Uinf values are in the original table
    label i = 0;
    UinfOnlyList_.append(UinfList_orig[i]);
    while (UinfList_orig[i] < UinfList_orig[i +  nR_orig])
    {
        i += nR_orig;
        UinfOnlyList_.append(UinfList_orig[i]);
    }
    nUinf_orig_ = UinfOnlyList_.size();

    // Count how many unique Uref values are in the original table
    UrefOnlyList_.append(Uref2List_orig[i]);
    for (label i = 1; i < table1size; i++)
    {
        if (Uref2List_orig[i] > Uref2List_orig[i - 1])
        {
            UrefOnlyList_.append(Uref2List_orig[i]);
        }
    }
    nUref_orig_ = UrefOnlyList_.size();

    // Set up new table lists with the corresponding size
    label table2size = nR_ * nUinf_orig_ * nUref_orig_;
    UdiList_.setSize(table2size);
    fnList_.setSize(table2size);
    ftList_.setSize(table2size);

    // Set up new list with radial positions of nodes from the old table that are before each actuator node.
    posrList_.setSize(nR_);
    for (label i = 0; i < nR_; i++)
    {
        label i1, i2;
        posInList(rListOnly_orig, rNodeList_[i], i1, i2);

        posrList_[i] = i1;
    }

    // Fill up the new table
    for (label i = 0; i < nUref_orig_; i++)
    {
        for (label j = 0; j < nUinf_orig_; j++)
        {
            for (label k = 0; k < nR_; k++)
            {
                label pos_new = i * (nUinf_orig_ * nR_) + j * nR_ + k;
                label pos_orig = i * (nUinf_orig_ * nR_orig) + j * nR_orig + posrList_[k];

                scalar ddr = 0;
                if (posrList_[k] + 1 < nR_orig)
                {
                    ddr = (rNodeList_[k] - rList_orig[posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);
                }

                UdiList_[pos_new] = UdiList_orig[pos_orig] + ddr * (UdiList_orig[pos_orig + 1] - UdiList_orig[pos_orig]) ;
                fnList_[pos_new] = fnList_orig[pos_orig] + ddr * (fnList_orig[pos_orig + 1] - fnList_orig[pos_orig]);
                if (rNodeList_[k] < VSMALL)
                {
                    ftList_[pos_new] = 0;
                }
                else
                {
                    ftList_[pos_new] = ftList_orig[pos_orig] + ddr * (ftList_orig[pos_orig + 1] - ftList_orig[pos_orig]);
                }
            }
        }
    }

    rNodePosList_.setSize(rotor_.nodesNumber());
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar dist = VGREAT;
        for(label posrNode = 0; posrNode < rNodeList_.size(); posrNode++)
        {
            scalar dist_new = fabs(rNode - rNodeList_[posrNode]);
            if(dist_new < dist)
            {
                dist = dist_new;
                rNodePosList_[node] = posrNode;
            }
        }
    }
}

void Foam::actuatorModel_adaptive::posInTableUdi
(
    const label &posUref,
    const scalar &posrNode,
    const scalar &magUnode,
    label &pos1,
    label &pos2
)
{
    pos1 = 0;
    pos2 = 0;
    if (nUinf_orig_ != 1)
    {
        while ((pos2 < nUinf_orig_) and
               (UdiList_[posUref * nUinf_orig_ * nR_ + pos2 * nR_ + posrNode] < magUnode))
        {
            pos2++;
        }

        if (pos2 == 0)
        {
            pos1 = pos2;
            pos2++;
        }
        else if (pos2 == nUinf_orig_)
        {
            pos2 = nUinf_orig_ - 1;
            pos1 = pos2 - 1;
        }
        else
        {
            pos1 = pos2 - 1;
        }
    }
}

// template <class RhoFieldType>
void Foam::actuatorModel_adaptive::applyForce(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{
    volTensorField *gradU = NULL;
    if (gradInterpolation_) gradU = new volTensorField("gradU", fvc::grad(U));

    // Search position below and above in Uref list
    label pos1, pos2;
    posInList(UrefOnlyList_, rotor_.Uref(), pos1, pos2);

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];
        vector coordNode = rotor_.nodesPosList()[node];

        vector Unode = rotor_.getNodeVelocity(node, U, gradU);
        scalar magUnode = mag(Unode);

        vector uniBladeDir, uniThetaDir;
        tensor iTransform;
        if (fabs(rNode) < VSMALL)
        {
            uniThetaDir = {0, 0, 1};
            iTransform = tensor::one;
        }
        else
        {
            ntrVectors(coordNode, uniBladeDir, uniThetaDir, iTransform);
        }

        //---- Calculate Uinf, fn and ft for Uref position 1 ----------------------------------
        label pos1Udi1, pos1Udi2;
        posInTableUdi(pos1, rNodePosList_[node], magUnode, pos1Udi1, pos1Udi2);

        label pos1UdiList1 = pos1 * nUinf_orig_ * nR_ + pos1Udi1 * nR_ + rNodePosList_[node];
        label pos1UdiList2 = pos1 * nUinf_orig_ * nR_ + pos1Udi2 * nR_ + rNodePosList_[node];

        scalar ddUdi1 = 0;
        if (nUinf_orig_ != 1)
        {
            ddUdi1 = (magUnode - UdiList_[pos1UdiList1])
                     / (UdiList_[pos1UdiList2] - UdiList_[pos1UdiList1]);
        }

        scalar fn1 = ddUdi1 * ((fnList_[pos1UdiList2]) - (fnList_[pos1UdiList1]))
                     + (fnList_[pos1UdiList1]);

        scalar ft1 = ddUdi1 * ((ftList_[pos1UdiList2]) - (ftList_[pos1UdiList1]))
                     + (ftList_[pos1UdiList1]);

        //---- Calculate Uinf, fn and ft for Uref position 2 ----------------------------------
        label pos2Udi1, pos2Udi2;
        posInTableUdi(pos2, rNodePosList_[node], magUnode, pos2Udi1, pos2Udi2);

        label pos2UdiList1 = pos2 * nUinf_orig_ * nR_ + pos2Udi1 * nR_ + rNodePosList_[node];
        label pos2UdiList2 = pos2 * nUinf_orig_ * nR_ + pos2Udi2 * nR_ + rNodePosList_[node];

        scalar ddUdi2 = 0;
        if (nUinf_orig_ != 1)
        {
            ddUdi2 = (magUnode - UdiList_[pos2UdiList1])
                     / (UdiList_[pos2UdiList2] - UdiList_[pos2UdiList1]);
        }

        scalar fn2 = ddUdi2 * ((fnList_[pos2UdiList2]) - (fnList_[pos2UdiList1]))
                     + (fnList_[pos2UdiList1]);

        scalar ft2 = ddUdi2 * ((ftList_[pos2UdiList2]) - (ftList_[pos2UdiList1]))
                     + (ftList_[pos2UdiList1]);

        //---- Interpolate Uinf, fn and ft from values of position 1 and 2 ----------------------------------
        scalar ddUref = 0;
        if (pos1 != pos2)
        {
            ddUref = ((rotor_.Uref() - UrefOnlyList_[pos1])
                     / (UrefOnlyList_[pos2] - UrefOnlyList_[pos1]));
        }

        scalar Faero_t = (ddUref * (ft2 - ft1) + ft1) * areaNode;
        scalar Faero_n = (ddUref * (fn2 - fn1) + fn1) * areaNode;

        rotor_.distributeActuatorForces(Usource,
                                        (Faero_t * uniThetaDir + Faero_n * rotor_.uniDiskDir()),
                                        node,
                                        iTransform);

        if ((Pstream::master()) and (flagWrite)) //In Master node
        {
            if (saveLevel_ > 1)
            {
                Thrust_nodes += Faero_n;
                Torque_nodes += (-((rotor_.diskPoint() - coordNode) ^ uniThetaDir) & rotor_.uniDiskDir()) * Faero_t;
            }
            if (saveNodeForces_)
            {
                rotor_.saveNodeForces(node, Unode, Faero_n / areaNode, Faero_t / areaNode);
            }
        }
    } // close node loop
    if(gradU != NULL) delete gradU;
}

// ************************************************************************* //
