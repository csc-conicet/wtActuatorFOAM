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

#include "actuatorModel_numeric.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_numeric, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_numeric, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_numeric::actuatorModel_numeric
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

void Foam::actuatorModel_numeric::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_numeric::init()
{
    Info << "    - initializing actuatorModel: numeric" << endl;
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
    List<scalar> fnList_orig(table1size);
    List<scalar> ftList_orig(table1size);

    forAll(forcesRadialDistribution_, i)
    {
        Uref2List_orig[i] = forcesRadialDistribution_[i][0];
        UinfList_orig[i] = forcesRadialDistribution_[i][1];
        rList_orig[i] = forcesRadialDistribution_[i][2];
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

    scalar Urefaux = Uref2List_orig[0]; 
    // scalar Uinfaux = UinfList_orig[0];
    UrefOnlyList_.append(Uref2List_orig[0]);
    DynamicList<label> posInUref;
    posInUref.append(0);

    for (label j = 0; j < Uref2List_orig.size(); j = j + nR_orig)
    {
        if (Uref2List_orig[j] > Urefaux)
        {
            Urefaux = Uref2List_orig[j];
            UrefOnlyList_.append(Urefaux);
            posInUref.append(j);
        }
    }
    posInUref.append(Uref2List_orig.size());
    nUref_orig_ = UrefOnlyList_.size();

    DynamicList<label> posUrefUinf;
    posUrefUinf.setSize(posInUref.size() - 1);
    for (label i = 0; i < posInUref.size() - 1; i++)
    {
        if (posInUref[i + 1] - posInUref[i] == nR_orig)
        {
            posUrefUinf[i] = posInUref[i];
        }
        else
        {
            scalar diffU = VGREAT;
            for (label j = posInUref[i]; j < posInUref[i + 1]; j = j + nR_orig)
            {
                scalar diffU_new = fabs(Uref2List_orig[j] - UinfList_orig[j]);
                if(diffU_new < diffU)
                {
                    diffU = diffU_new;
                    posUrefUinf[i] = j;
                }
            }
        }
    }

    // Set up new table lists with the corresponding size
    label table2size = nR_ * posUrefUinf.size();
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
    for (label i = 0; i < posUrefUinf.size(); i++)
    {
        for (label k = 0; k < nR_; k++)
        {
            label pos_new = i * nR_ + k;
            // label pos_orig = i * (nUinf_orig_ * nR_orig) + posUrefInUinf_[i] * nR_orig + posrList_[k];
            label pos_orig = posUrefUinf[i] + posrList_[k];

            scalar ddr = 0;
            if (posrList_[k] + 1 < nR_orig)
            {
                ddr = (rNodeList_[k] - rList_orig[posrList_[k]]) / (rList_orig[posrList_[k] + 1] - rList_orig[posrList_[k]]);
            }

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

    sumFaero_t_.setSize(nUref_orig_, 0.0);
    sumFaero_n_.setSize(nUref_orig_, 0.0);
    rNodePosList_.setSize(rotor_.nodesNumber());
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];
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

        for (label i = 0; i < nUref_orig_; i++)
        {
            label pos_new = i * nR_ + rNodePosList_[node];
            sumFaero_n_[i] += fnList_[pos_new] * areaNode;
            sumFaero_t_[i] += ftList_[pos_new] * rNode * areaNode;
        }
    }
}

// template <class RhoFieldType>
void Foam::actuatorModel_numeric::applyForce(
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

        //---- Interpolate Uinf, fn and ft from values of position 1 and 2 ----------------------------------
        scalar ddUref = 0;
        if (pos1 != pos2)
        {
            ddUref = ((rotor_.Uref() - UrefOnlyList_[pos1])
                     / (UrefOnlyList_[pos2] - UrefOnlyList_[pos1]));
        }

        scalar fn1 = fnList_[pos1 * nR_ + rNodePosList_[node]] / sumFaero_n_[pos1];
        scalar fn2 = fnList_[pos2 * nR_ + rNodePosList_[node]] / sumFaero_n_[pos2];

        scalar ft1 = ftList_[pos1 * nR_ + rNodePosList_[node]] / sumFaero_t_[pos1];
        scalar ft2 = ftList_[pos2 * nR_ + rNodePosList_[node]] / sumFaero_t_[pos2];

        scalar Faero_n = (ddUref * (fn2 - fn1) + fn1) * areaNode * rotor_.Thrust();
        scalar Faero_t = (ddUref * (ft2 - ft1) + ft1) * areaNode * rotor_.Torque();

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
