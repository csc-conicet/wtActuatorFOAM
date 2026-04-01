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

#include "actuatorModel_genanalytic.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"
#include "fvCFD.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_genanalytic, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_genanalytic, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_genanalytic::actuatorModel_genanalytic
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

void Foam::actuatorModel_genanalytic::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_genanalytic::init()
{
    Info << "    - initializing actuatorModel: genanalytic" << endl;
}


scalar Foam::actuatorModel_genanalytic::induction1(scalar &UdDir)
{
    return 2 * UdDir / (1 + sqrt(1 - rotor_.Ct()));
}

scalar Foam::actuatorModel_genanalytic::induction_1(scalar Uref)
{
    return Uref * (1 + sqrt(1 - rotor_.Ct())) / 2;
}

// template <class RhoFieldType>
void Foam::actuatorModel_genanalytic::applyForce(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{
    volTensorField *gradU = NULL;
    if (gradInterpolation_) gradU = new volTensorField("gradU", fvc::grad(U));

    DynamicList<vector> Unode_list;
    DynamicList<vector> uniThetaDir_list;
    DynamicList<tensor> iTransform_list;
    DynamicList<scalar> tiprootfactor_list;
    DynamicList<scalar> Udi_list;
    DynamicList<scalar> Uinfi_list;

    scalar sum_a1factor = 0;
    scalar sum_a21factor = 0;
    scalar sum_a22factor = 0;
    scalar sum_a01factor = 0;
    scalar sum_a02factor = 0;
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];

        if (fabs(rNode) < VSMALL * rotor_.maxR())
        {
            vector Unode = rotor_.getNodeVelocity(node, U, gradU);
            scalar tiprootfactor = (this->*rootfactor_)(rNode, M_PI, rDist_);

            scalar a21factor = tiprootfactor * areaNode;
            scalar a22factor = tiprootfactor * tiprootfactor * areaNode;
            sum_a21factor += a21factor;
            sum_a22factor += a22factor;

            uniThetaDir_list.append(vector(0, 0, 0));
            iTransform_list.append(tensor::one);
            tiprootfactor_list.append(tiprootfactor);
            scalar Udi = Unode & rotor_.uniDiskDir();
            Udi_list.append(Udi);
            Uinfi_list.append(induction1(Udi));
            if (saveNodeForces_)
            {
                Unode_list.append(Unode);
            }
        }
        else
        {
            vector coordNode = rotor_.nodesPosList()[node];

            vector Unode = rotor_.getNodeVelocity(node, U, gradU);

            vector uniBladeDir, uniThetaDir;
            tensor iTransform;
            ntrVectors(coordNode, uniBladeDir, uniThetaDir, iTransform);
            vector Unode_ntr = - iTransform & Unode; // Change from global to ntr coordinates

            scalar phi;
            vector Urel;
            nodeUrel(Unode_ntr, rNode, Urel, phi);

            scalar tiprootfactor = (this->*tipfactor_)(rNode, rotor_.lambda(), phi);
            tiprootfactor *= (this->*rootfactor_)(rNode, phi, rDist_);

            scalar a1factor = tiprootfactor * tiprootfactor * areaNode / (rNode * rNode);
            scalar a21factor = tiprootfactor * areaNode;
            scalar a22factor = tiprootfactor * tiprootfactor * areaNode;
            scalar a01factor = tiprootfactor * rNode * rNode * areaNode;
            scalar a02factor = tiprootfactor * tiprootfactor * rNode * rNode * areaNode;
            sum_a1factor += a1factor;
            sum_a21factor += a21factor;
            sum_a22factor += a22factor;
            sum_a01factor += a01factor;
            sum_a02factor += a02factor;

            uniThetaDir_list.append(uniThetaDir);
            iTransform_list.append(iTransform);
            tiprootfactor_list.append(tiprootfactor);
            scalar Udi = Unode & rotor_.uniDiskDir();
            Udi_list.append(Udi);
            Uinfi_list.append(induction1(Udi));
            if (saveNodeForces_)
            {
                Unode_list.append(Unode);
            }
        }
    }
    if(gradU != NULL) delete gradU;

    scalar Ctdiff = (rotor_.Ct_rated() - rotor_.Ct()) / rotor_.Ct_rated();

    scalar S0;
    if(Ctdiff > 0)
    {
        S0 = 0.08 * pow(Ctdiff, 3);
    }
    else
    {
        S0 = 0.05 * Ctdiff;
    }

    scalar q0 = ((sqrt(pow(S0 * sum_a22factor - rotor_.lambda() * sum_a21factor, 2) +
                       sum_a1factor * (2 * rotor_.lambda() * S0 * sum_a01factor -
                                       pow(S0, 2) * sum_a02factor +
                                       pow(rotor_.maxR(), 2) * rotor_.diskArea() * rotor_.Ct()))) +
                 (S0 * sum_a22factor - rotor_.lambda() * sum_a21factor)) /
                (pow(rotor_.maxR(), 2) * sum_a1factor);

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];
        vector coordNode = rotor_.nodesPosList()[node];

        scalar Faero_n = 0;
        scalar Faero_t = 0;
        scalar u_theta = 0;
        if (fabs(rNode) > VSMALL * rotor_.maxR())
        {
            u_theta = Uinfi_list[node] * tiprootfactor_list[node] * (q0 * rotor_.maxR() / rNode - S0 * rNode / rotor_.maxR());
        }

        Faero_n = areaNode * u_theta * (rotor_.lambda() * Uinfi_list[node] * rNode / rotor_.maxR() + u_theta / 2);

        Faero_t = areaNode * Udi_list[node] * u_theta;
    
        rotor_.distributeActuatorForces(Usource,
                                        (Faero_t * uniThetaDir_list[node] + Faero_n * rotor_.uniDiskDir()),
                                        node,
                                        iTransform_list[node]);

        if ((Pstream::master()) and (flagWrite)) //In Master node
        {
            if (saveLevel_ > 1)
            {
                Thrust_nodes += Faero_n;
                Torque_nodes += (-((rotor_.diskPoint() - coordNode) ^ uniThetaDir_list[node]) & rotor_.uniDiskDir()) * Faero_t;
            }
            if (saveNodeForces_)
            {
                rotor_.saveNodeForces(node, Unode_list[node], Faero_n / areaNode, Faero_t / areaNode);
            }
        }
    } // close node loop
}

// ************************************************************************* //
