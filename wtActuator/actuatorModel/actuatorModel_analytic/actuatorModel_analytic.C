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

#include "actuatorModel_analytic.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"
#include "fvCFD.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_analytic, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_analytic, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_analytic::actuatorModel_analytic
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

void Foam::actuatorModel_analytic::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_analytic::init()
{
    Info << "    - initializing actuatorModel: analytic" << endl;
}


scalar Foam::actuatorModel_analytic::induction1(scalar &UdDir)
{
    return 2 * UdDir / (1 + sqrt(1 - rotor_.Ct()));
}


// template <class RhoFieldType>
void Foam::actuatorModel_analytic::applyForce(
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
    scalar sum_a2factor = 0;
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];

        if (abs(rNode) < VSMALL)
        {
            vector Unode = rotor_.getNodeVelocity(node, U, gradU);
            scalar tiprootfactor = (this->*rootfactor_)(rNode, M_PI, rDist_);

            scalar a2factor = tiprootfactor * areaNode;
            sum_a2factor += a2factor;

            uniThetaDir_list.append(vector(0, 0, 0));
            iTransform_list.append(tensor::one);
            tiprootfactor_list.append(tiprootfactor);
            scalar Udi = mag(Unode);
            Udi_list.append(Udi);
            Uinfi_list.append(induction1(Udi));
            if (saveNodeForces_)
            {
                Unode_list.append(Unode);
            }
        }
        else
        {
            scalar titaNode_rad = rotor_.nodesList()[node][1];
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
            scalar a2factor = tiprootfactor * areaNode;
            sum_a1factor += a1factor;
            sum_a2factor += a2factor;

            uniThetaDir_list.append(uniThetaDir);
            iTransform_list.append(iTransform);
            tiprootfactor_list.append(tiprootfactor);
            scalar Udi = mag(Unode);
            Udi_list.append(Udi);
            Uinfi_list.append(induction1(Udi));
            if (saveNodeForces_)
            {
                Unode_list.append(Unode);
            }
        }
    }
    if(gradU != NULL) delete gradU;

    scalar q0 = ((sqrt(pow(rotor_.lambda() * sum_a2factor, 2) +
                       pow(rotor_.maxR(), 2) * rotor_.diskArea() * rotor_.Ct() * sum_a1factor)) -
                 rotor_.lambda() * sum_a2factor) /
                (pow(rotor_.maxR(), 2) * sum_a1factor);

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];
        vector coordNode = rotor_.nodesPosList()[node];

        scalar Faero_n = 0;
        scalar Faero_t = 0;
        if (abs(rNode) < VSMALL)
        {
            Faero_n = areaNode * pow(Uinfi_list[node], 2) * (q0 * rotor_.lambda() * tiprootfactor_list[node]);
        }
        else
        {
            Faero_n = areaNode * pow(Uinfi_list[node], 2)
                            * (q0 * rotor_.lambda() * tiprootfactor_list[node]
                               + pow(q0 * tiprootfactor_list[node] * rotor_.maxR() / rNode, 2) / 2);

            Faero_t = areaNode * Uinfi_list[node] * Udi_list[node] * q0 * tiprootfactor_list[node] * rotor_.maxR() / rNode;
        }

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
