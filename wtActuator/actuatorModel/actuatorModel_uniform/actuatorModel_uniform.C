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

#include "actuatorModel_uniform.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_uniform, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_uniform, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_uniform::actuatorModel_uniform
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

void Foam::actuatorModel_uniform::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_uniform::init()
{
    Info << "    - initializing actuatorModel: uniform" << endl;
}

// template <class RhoFieldType>
void Foam::actuatorModel_uniform::applyForce(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{

    if ((rootFactorSelect != rootFactorName::ROOT_NONE) or (tipFactorSelect != tipFactorName::TIP_NONE))
    {
        applyForce_tiproot(U, rho, Usource, Thrust_nodes, Torque_nodes, flagWrite);
    }
    else
    {
        applyForce_no_tiproot(U, rho, Usource, Thrust_nodes, Torque_nodes, flagWrite);
    }

}

void Foam::actuatorModel_uniform::applyForce_tiproot(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{
    scalar center_area = 0.0;
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        if (abs(rNode) < VSMALL)
        {
            center_area += rotor_.nodesList()[node][2];
        }
    }

    scalar Thrust_area = rotor_.Thrust() / rotor_.diskArea();
    scalar Torque_area = rotor_.Torque() / (rotor_.diskArea() - center_area);

    volTensorField *gradU = NULL;
    if (gradInterpolation_) gradU = new volTensorField("gradU", fvc::grad(U));

    DynamicList<vector> Unode_list;
    DynamicList<vector> uniThetaDir_list;
    DynamicList<scalar> F_ni_list;
    DynamicList<scalar> F_ti_list;
    DynamicList<tensor> iTransform_list;

    scalar sumF_nXfactor = 0;
    scalar sumF_t = 0;
    scalar sumF_tXfactor = 0;

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];

        if (abs(rNode) < VSMALL)
        {
            scalar tiprootfactor = (this->*rootfactor_)(rNode, M_PI, rDist_);

            scalar F_niXfactor = tiprootfactor * areaNode;
            sumF_nXfactor += F_niXfactor;

            uniThetaDir_list.append(vector(0, 0, 0));
            F_ni_list.append(F_niXfactor);
            F_ti_list.append(0);
            iTransform_list.append(tensor::one);
            if (saveNodeForces_)
            {
                vector Unode = rotor_.getNodeVelocity(node, U, gradU);
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

            scalar F_niXfactor = tiprootfactor * areaNode;
            scalar F_tiXfactor = tiprootfactor * areaNode / rNode;
            sumF_nXfactor += F_niXfactor;
            sumF_tXfactor += F_tiXfactor;
            sumF_t += areaNode / rNode;

            uniThetaDir_list.append(uniThetaDir);
            F_ni_list.append(F_niXfactor);
            F_ti_list.append(F_tiXfactor);
            iTransform_list.append(iTransform);
            if (saveNodeForces_)
            {
                Unode_list.append(Unode);
            }

        }
    }
    if(gradU != NULL) delete gradU;

    scalar scaleFactor_n = rotor_.diskArea() / sumF_nXfactor;
    scalar scaleFactor_t = sumF_t / sumF_tXfactor;

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        vector coordNode = rotor_.nodesPosList()[node];
        scalar areaNode = rotor_.nodesList()[node][2];

        scalar Faero_n = scaleFactor_n * F_ni_list[node] * Thrust_area;
        scalar Faero_t = scaleFactor_t * F_ti_list[node] * Torque_area;

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

void Foam::actuatorModel_uniform::applyForce_no_tiproot(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{
    scalar center_area = 0.0;
    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        if (abs(rNode) < VSMALL)
        {
            center_area += rotor_.nodesList()[node][2];
        }
    }

    scalar Thrust_area = rotor_.Thrust() / rotor_.diskArea();
    scalar Torque_area = rotor_.Torque() / (rotor_.diskArea() - center_area);

    volTensorField *gradU = NULL;
    if (saveNodeForces_ and gradInterpolation_) gradU = new volTensorField("gradU", fvc::grad(U));

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        scalar areaNode = rotor_.nodesList()[node][2];

        scalar Faero_n = areaNode * Thrust_area;

        if ((abs(rNode) < VSMALL) or (rotor_.omega() < VSMALL))
        {
            rotor_.distributeActuatorForces(Usource, (Faero_n * rotor_.uniDiskDir()), node, tensor::one);

            if (flagWrite)
            {
                if (saveNodeForces_)
                {
                    vector Unode;
                    Unode = rotor_.getNodeVelocity(node, U, gradU);
                    if (Pstream::master()) //In Master proc
                    {
                        rotor_.saveNodeForces(node, Unode, Faero_n / areaNode, 0.0);
                    }
                }

                if (Pstream::master() and (saveLevel_ > 1)) //In Master proc
                    {
                        Thrust_nodes += Faero_n;
                    }
            }
        }
        else
        {
            scalar Faero_t = areaNode * Torque_area / rNode;

            vector coordNode = rotor_.nodesPosList()[node];

            vector uniBladeDir, uniThetaDir;
            tensor iTransform;
            ntrVectors(coordNode, uniBladeDir, uniThetaDir, iTransform);

            rotor_.distributeActuatorForces(Usource, (Faero_t * uniThetaDir + Faero_n * rotor_.uniDiskDir()), node, iTransform);

            if (flagWrite)
            {
                if (saveNodeForces_)
                {
                    vector Unode;
                    Unode = rotor_.getNodeVelocity(node, U, gradU);
                    if (Pstream::master()) //In Master proc
                    {
                        rotor_.saveNodeForces(node, Unode, Faero_n / areaNode, Faero_t / areaNode);
                    }
                }

                if (Pstream::master() and (saveLevel_ > 1)) //In Master proc
                {
                    Thrust_nodes += Faero_n;
                    Torque_nodes += mag(rotor_.diskPoint() - coordNode) * Faero_t;
                }
            }
        }
    }
    if(gradU != NULL) delete gradU;
}

// ************************************************************************* //
