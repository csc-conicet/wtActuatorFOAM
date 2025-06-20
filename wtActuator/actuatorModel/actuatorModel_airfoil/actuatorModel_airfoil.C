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

#include "actuatorModel_airfoil.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "fvc.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel_airfoil, 0);
    addToRunTimeSelectionTable(actuatorModel, actuatorModel_airfoil, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel_airfoil::actuatorModel_airfoil
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

void Foam::actuatorModel_airfoil::read(const dictionary& dict)
{
    actuatorModel::read(dict);
}

void Foam::actuatorModel_airfoil::init()
{
    Info << "    - initializing actuatorModel: airfoil" << endl;

    blade_ = new bladeModel(rotor_.getCoeffsDict().subDict("blade"));

    profiles_ = new profileModelList(rotor_.getCoeffsDict().subDict("profiles"));
    profiles_->connectBlades(blade_->profileName(), blade_->profileID());
}

void Foam::actuatorModel_airfoil::airfoilAnglesCoeffs(scalar &rNode, scalar &phi, scalar &chord, scalar &twist,
    scalar &alpha, scalar &beta, scalar &Cl, scalar &Cd
)
{
    label i1 = -1;
    label i2 = -1;
    scalar invDr = 0.0;
    blade_->interpolate(rNode, twist, chord, i1, i2, invDr);

    beta = rotor_.pitch() + twist; //AO: pitch debería ser global

    if (rotor_.omega() < 0)
    {
        beta = M_PI - beta;
    }

    alpha = phi - beta;

    // Determine profile data for this radius and angle of attack
    const label profile1 = blade_->profileID()[i1];
    const label profile2 = blade_->profileID()[i2];

    scalar Cd1 = 0.0;
    scalar Cl1 = 0.0;

    (*profiles_)[profile1].Cdl(alpha, Cd1, Cl1);

    scalar Cd2 = 0.0;
    scalar Cl2 = 0.0;
    (*profiles_)[profile2].Cdl(alpha, Cd2, Cl2);

    Cd = invDr * (Cd2 - Cd1) + Cd1;
    Cl = invDr * (Cl2 - Cl1) + Cl1;
}

// template <class RhoFieldType>
void Foam::actuatorModel_airfoil::applyForce(
    const volVectorField &U,
    // const RhoFieldType &rho, // TODO: use when templatized
    const geometricOneField &rho, // TODO: remove when templatized
    vectorField &Usource,
    scalar &Thrust_nodes, scalar &Torque_nodes, bool flagWrite
)
{
    volTensorField *gradU = NULL;
    if (gradInterpolation_) gradU = new volTensorField("gradU", fvc::grad(U));

    for(label node = 0; node < rotor_.nodesNumber(); node++)
    {
        scalar rNode = rotor_.nodesList()[node][0];
        if (abs(rNode) < VSMALL) continue; // No forces at central node

        scalar titaNode_rad = rotor_.nodesList()[node][1];
        scalar areaNode = rotor_.nodesList()[node][2];
        vector coordNode = rotor_.nodesPosList()[node];

        vector Unode = rotor_.getNodeVelocity(node, U, gradU);

        vector uniBladeDir, uniThetaDir;
        tensor iTransform;
        ntrVectors(coordNode, uniBladeDir, uniThetaDir, iTransform);
        vector Unode_ntr = - iTransform & Unode; // Change from global to ntr coordinates

        //----BEM CALCULATIONS-------------------
        scalar phi;
        vector Urel;
        nodeUrel(Unode_ntr, rNode, Urel, phi);

        scalar chord, twist, alpha, beta, Cl, Cd;
        airfoilAnglesCoeffs(rNode, phi, chord, twist, alpha, beta, Cl, Cd );

        // Reference aerodynamical force [N]
        scalar solidity;
        if (abs(rNode) < VSMALL) solidity = 1;
        else solidity = chord * rotor_.bladesNumber() / (2 * M_PI * rNode);

        scalar Faero = 0.5 * pow(mag(Urel), 2) * solidity * areaNode;

        scalar tiprootfactor = (this->*tipfactor_)(rNode, rotor_.lambda(), phi);
        tiprootfactor *= (this->*rootfactor_)(rNode, phi, rDist_);

        // Tangential and axial forces to be distributed
        scalar Faero_t = tiprootfactor * Faero * (Cl * sin(phi) - Cd * cos(phi));
        scalar Faero_n = tiprootfactor * Faero * (Cl * cos(phi) + Cd * sin(phi));

        rotor_.distributeActuatorForces(Usource, (Faero_t * uniThetaDir + Faero_n * rotor_.uniDiskDir()), node, iTransform);

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
