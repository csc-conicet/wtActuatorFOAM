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

#include "actuatorModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuatorModel, 0);
    defineRunTimeSelectionTable(actuatorModel, dictionary);
}

const Foam::Enum<Foam::actuatorModel::tipFactorName>
Foam::actuatorModel::tipFactorName_
({
    { tipFactorName::TIP_NONE, "none" },
    { tipFactorName::TIP_SHEN, "shen" },
    { tipFactorName::TIP_GLAUERT, "glauert" },
    { tipFactorName::TIP_PRANDTL, "prandtl" },
});

const Foam::Enum<Foam::actuatorModel::rootFactorName>
Foam::actuatorModel::rootFactorName_
({
    { rootFactorName::ROOT_NONE, "none" },
    { rootFactorName::ROOT_GLAUERT, "glauert" },
    { rootFactorName::ROOT_SORENSEN, "sorensen" },
});

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuatorModel::actuatorModel
(
    const fv::wtActuatorSource& rotor,
    const dictionary& dict,
    const word& name
)
:
    rotor_(rotor),
    name_(name),
    coeffs_()
{
    read(dict);

    coeffs_.readIfPresent("rootfactordistance", rDist_);
    coeffs_.readIfPresent("gradInterpolation", gradInterpolation_);
    saveLevel_ = coeffs_.lookupOrDefault("saveLevel", 1);
    saveNodeForces_ = coeffs_.lookupOrDefault("saveNodeForces", false);

    if (tipFactorName_.readIfPresent("tipfactor", coeffs_, tipFactorSelect))
    {
        switch (tipFactorSelect)
        {
            case TIP_SHEN:
            {
                Info << "    - setting Shen's tip factor" << endl;
                tipfactor_ = &Foam::actuatorModel::tipfactor_Shen;
                break;
            }
            case TIP_GLAUERT:
            {
                Info << "    - setting Glauert's tip factor" << endl;
                tipfactor_ = &Foam::actuatorModel::tipfactor_Glauert;
                break;
            }
            case TIP_PRANDTL:
            {
                Info << "    - setting Prandtl's tip factor" << endl;
                tipfactor_ = &Foam::actuatorModel::tipfactor_Prandtl;
                break;
            }
            case TIP_NONE:
            {
                Info << "    - setting no tip factor" << endl;
                tipfactor_ = &Foam::actuatorModel::tipfactor_none;
            }
        }
    }
    else
    {
        Info << "    - defaulting to no tip factor" << endl;
        tipFactorSelect = tipFactorName::TIP_NONE;
        tipfactor_ = &Foam::actuatorModel::tipfactor_none;
    }

    if (rootFactorName_.readIfPresent("rootfactor", coeffs_, rootFactorSelect))
    {
        switch (rootFactorSelect)
        {
            case ROOT_GLAUERT:
            {
                Info << "    - setting Glauert's root factor" << endl;
                rootfactor_ = &Foam::actuatorModel::rootfactor_Glauert;
                break;
            }
            case ROOT_SORENSEN:
            {
                Info << "    - setting Sorensen's root factor" << endl;
                rootfactor_ = &Foam::actuatorModel::rootfactor_Sorensen;
                break;
            }
            case ROOT_NONE:
            {
                Info << "    - setting no root factor" << endl;
                rootfactor_ = &Foam::actuatorModel::rootfactor_none;
            }
        }
    }
    else
    {
        Info << "    - defaulting to no root factor" << endl;
        rootFactorSelect = rootFactorName::ROOT_NONE;
        rootfactor_ = &Foam::actuatorModel::rootfactor_none;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::actuatorModel::read(const dictionary& dict)
{
    coeffs_ = dict.optionalSubDict(name_ + "Coeffs");
}

scalar Foam::actuatorModel::tipfactor_Shen(
    const scalar rNode, const scalar tsr, const scalar phi
)
{
    // Shen tip correction factor:
    scalar c1 = 0.125;
    scalar c2 = 27;
    scalar c3 = 0.1;

    scalar tipfactor = 1;
    // f = 1/0 -> infty, then exp(-infty)=0, then acos(0) = pi/2, then F=1
    if (phi > VSMALL and rNode > VSMALL) {
        scalar tipfactor_f = (rotor_.bladesNumber() / 2) * (rotor_.maxR() - rNode) / (rNode * sin(phi));
        scalar g = exp(-c1 * (rotor_.bladesNumber() * tsr - c2)) + c3;
        scalar g2 = exp(-g * tipfactor_f);
        if ((g2 > -1) and (g2 < 1))
        {
            tipfactor = 2 / M_PI * acos(g2);
        }
    }
    return tipfactor;
}

scalar Foam::actuatorModel::tipfactor_Glauert(
    const scalar rNode, const scalar tsr, const scalar phi
)
{
    // Glauert tip correction factor:

    scalar tipfactor = 1;
    // f = 1/0 -> infty, then exp(-infty)=0, then acos(0) = pi/2, then F=1
    if (phi > VSMALL and rNode > VSMALL) {
        scalar tipfactor_f = (rotor_.bladesNumber() / 2) * (rotor_.maxR() - rNode) / (rNode * sin(phi));
        tipfactor = 2 / M_PI * acos(exp(-tipfactor_f));
    }
    return tipfactor;
}

scalar Foam::actuatorModel::tipfactor_Prandtl(
    // Prandtl tip correction factor:
    const scalar rNode, const scalar tsr, const scalar phi
)
{
    scalar tipfactor_f = (rotor_.bladesNumber() / 2) * (1 - rNode / rotor_.maxR()) * sqrt(1 + pow(tsr, 2));
    scalar tipfactor = 2 / M_PI * acos(exp(-tipfactor_f));
    return tipfactor;
}

scalar Foam::actuatorModel::tipfactor_none(
    const scalar rNode, const scalar tsr, const scalar phi
)
{
    return 1.0;
}


scalar Foam::actuatorModel::rootfactor_Glauert(
    const scalar rNode, const scalar phi, const scalar rDist_
)
{
    // Glauert root correction factor:

    scalar rootFactor = 0;
    if ((rNode / rotor_.maxR()) > rDist_)
    {
        if (phi < VSMALL) rootFactor = 1;
        else
        {
            scalar rootfactor_f = (rotor_.bladesNumber() / 2) * (1 - rDist_ * (rotor_.maxR() / rNode)) / sin(phi);
            if (rootfactor_f > VSMALL)
            {
                rootFactor = 2 / M_PI * acos(exp(-rootfactor_f));
            }
        }
    }
    return rootFactor;
}

scalar Foam::actuatorModel::rootfactor_Sorensen(
    const scalar rNode, const scalar phi, const scalar rDist_
)
{
    // Sorensen root correction factor:

    scalar a = 2.335;
    scalar b = 4;

    return (1 - exp(-a * pow((rNode / (rotor_.maxR() * rDist_)), b)));
}

scalar Foam::actuatorModel::rootfactor_none(
    const scalar rNode, const scalar tsr, const scalar phi
)
{
    return 1.0;
}

void Foam::actuatorModel::ntrVectors(vector &coordNode, vector &uniBladeDir, vector &uniThetaDir, tensor &iTransform)
{
    uniBladeDir = coordNode - rotor_.diskPoint();
    uniBladeDir /= mag(uniBladeDir);

    uniThetaDir =  rotor_.uniDiskDir() ^ uniBladeDir;
    uniThetaDir /=  mag(uniThetaDir);

    // change of coordinate system
    tensor Transform(-rotor_.uniDiskDir()[0], uniThetaDir[0], uniBladeDir[0],
                     -rotor_.uniDiskDir()[1], uniThetaDir[1], uniBladeDir[1],
                     -rotor_.uniDiskDir()[2], uniThetaDir[2], uniBladeDir[2]);

    iTransform = inv(Transform);
}

void Foam::actuatorModel::nodeUrel(vector &Unode_ntr, scalar &rNode, vector &Urel, scalar &phi)
{
    // velocities in the profile coordinates
    scalar U_n = Unode_ntr[0];
    scalar U_t = Unode_ntr[1];
    scalar U_r = Unode_ntr[2];

    // scalar phi;
    if (abs(rNode * rotor_.omega() - U_t) < VSMALL) phi = M_PI;
    else (phi = Foam::atan(U_n / (rNode * rotor_.omega() - U_t)));

    Urel = vector(U_n, U_t - rNode * rotor_.omega(), 0.0);

    if (((U_t - rNode * rotor_.omega()) <= 0) and (U_n >= 0))
    {
    }

    if (((U_t - rNode * rotor_.omega()) <= 0) and (U_n < 0))
    {
    }

    if (((U_t - rNode * rotor_.omega()) > 0) and (U_n >= 0))
    {
        phi = phi + M_PI;
    }

    if (((U_t - rNode * rotor_.omega()) > 0) and (U_n < 0))
    {
        phi = phi - M_PI;
    }
}

void Foam::actuatorModel::posInList
(
    const List<scalar> &inputList,
    const scalar &searchValue,
    label &pos1,
    label &pos2
)
{
    pos2 = 0;
    label nElem = inputList.size();

    if (nElem == 1)
    {
        pos1 = pos2;
    }
    else
    {
        while ((pos2 < nElem) and (inputList[pos2] < searchValue))
        {
            pos2++;
        }

        if (pos2 == 0)
        {
            pos1 = 0;
            pos2 = 1;
        }
        else if (pos2 == nElem)
        {
            pos2 = nElem - 1;
            pos1 = pos2 - 1;
        }
        else
        {
            pos1 = pos2 - 1;
        }
    }
}


// ************************************************************************* //
