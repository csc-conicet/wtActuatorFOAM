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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline Foam::scalar Foam::fv::wtActuatorSource::maxR() const
{
    return maxR_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::diskArea() const
{
    return diskArea_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::Thrust() const
{
    return Thrust_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::Ct() const
{
    return Ct_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::Torque() const
{
    return Torque_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::Uref() const
{
    return Uref_;
}

inline Foam::dictionary Foam::fv::wtActuatorSource::getCoeffsDict() const
{
    return coeffs_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::nodesNumber() const
{
    return nodesNumber_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::bladesNumber() const
{
    return bladesNumber_;
}

inline Foam::DynamicList<Foam::vector> Foam::fv::wtActuatorSource::nodesList() const
{
    return nodesList_;
}

inline Foam::DynamicList<Foam::vector> Foam::fv::wtActuatorSource::nodesPosList() const
{
    return nodesPosList_;
}

inline Foam::DynamicList<Foam::scalar> Foam::fv::wtActuatorSource::rNodeList() const
{
    return rNodeList_;
}

inline Foam::vector Foam::fv::wtActuatorSource::uniDiskDir() const
{
    return uniDiskDir_;
}

inline Foam::point Foam::fv::wtActuatorSource::diskPoint() const
{
    return diskPoint_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::omega() const
{
    return omega_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::pitch() const
{
    return pitch_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::lambda() const
{
    return lambda_;
}

inline Foam::scalar Foam::fv::wtActuatorSource::density() const
{
    return density_;
}

// ************************************************************************* //
