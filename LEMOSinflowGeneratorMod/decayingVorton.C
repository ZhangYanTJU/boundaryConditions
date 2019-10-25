/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#include "decayingVorton.H"
#include "Random.H"

extern Foam::Random ranGen;

namespace Foam
{

decayingVorton::decayingVorton()
:
    length_  (0.0),
    location_(pTraits<vector>::zero),
    omega_   (pTraits<vector>::zero),
    velocity_(pTraits<vector>::zero),
    xmax_    (0.0)
{
}

decayingVorton::decayingVorton(Istream& s)
:
    length_  (readScalar(s)),
    location_(s),
    omega_   (s),
    velocity_(s),
    xmax_    (readScalar(s))
{
}

decayingVorton::decayingVorton(scalar length, const vector& location, const vector& velocity, scalar xmax)
:
    length_  (length),
    location_(location),
    omega_   (pTraits<vector>::zero),
    velocity_(velocity),
    xmax_    (xmax)
{
    //vector omg = 2*ranGen.vector01() - pTraits<vector>::one;
    //vector omg = 2*vector(ranGen.scalar01(),ranGen.scalar01(),ranGen.scalar01()) - pTraits<vector>::one;
    /*vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = ranGen.scalar01();
    }
    vector omg = 2*rndVec - pTraits<vector>::one;*/
    vector omg = 2*ranGen.sample01<vector>() - pTraits<vector>::one;
    omega_     = omg/mag(omg);
}

autoPtr<decayingVorton> decayingVorton::New(Istream& s)
{
    Info<<"reading"<<endl;
    return autoPtr<decayingVorton>(new decayingVorton(s));
}

void decayingVorton::operator=(const decayingVorton& vt)
{
    if (this == &vt)
        FatalErrorIn("decaying::operator=(const decaying&)")<< "Attempted assignment to self"<<abort(FatalError);

    length_   = vt.length_;
    location_ = vt.location_;
    omega_    = vt.omega_;
    velocity_ = vt.velocity_;
    xmax_     = vt.xmax_;
}

bool decayingVorton::operator!=(const decayingVorton& vt) const
{
    return (length_   != vt.length_)   ||
           (location_ != vt.location_) ||
           (omega_    != vt.omega_)    ||
           (velocity_ != vt.velocity_) ||
           (xmax_     != vt.xmax_);
}

Ostream& operator<<(Ostream& s, const decayingVorton& vt)
{
    s<<vt.length_<<endl;
    s<<vt.location_<<endl;
    s<<vt.omega_<<endl;
    s<<vt.velocity_<<endl;
    s<<vt.xmax_<<endl;

    return s;
}

Istream& operator>>(Istream& s, decayingVorton& vt)
{
    scalar len(readScalar(s));
    vector loc(s);
    vector omg(s);
    vector vel(s);
    scalar xmax(readScalar(s));

    vt.length_   = len;
    vt.location_ = loc;
    vt.omega_    = omg;
    vt.velocity_ = vel;
    vt.xmax_     = xmax;

    return s;
}

}
