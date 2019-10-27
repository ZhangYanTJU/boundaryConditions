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

Class
    Foam::decayingTurbulenceFvPatchVectorField_ZY

Description
    Generate a valid (in term of Reynolds stress and lenght scales) turbulent
    inflow for LES, DES or DNS. This method use the Lund transformation
    (Lund, Wu and Squires, 1998) to construct a fluctuating velocity field.
    Check if your R field is valid against this transformation! If not, use the
    robust Lund transformation. It sets Rij=0 (with i not equal to j).

    Parameters:
    * direction: It is +1 or −1 depending on wheather the flow goes in positive or negative x-direction,respectively.
    * writeSourceFields: () write the interpolated ref, R, L and Lund fields if true. Debug purpose.
    * robustLundTransformation: (default=false) Use the robust Lund transformation algorithm.
    * LField: The integral length field.
    * refField: The external velocity field.
    * RField: The Reynolds stress tensor field.
    * value: The value of the field; at the initial time it should be equal to the refField field.

    inlet
    {
        type decayingTurbulenceInflowGeneratorMod;
        direction 1;
        LField nonuniform List<scalar> ...;
        refField nonuniform List<vector> ...;
        RField nonuniform List<symmTensor> ...;
        value nonuniform List<vector> ...;
    }

SourceFiles
    decayingTurbulenceFvPatchVectorField_ZY.C

\*---------------------------------------------------------------------------*/

#ifndef decayingTurbulenceFvPatchVectorField_H
#define decayingTurbulenceFvPatchVectorField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "SLList.H"
#include "pointToPointPlanarInterpolation.H"

#include "decayingVorton.H"

namespace Foam
{


class decayingTurbulenceFvPatchVectorField_ZY : public fixedValueFvPatchField<vector>
{
private:

    Field<scalar> LField_;
    Field<vector> refField_;
    Field<symmTensor> RField_;
    label curTimeIndex_;
    SLList<decayingVorton> vortons_;
    Field<tensor> Lund_;
    Field<symmTensor> R_;
    label ind_;
    label direction_;
    Field<scalar> nVortField_; //debug field, counts number of vortons
    label inletShape_;

    vector Umean_;
    scalar width_;
    scalar midRadius_;
    point center_;

    scalar Radius_;

    //- 2D interpolation
    //autoPtr<pointToPointPlanarInterpolation> mapperPtr_;
    //- Fraction of perturbation (fraction of bounding box) to add
    scalar perturb_;

public:

    TypeName("decayingTurbulenceInflowGeneratorMod_ZY");

    decayingTurbulenceFvPatchVectorField_ZY
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    decayingTurbulenceFvPatchVectorField_ZY
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    decayingTurbulenceFvPatchVectorField_ZY
    (
        const decayingTurbulenceFvPatchVectorField_ZY&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    decayingTurbulenceFvPatchVectorField_ZY
    (
        const decayingTurbulenceFvPatchVectorField_ZY&
    );

    decayingTurbulenceFvPatchVectorField_ZY
    (
        const decayingTurbulenceFvPatchVectorField_ZY&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchField<vector> > clone() const
    {
        return tmp<fvPatchField<vector> >(new decayingTurbulenceFvPatchVectorField_ZY(*this));
    }

    virtual tmp<fvPatchField<vector> > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchField<vector> >(new decayingTurbulenceFvPatchVectorField_ZY(*this, iF));
    }

    virtual void autoMap(const fvPatchFieldMapper&);

    virtual void rmap(const fvPatchField<vector>&, const labelList&);

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;

private:

    void doUpdate();
    void readInletData();
    void writeSourceFields();

};
}

#endif
