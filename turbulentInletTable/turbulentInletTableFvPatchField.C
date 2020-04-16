/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "turbulentInletTableFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::turbulentInletTableFvPatchField<Type>::turbulentInletTableFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    ranGen_(label(0)),
    fluctuationScale_(Zero),
    referenceField_(p.size()),
    alpha_(0.1),
    curTimeIndex_(-1),
    UmeanTable_()
{
}


template<class Type>
Foam::turbulentInletTableFvPatchField<Type>::turbulentInletTableFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    ranGen_(label(0)),
    fluctuationScale_(pTraits<Type>(dict.lookup("fluctuationScale"))),
    referenceField_("referenceField", dict, p.size()),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1),
    UmeanTable_(Function1<scalar>::New("UmeanTable", dict))
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(referenceField_);
    }
}


template<class Type>
Foam::turbulentInletTableFvPatchField<Type>::turbulentInletTableFvPatchField
(
    const turbulentInletTableFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(mapper(ptf.referenceField_)),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    UmeanTable_(ptf.UmeanTable_, false)
{
}


template<class Type>
Foam::turbulentInletTableFvPatchField<Type>::turbulentInletTableFvPatchField
(
    const turbulentInletTableFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    UmeanTable_(ptf.UmeanTable_, false)
{	Info << "333333333333333333333333333333333"<< endl;}


template<class Type>
Foam::turbulentInletTableFvPatchField<Type>::turbulentInletTableFvPatchField
(
    const turbulentInletTableFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    UmeanTable_(ptf.UmeanTable_, false)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::turbulentInletTableFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    m(referenceField_, referenceField_);
}


template<class Type>
void Foam::turbulentInletTableFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const turbulentInletTableFvPatchField<Type>& tiptf =
        refCast<const turbulentInletTableFvPatchField<Type>>(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void Foam::turbulentInletTableFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        Field<Type> randomField(this->size());

        forAll(patchField, facei)
        {
            randomField[facei] = ranGen_.sample01<Type>();
        }

        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
        scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

		const scalar t = this->db().time().timeOutputValue();
		scalar amplitude = UmeanTable_->value(t);///UmeanTable_.max();
		Info << "amplitude is "<< amplitude << endl;

        patchField = amplitude*
        (
            (1 - alpha_)*patchField
          + alpha_*
            (
                referenceField_
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<Type>::one,
                    fluctuationScale_
                )*mag(referenceField_)
            )
        );
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::turbulentInletTableFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "fluctuationScale", fluctuationScale_);
    writeEntry(os, "referenceField", referenceField_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "value", *this);
    writeEntry(os, UmeanTable_());
}


// ************************************************************************* //
