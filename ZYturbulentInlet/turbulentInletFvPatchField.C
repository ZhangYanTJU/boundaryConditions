/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "turbulentInletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::ZYturbulentInletFvPatchField::ZYturbulentInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    ranGen_(label(0)),
    fluctuationScale_(Zero),
    referenceField_(p.size()),
    Umean_   (0,0,0),
    width_   (0),
    midRadius_   (0),
    center_   (0,0,1)      
{}



Foam::ZYturbulentInletFvPatchField::ZYturbulentInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    ranGen_(label(0)),
    fluctuationScale_(pTraits<vector>(dict.lookup("fluctuationScale"))),
    referenceField_("referenceField", dict, p.size()),
    Umean_   (0,0,0),
    width_   (0),
    midRadius_   (0),
    center_   (0,0,1)    
{
    Umean_ = dict.lookup("Umean");
    width_ = readScalar(dict.lookup("width"));
    midRadius_ = readScalar(dict.lookup("midRadius"));
    center_ = dict.lookup("center");
    
    scalar delta_ = width_/2;

    const vectorField& Cf = this->patch().Cf();
    Info<<"width_======"<<width_<<endl;
    Info<<"delta_======"<<delta_<<endl;
    
    referenceField_ = 1.218*Umean_*pow( 1. - mag(mag(Cf - center_) - midRadius_)/1.01/delta_, 1./7.);
    if (dict.found("value"))
    {
        fixedValueFvPatchField<vector>::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<vector>::operator==(referenceField_);
    }
}



Foam::ZYturbulentInletFvPatchField::ZYturbulentInletFvPatchField
(
    const ZYturbulentInletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(mapper(ptf.referenceField_)),
    Umean_   (ptf.Umean_),
    width_   (ptf.width_),
    midRadius_   (ptf.midRadius_),
    center_   (ptf.center_)    
{}



Foam::ZYturbulentInletFvPatchField::ZYturbulentInletFvPatchField
(
    const ZYturbulentInletFvPatchField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    Umean_   (ptf.Umean_),
    width_   (ptf.width_),
    midRadius_   (ptf.midRadius_),
    center_   (ptf.center_)     
{}



Foam::ZYturbulentInletFvPatchField::ZYturbulentInletFvPatchField
(
    const ZYturbulentInletFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    Umean_   (ptf.Umean_),
    width_   (ptf.width_),
    midRadius_   (ptf.midRadius_),
    center_   (ptf.center_)     
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::ZYturbulentInletFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
    //referenceField_.autoMap(m);
}



void Foam::ZYturbulentInletFvPatchField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const ZYturbulentInletFvPatchField& tiptf =
        refCast<const ZYturbulentInletFvPatchField>(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}



void Foam::ZYturbulentInletFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<vector>& patchField = *this;

        Field<vector> randomField(this->size());

        forAll(patchField, facei)
        {
            randomField[facei] = ranGen_.sample01<vector>();
        }

        //patchField = referenceField_+randomField*referenceField_*fluctuationScale_;
        patchField = referenceField_+cmptMultiply(randomField,fluctuationScale_)*mag(referenceField_);
        
        //OF original
        /*
        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
        scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

        patchField =
            (1 - alpha_)*patchField
          + alpha_*
            (
                referenceField_
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<Type>::one,
                    fluctuationScale_
                )*mag(referenceField_)
            );  
        */            
        
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}



void Foam::ZYturbulentInletFvPatchField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("fluctuationScale")
        << fluctuationScale_ << token::END_STATEMENT << nl;
    writeEntry(os, "referenceField", referenceField_);
    writeEntry(os, "value", *this);
    os.writeKeyword("Umean")<<Umean_<<token::END_STATEMENT<<nl;
    os.writeKeyword("width")<<width_<<token::END_STATEMENT<<nl;
    os.writeKeyword("midRadius")<<midRadius_<<token::END_STATEMENT<<nl;
    os.writeKeyword("center")<<center_<<token::END_STATEMENT<<nl;    
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        ZYturbulentInletFvPatchField
    );
}

// ************************************************************************* //
