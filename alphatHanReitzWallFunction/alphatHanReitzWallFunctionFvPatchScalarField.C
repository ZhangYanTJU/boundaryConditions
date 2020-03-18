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

#include "alphatWallFunctionFvPatchScalarField.H"
#include "compressibleTurbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatWallFunctionFvPatchScalarField::alphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85)
{}


alphatWallFunctionFvPatchScalarField::alphatWallFunctionFvPatchScalarField
(
    const alphatWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_)
{}


alphatWallFunctionFvPatchScalarField::alphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85))
{}


alphatWallFunctionFvPatchScalarField::alphatWallFunctionFvPatchScalarField
(
    const alphatWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_)
{}


alphatWallFunctionFvPatchScalarField::alphatWallFunctionFvPatchScalarField
(
    const alphatWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const compressibleTurbulenceModel& turbModel =
        db().lookupObject<compressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                compressibleTurbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);
        
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<scalarField> tnutw = turbModel.nut(patchi);

    // original wall function
    //operator==(rhow*tnutw/Prt_);
    
    // Han Reitz wall function
    scalarField& alphatw = *this;
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
        const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const volScalarField& alpha = db().lookupObject<volScalarField> ("thermo:alpha");
    alpha.writeMinMax(Info);
    const scalarField& alphaw = alpha.boundaryField()[patchi];
    //const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    //const scalarField& alphaw = talphaw();    
    const scalar Cmu25 = pow025(nutw.Cmu());
    const volScalarField& T = db().lookupObject<volScalarField> ("T");
    const scalarField& Tw = T.boundaryField()[patchi];
    scalarField qwByCp(this->size(),pTraits<scalar>::zero); 
    
    forAll(alphatw, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        scalar uTau = Cmu25*sqrt(k[cellI]);
        scalar yPlus = uTau*y[faceI]/(nuw[faceI]);
            
        if (yPlus > nutw.yPlusLam())
        {
            qwByCp[faceI] = rhow[faceI]*uTau*T[cellI]*log(T[cellI]/Tw[faceI])/(2.1*log(yPlus)+2.513);
            alphatw[faceI] = qwByCp[faceI]*y[faceI]/(T[cellI]-Tw[faceI]+SMALL)-alphaw[faceI];
        }
    }    

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "Prt", Prt_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
