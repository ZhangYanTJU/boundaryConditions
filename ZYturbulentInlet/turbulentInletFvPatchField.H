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

Class
    Foam::ZYturbulentInletFvPatchField

Description
    This boundary condition generates a fluctuating inlet condition by adding
    a random component to a reference (mean) field.

    \f[
        x_p = x_{ref} + random*x_{ref}*s
    \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | reference patch values
        s       | fluctuation scale
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        fluctuationScale | RMS fluctuation scale (fraction of mean) | yes |
        referenceField | reference (mean) field | yes        |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            ZYturbulentInlet;
        fluctuationScale 0.1;
        referenceField  uniform (10 0 0);
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    ZYturbulentInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ZYturbulentInletFvPatchField_H
#define ZYturbulentInletFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ZYturbulentInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/


class ZYturbulentInletFvPatchField
:
    public fixedValueFvPatchField<vector>
{
    // Private data

        //- Random number generator
        Random ranGen_;

        //- Fluctuation scale
        vector fluctuationScale_;

        //- Reference field
        Field<vector> referenceField_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        vector Umean_;
        scalar width_;
        scalar midRadius_;
        point center_;

public:

    //- Runtime type information
    TypeName("ZYturbulentInlet");


    // Constructors

        //- Construct from patch and internal field
        ZYturbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ZYturbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given ZYturbulentInletFvPatchField
        //  onto a new patch
        ZYturbulentInletFvPatchField
        (
            const ZYturbulentInletFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ZYturbulentInletFvPatchField
        (
            const ZYturbulentInletFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new ZYturbulentInletFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ZYturbulentInletFvPatchField
        (
            const ZYturbulentInletFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new ZYturbulentInletFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const vector& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            vector& fluctuationScale()
            {
                return fluctuationScale_;
            }

            //- Return the reference field
            const Field<vector>& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<vector>& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*#ifdef NoRepository
    #include "turbulentInletFvPatchField.C"
#endif
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
