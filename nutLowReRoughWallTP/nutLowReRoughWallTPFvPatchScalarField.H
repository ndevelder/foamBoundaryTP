/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::
    nutLowReRoughWallTPFvPatchScalarField

Description
    Wall function boundary condition for rough walls

SourceFiles
    nutSpalartAllmarasStandardWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutLowReRoughWallTPFvPatchScalarField_H
#define nutLowReRoughWallTPFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
Class nutLowReRoughWallTPFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutLowReRoughWallTPFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private data

        // Roughness model parameters

            //- Height
            scalar ks_;

            label bz_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFudgeFactor_;
			
			//- Choose nut wall condition
			word nutExp_;


    // Protected member functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutLowReRoughWallTP");


    // Constructors

        //- Construct from patch and internal field
        nutLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutLowReRoughWallTPFvPatchScalarField
        //  onto a new patch
        nutLowReRoughWallTPFvPatchScalarField
        (
            const nutLowReRoughWallTPFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutLowReRoughWallTPFvPatchScalarField
        (
            const nutLowReRoughWallTPFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutLowReRoughWallTPFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nutLowReRoughWallTPFvPatchScalarField
        (
            const nutLowReRoughWallTPFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutLowReRoughWallTPFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the roughness height
            const scalar& roughnessHeight() const
            {
                return ks_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight()
            {
                return ks_;
            }


            //- Return the roughness constant scale
            const scalar& roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            const scalar& roughnessFudgeFactor() const
            {
                return roughnessFudgeFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //  adjustment
            scalar& roughnessFudgeFactor()
            {
                return roughnessFudgeFactor_;
            }

			//- Return reference to the roughness scale factor to allow
            word& nutExp()
            {
                return nutExp_;
            }

    // I-O

        // Evaluation functions

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
