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
    Foam::incompressible::RASModels::epsilonLowReRoughWallTPFvPatchScalarField

Description
    Boundary condition for epsilon when using wall functions
    - calculates epsilon and G
    - epsilon values added directly into the matrix to act as a constraint

SourceFiles
    epsilonLowReRoughWallTPFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonLowReRoughWallTPFvPatchScalarField_H
#define epsilonLowReRoughWallTPFvPatchScalarField_H

#include "fixedInternalValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class epsilonLowReRoughWallTPFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonLowReRoughWallTPFvPatchScalarField
:
    public fixedInternalValueFvPatchScalarField
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulence kinetic energy field
        word kName_;

        //- Name of turbulence generation field
        word GName_;

        //- Name of laminar viscosity field
        word nuName_;

        //- Name of turbulent viscosity field
        word nutName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;
		
		scalar ks_;
		
		scalar sr_;
		
		word epsType_;


    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("epsilonLowReRoughWallTP");


    // Constructors

        //- Construct from patch and internal field
        epsilonLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonLowReRoughWallTPFvPatchScalarField
        //  onto a new patch
        epsilonLowReRoughWallTPFvPatchScalarField
        (
            const epsilonLowReRoughWallTPFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonLowReRoughWallTPFvPatchScalarField
        (
            const epsilonLowReRoughWallTPFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonLowReRoughWallTPFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonLowReRoughWallTPFvPatchScalarField
        (
            const epsilonLowReRoughWallTPFvPatchScalarField&,
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
                new epsilonLowReRoughWallTPFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes  = Pstream::blocking
            );


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
