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
    Foam::incompressible::RASModels::tpphiLowReRoughWallTPFvPatchScalarField

Description
    Boundary condition for turbulent (kinematic) viscosity when using wall
    functions, based on turbulence kinetic energy.
    - replicates OpenFOAM v1.5 (and earlier) behaviour

SourceFiles
    tpphiLowReRoughWallTPFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef tpphiLowReRoughWallTPFvPatchScalarField_H
#define tpphiLowReRoughWallTPFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
             Class tpphiLowReRoughWallTPFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class tpphiLowReRoughWallTPFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;
		
		//- Ks roughness sandgrain height
		scalar ks_;

		//- Phi/K ratio for rough surfaces
		scalar pkC_;

        label bz_;
		
        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;
		
		//- Choose tpphi wall condition
		word rType_;


    // Protected member functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the Y+ at the edge of the laminar sublayer
        virtual scalar calcYPlusLam(const scalar kappa, const scalar E) const;


        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("tpphiLowReRoughWallTP");


    // Constructors

        //- Construct from patch and internal field
        tpphiLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tpphiLowReRoughWallTPFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tpphiLowReRoughWallTPFvPatchScalarField
        //  onto a new patch
        tpphiLowReRoughWallTPFvPatchScalarField
        (
            const tpphiLowReRoughWallTPFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tpphiLowReRoughWallTPFvPatchScalarField
        (
            const tpphiLowReRoughWallTPFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new tpphiLowReRoughWallTPFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tpphiLowReRoughWallTPFvPatchScalarField
        (
            const tpphiLowReRoughWallTPFvPatchScalarField&,
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
                new tpphiLowReRoughWallTPFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
			
			//- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );
			
			//- Return reference to tpphi rough type calc
            word& rType()
            {
                return rType_;
            }

			//- Return reference to rough phi/k ratio
            scalar& pkC()
            {
                return pkC_;
            }


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
