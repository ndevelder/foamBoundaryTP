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

\*---------------------------------------------------------------------------*/

#include "nutLowReRoughWallTPFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
nutLowReRoughWallTPFvPatchScalarField::calcNut() const
{
	// Get patch/mesh indices
    const label patchI = patch().index();
	const fvMesh& mesh = patch().boundaryMesh().mesh();
    
	// Load Turbulence Model object
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");	
    const scalarField& y = rasModel.y()[patchI];
		
	const dictionary& rasDictionary = db().lookupObject<IOdictionary>("RASProperties");
	dictionary tpCoeffDict(rasDictionary.subDict("turbulentPotentialCoeffs"));
	word tslimiter(tpCoeffDict.lookup("tslimiter"));
	scalar nRMax = readScalar(tpCoeffDict.lookup("nutRatMax"));
	
	const volScalarField& kr = mesh.lookupObject<volScalarField>("k");
	const volScalarField& epsr = mesh.lookupObject<volScalarField>("epsilon");
	const volScalarField& tpr = mesh.lookupObject<volScalarField>("tpphi");

	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");	
	
	const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
	const scalarField magGradUw = mag(Uw.snGrad());
    
	tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
	scalarField& nutw = tnutw();
	
	scalar T = 0.0;

    forAll(kr.boundaryField()[patchI], faceI)
    {
		label faceCellI = patch().faceCells()[faceI];
		if(nutExp_ == "default"){
			
		   if(tslimiter == "true"){
			T = max(kr[faceCellI]/(epsr[faceCellI]+SMALL), 6.0*sqrt(nuw[faceI]/(epsr[faceCellI]+SMALL)));
		   }
		   
		   if(tslimiter == "false"){
			T = kr[faceCellI]/epsr[faceCellI];
		   }           
		   
		   nutw[faceI] = 0.21*kr.boundaryField()[patchI][faceI]*tpr.boundaryField()[patchI][faceI]*T;
		   nutw[faceI] = min(nutw[faceI],nRMax*nuw[faceI]);
		}
		
		if(nutExp_ == "ksquared"){
           nutw[faceI] = 0.09*kr.boundaryField()[patchI][faceI]*kr[faceCellI]/(epsr[faceCellI]+SMALL);
		}		
    }

	return tnutw;
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    roughnessHeight_(pTraits<scalar>::zero),
    roughnessConstant_(pTraits<scalar>::zero),
    roughnessFudgeFactor_(pTraits<scalar>::zero),
	nutExp_("default")
{}


nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const nutLowReRoughWallTPFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    roughnessHeight_(ptf.roughnessHeight_),
    roughnessConstant_(ptf.roughnessConstant_),
    roughnessFudgeFactor_(ptf.roughnessFudgeFactor_),
	nutExp_(ptf.nutExp_)
{}


nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    roughnessHeight_(readScalar(dict.lookup("roughnessHeight"))),
    roughnessConstant_(readScalar(dict.lookup("roughnessConstant"))),
    roughnessFudgeFactor_(readScalar(dict.lookup("roughnessFudgeFactor"))),
	nutExp_(dict.lookup("nutExp"))
{}


nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const nutLowReRoughWallTPFvPatchScalarField& rwfpsf
)
:
    nutWallFunctionFvPatchScalarField(rwfpsf),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFudgeFactor_(rwfpsf.roughnessFudgeFactor_),
	nutExp_(rwfpsf.nutExp_)
{}


nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const nutLowReRoughWallTPFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(rwfpsf, iF),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFudgeFactor_(rwfpsf.roughnessFudgeFactor_),
	nutExp_(rwfpsf.nutExp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void nutLowReRoughWallTPFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeKeyword("roughnessHeight")
        << roughnessHeight_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessConstant")
        << roughnessConstant_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessFudgeFactor")
        << roughnessFudgeFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("nutExp")
        << nutExp_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutLowReRoughWallTPFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
