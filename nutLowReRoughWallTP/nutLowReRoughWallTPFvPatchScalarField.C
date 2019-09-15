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
	dictionary bCoeffDict(rasDictionary.subDict("boundaryCoeffs"));
	word tslimiter(bCoeffDict.lookup("tslimiter"));
	scalar nRMax = readScalar(bCoeffDict.lookup("nutRatMax"));
	scalar cMuBc = readScalar(bCoeffDict.lookup("cMu"));
	
	const volScalarField& kr = db().lookupObject<volScalarField>("k");
    const volScalarField& ksr = db().lookupObject<volScalarField>("kSqrt");
	const volScalarField& epsr = db().lookupObject<volScalarField>("epsilon");
	const volScalarField& epsh = db().lookupObject<volScalarField>("epsHat");
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
    const volVectorField& tpsr = db().lookupObject<volVectorField>("tppsi");
    const volScalarField& lam = db().lookupObject<volScalarField>("lambda");
	
	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");	
	
	const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
	const scalarField magGradUw = mag(Uw.snGrad());

    const volVectorField& tppsiw = db().lookupObject<volVectorField>("tppsi");
    const scalarField magPsiSqrw = (tppsiw.boundaryField()[patchI] & tppsiw.boundaryField()[patchI])*sqr(kr.boundaryField()[patchI]);
    const scalarField magPsiw = sqrt(magPsiSqrw + SMALL);
    const scalarField psiDV = magPsiw/(magGradUw + SMALL);
    
	tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
	scalarField& nutw = tnutw();
	
	scalar T = 0.0;
    scalar minT = 0.0;
    scalar pSqr = 0.0;
	scalar nW = 0.0;
	scalar phW = 0.0;

    forAll(kr.boundaryField()[patchI], faceI)
    {
		label faceCellI = patch().faceCells()[faceI];
        scalar nuEffw = nuw[faceI]+nutw[faceI];
        scalar nuPsiw = nuw[faceI]+psiDV[faceI];

        scalar utauw = sqrt(nuPsiw*magGradUw[faceI] + SMALL);
        scalar kPlus = ks_*utauw/nuw[faceI];
		
		if(nutExp_ == "default"){			
           minT = 6.0*sqrt(nuw[faceI]/epsr[faceCellI]);
		   T = max(kr[faceCellI]/(epsr[faceCellI]+SMALL),minT);
		   nutw[faceI] = cMuBc*kr.boundaryField()[patchI][faceI]*tpr.boundaryField()[patchI][faceI]*T;
		   nutw[faceI] = min(nutw[faceI],nRMax*nuw[faceI]);
		}

        if(nutExp_ == "lam"){
           pSqr = tpsr.boundaryField()[patchI][faceI] & tpsr.boundaryField()[patchI][faceI];
           nutw[faceI] = (0.6*(0.12 + 0.37*lam.boundaryField()[patchI][faceI]) + 0.4*pSqr)*kr.boundaryField()[patchI][faceI]*kr[faceCellI]/(epsh[faceCellI]+SMALL);
        }
		
		if(nutExp_ == "sdiv"){
		   nutw[faceI] = epsr.boundaryField()[patchI][faceI]/sqr(magGradUw[faceI]);
		}
		
		if(nutExp_ == "psi"){
		   nutw[faceI] = mag(tpsr.boundaryField()[patchI][faceI])*kr.boundaryField()[patchI][faceI]/magGradUw[faceI];
		}

        if(nutExp_ == "ks"){
           nutw[faceI] = pow(0.09,0.25)*ksr.boundaryField()[patchI][faceI]*0.41*0.03*ks_*min(1.0, kPlus/90.0);
        }

        if(nutExp_ == "kss"){
           nutw[faceI] = pow(0.09,0.25)*ksr.boundaryField()[patchI][faceI]*0.41*0.03*ks_*pow(min(1.0, kPlus/90.0),2.0);
        }
    
	
		
        //nutw[faceI] = 1e-10;
        nW = nutw[faceI];
        phW = tpr.boundaryField()[patchI][faceI];		
		
    }
	
	//Info << "nutw: " << nW << " | nutT: " << T << " | phW: " << phW << endl;

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
    ks_(pTraits<scalar>::zero),
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
    ks_(ptf.ks_),
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
    ks_(dict.lookupOrDefault<scalar>("ks",0.0)),
    roughnessConstant_(dict.lookupOrDefault<scalar>("roughnessConstant",0.0)),
    roughnessFudgeFactor_(dict.lookupOrDefault<scalar>("roughnessFudgeFactor",0.0)),
	nutExp_(dict.lookupOrDefault<word>("nutExp","default"))
{}


nutLowReRoughWallTPFvPatchScalarField::
nutLowReRoughWallTPFvPatchScalarField
(
    const nutLowReRoughWallTPFvPatchScalarField& rwfpsf
)
:
    nutWallFunctionFvPatchScalarField(rwfpsf),
    ks_(rwfpsf.ks_),
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
    ks_(rwfpsf.ks_),
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
    os.writeKeyword("ks")
        << ks_ << token::END_STATEMENT << nl;
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
