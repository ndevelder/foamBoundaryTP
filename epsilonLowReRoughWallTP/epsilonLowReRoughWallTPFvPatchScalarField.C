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

#include "epsilonLowReRoughWallTPFvPatchScalarField.H"
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonLowReRoughWallTPFvPatchScalarField::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorIn("epsilonLowReRoughWallTPFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonLowReRoughWallTPFvPatchScalarField::epsilonLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
	ks_(0.0),
	sr_(0.235),
	epsType_("rough")
{
    checkType();
}


epsilonLowReRoughWallTPFvPatchScalarField::epsilonLowReRoughWallTPFvPatchScalarField
(
    const epsilonLowReRoughWallTPFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
	ks_(ptf.ks_),
	sr_(ptf.sr_),
	epsType_(ptf.epsType_)
{
    checkType();
}


epsilonLowReRoughWallTPFvPatchScalarField::epsilonLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
	ks_(dict.lookupOrDefault<scalar>("ks", 0.0)),
	sr_(dict.lookupOrDefault<scalar>("sr", 0.235)),
	epsType_(dict.lookupOrDefault<word>("epsType", "rough"))
{
    checkType();
}


epsilonLowReRoughWallTPFvPatchScalarField::epsilonLowReRoughWallTPFvPatchScalarField
(
    const epsilonLowReRoughWallTPFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
	ks_(ewfpsf.ks_),
	sr_(ewfpsf.sr_),
	epsType_(ewfpsf.epsType_)
{
    checkType();
}


epsilonLowReRoughWallTPFvPatchScalarField::epsilonLowReRoughWallTPFvPatchScalarField
(
    const epsilonLowReRoughWallTPFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
	ks_(ewfpsf.ks_),
	sr_(ewfpsf.sr_),
	epsType_(ewfpsf.epsType_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonLowReRoughWallTPFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get patch indices
    const label patchI = patch().index();
	const fvMesh& mesh = patch().boundaryMesh().mesh();
    
	// Load Turbulence Model object
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");	
    const scalarField& y = rasModel.y()[patchI];
	
	const volScalarField& kr = db().lookupObject<volScalarField>("k");
	const volScalarField& kSqrtr = db().lookupObject<volScalarField>("kSqrt");
	const scalarField gradkSqrt = kSqrtr.boundaryField()[patchI].snGrad();
	
	const volScalarField& prodr = db().lookupObject<volScalarField>("tpProd");
	
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
	volScalarField phir = tpr*kr;
	scalarField gradphiSqrt = phir.boundaryField()[patchI].snGrad();
	
	const volVectorField& vort = db().lookupObject<volVectorField>("vorticity");
	const scalarField magVort = mag(vort.boundaryField()[patchI]);
	
	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
	const fvPatchScalarField& nutw = lookupPatchField<volScalarField, scalar>("nut");
	
	const fvPatchVectorField& Uw = lookupPatchField<volVectorField, vector>("U");
	const scalarField magGradUw = mag(Uw.snGrad());
    
   	const scalar Cmu25 = pow(Cmu_, 0.25);
	scalar epsC = 0.257;

    scalarField& epsw = refValue();
	
	
	if(epsType_ == "rough"){

		const scalar cB = 0.05;
		const scalar sM = 0.0025;
		const scalar cUp = 0.5;

    // Started with cUp = 0.26  
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
			scalar utauw = sqrt(nuw[faceI]*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/(nuw[faceI]);
			
			// Use epsilon constant region formula
			if(kPlus<=5.0){
				epsC = cUp;
			}else if(kPlus>5.0 && kPlus<=50.0){
				epsC = (cB+sM) + ((cUp-(cB+sM))/(0.73))*pow((1.0-kPlus/50.0),3.0);
			}else if(kPlus>50.0 && kPlus<=100.0){
				epsC = (cB+sM) - (sM/50.0)*(kPlus-50.0);
			}else{
				epsC = cB;
			}
			
			epsw[faceI] = epsC*pow((nuw[faceI]+nutw[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]);
		}
	}
	
	if(epsType_ == "roughtwo"){

		const scalar cB = 0.05;
		const scalar sM = 0.0025;
		const scalar cUp = 0.26;

    // Started with cUp = 0.26  
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
			scalar utauw = sqrt(nuw[faceI]*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/(nuw[faceI]);
			
			// Use epsilon constant region formula
			if(kPlus<=5.0){
				epsC = cUp;
			}else if(kPlus>5.0 && kPlus<=90.0){
				epsC = cUp + (kPlus/270.0);
			}else{
				epsC = 0.6;
			}
			
			epsw[faceI] = epsC*pow((nuw[faceI]+nutw[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]);
		}
	}
	
	if(epsType_ == "smooth"){
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
		
			epsw[faceI] = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
		}
		
	}
	
	if(epsType_ == "smoothphi"){
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
		
			epsw[faceI] = 2.0*nuw[faceI]*sqr(mag(gradphiSqrt[faceI]));
		}
		
	}
	
	if(epsType_ == "srfixed"){
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
		
			epsw[faceI] = sr_*pow((nuw[faceI]+nutw[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]);
		}
		
	}
	
	if(epsType_ == "pequal"){
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];		
			scalar utauw = sqrt(nuw[faceI]*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/(nuw[faceI]);		
		
			epsw[faceI] = 2.0*nuw[faceI]*sqr(mag(gradkSqrt[faceI])) + 0.786*prodr[faceCellI]*kr[faceCellI];
		}
		
	}

    fixedInternalValueFvPatchScalarField::updateCoeffs();
}


void epsilonLowReRoughWallTPFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchScalarField::evaluate(commsType);
}


void epsilonLowReRoughWallTPFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
	os.writeKeyword("ks") << ks_ << token::END_STATEMENT << nl;
	os.writeKeyword("sr") << sr_ << token::END_STATEMENT << nl;
	os.writeKeyword("epsType") << epsType_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonLowReRoughWallTPFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
