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
	bz_(0),
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
	bz_(ptf.bz_),
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
	sr_(dict.lookupOrDefault<scalar>("sr", 1.0)),
	bz_(dict.lookupOrDefault<scalar>("bz", 0)),
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
	bz_(ewfpsf.bz_),
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
	bz_(ewfpsf.bz_),
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
	
	const dictionary& rasDictionary = db().lookupObject<IOdictionary>("RASProperties");
	dictionary bCoeffDict(rasDictionary.subDict("boundaryCoeffs"));
	const scalar& cMuw = readScalar(bCoeffDict.lookup("cMu"));
	const scalar& betaKw = readScalar(bCoeffDict.lookup("betaK"));
	
	const volScalarField& kr = db().lookupObject<volScalarField>("k");
	const volScalarField& ksr = db().lookupObject<volScalarField>("kSqrt");

	scalarField gradkSqrt = ksr.boundaryField()[patchI].snGrad();
	
	const volScalarField& prodr = db().lookupObject<volScalarField>("tpProd");
	const volScalarField& epsr = db().lookupObject<volScalarField>("epsilon");
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
	volScalarField phir = tpr*kr;
	volScalarField phirsqrt = sqrt(phir);
	scalarField gradphiSqrt = phir.boundaryField()[patchI].snGrad();
	
	const volVectorField& vort = db().lookupObject<volVectorField>("vorticity");
	const scalarField magVort = mag(vort.boundaryField()[patchI]);
	
	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
	const fvPatchScalarField& nutw = lookupPatchField<volScalarField, scalar>("nut");
	
	
	const fvPatchVectorField& Uw = lookupPatchField<volVectorField, vector>("U");
	const scalarField magGradUw = mag(Uw.snGrad());
	
	const volVectorField& tppsiw =  db().lookupObject<volVectorField>("tppsi");
    volVectorField Psiw = (tppsiw*kr);
	const scalarField magPsiw = mag(Psiw.boundaryField()[patchI]);
    const scalarField psiDV = magPsiw/(magGradUw + SMALL);

    const scalarField pde = (Psiw.boundaryField()[patchI] & Psiw.boundaryField()[patchI])/(cMuw*phir.boundaryField()[patchI]*kr.boundaryField()[patchI]);
	
   	const scalar Cmu25 = pow(Cmu_, 0.25);
   	const scalar bk25 = pow(0.09,0.25);
	scalar epsC = 0.257;

    scalarField& epsw = refValue();
	
	volScalarField& GdKw = const_cast<volScalarField&>(db().lookupObject<volScalarField>("tpProd"));
	
	scalar epW = 0.0;
	scalar epWS = 0.0;
	scalar epWR = 0.0;
	scalar kpout = 0.0;
	scalar pdeout = 0.0;
	
	
	if(epsType_ == "hellsten"){  

        //Info << "Using Hellsten Epsilon BC" << endl;
		
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			//label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/nuw[faceI];	
			
			
			// Use epsilon constant region formula
			if(kPlus<=5.0){
				epsw[faceI] = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			}else if(kPlus>5.0 && kPlus<=25.0){
				epsC = pow(50.0/kPlus,2.0);
				epsw[faceI] = 0.09*kr.boundaryField()[patchI][faceI]*epsC*(nuPsiw*magGradUw[faceI])/(nuw[faceI]);
			}else{
				epsC = max(100.0/kPlus,1.0);
				epsw[faceI] = 0.09*kr.boundaryField()[patchI][faceI]*epsC*(nuPsiw*magGradUw[faceI])/(nuw[faceI]);
			}
			
			//Info << epsw[faceI] << "--" << kr.boundaryField()[patchI][faceI]<< "--"  << epsC << "--" << ((nuw[faceI]+nutw[faceI])*magGradUw[faceI]) << endl;			
		}
	}
	
	

	if(epsType_ == "knoppthree" || epsType_ == "knopp" || epsType_ == "rough"){
		
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			//label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/nuw[faceI];	

			label psize = nutw.size();

			if(faceI < bz_){
				kPlus = kPlus*exp(-1.0*((bz_-1)-faceI));
			}

			if(faceI > ((psize-1)-bz_)){
				kPlus = kPlus*exp(-1.0*((bz_-1)-((psize-1)-faceI)));
			}
		
			
			scalar dzm = min(1,pow(kPlus/30.0,0.67))*min(1,pow(kPlus/45.0,0.25))*min(1,pow(kPlus/60.0,0.25));
			scalar dz = dzm*0.03*ks_;

			
			// Use epsilon constant region formula
			if(kPlus<=5.5){
				epsw[faceI] = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			}else{
				epsw[faceI] = pow(min(1,kPlus/90.0),1.5)*0.3*pow((nuw[faceI]+psiDV[faceI])*magGradUw[faceI],1.5)/(0.41*dz + SMALL);
				//epsw[faceI] = 0.0493*pow(kr[faceCellI],1.5)/(0.41*dz + SMALL);
			}

			pdeout = pde[faceI];
			epW = epsw[faceI];
			epWS = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			kpout = kPlus;
			
		}
		
		//Info << "epw: " << epW << " epwS: " << epWS << " pde: " << pdeout << " kP: " << kpout  <<endl;
	}


	if(epsType_ == "smooth"){
		
		forAll(nutw, faceI)
		{
			//label faceCellI = patch().faceCells()[faceI];
		
			epsw[faceI] = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
		}
		
	}

	
	if(epsType_ == "smoothdelta"){
		
		forAll(nutw, faceI)
		{
			label faceCellI = patch().faceCells()[faceI];
		
			epsw[faceI] = 2.0*nuw[faceI]*kr[faceCellI]/sqr(y[faceI]);
		}
		
	}

	
	if(epsType_ == "utau4"){
		
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			//label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/nuw[faceI];		
			
			//scalar epsCalc = min(2.0,(0.95 + kPlus/90.0))*0.229*pow((nuw[faceI]+nutw[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]);
			
			//scalar epsMult = min( 0.229*(1.0 + pow(kPlus/118.0, 1.5)) , 0.5);
			
			scalar epsMult = min((kPlus/90.0),1.0)*(pow(min(1.0,kPlus/90.0), 0.04));
						
			if(kPlus<=5.5){
				epsw[faceI] = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			}else{
				epsw[faceI] = 0.31*epsMult*(0.8 + 0.2*pde[faceI])*pow((nuw[faceI]+psiDV[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]+psiDV[faceI]);
			}
			
			epW = epsw[faceI];
			epWS = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			epWR = 0.31*epsMult*pow((nuw[faceI]+psiDV[faceI])*magGradUw[faceI],2.0)/(nuw[faceI]+psiDV[faceI]);
			
		}
		
		Info << "epw: " << epW << " epwS: " << epWS << " epwR: " << epWR << endl;
		
	}
	
	
	scalar tOne = 0.0;
	scalar tTwo = 0.0;


	if(epsType_ == "utau3"){
		
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			//label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/nuw[faceI];			
			
			scalar utausqr = sqr(utauw);
			
			scalar epsMult = 0.23*min(pow(5.0/kPlus,1.0),1.0)*utausqr + magPsiw[faceI];
			
			scalar epsCalc = epsMult*utausqr/nuPsiw;
			
			if(kPlus<=5.0){
				epsw[faceI] = min(2.0*nuw[faceI]*sqr(gradkSqrt[faceI]),epsCalc);
			}else{
				epsw[faceI] = epsCalc;
			}
		}
		
	}
	
	scalar eMp = 0.0;
	scalar eMu = 0.0;
	scalar kpo = 0.0;
		
	
	if(epsType_ == "psi4"){
		
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar utauvw = sqrt(nuw[faceI]*magGradUw[faceI]); 
			scalar kPlus = ks_*utauw/nuw[faceI];
			
			scalar epsMult = 1.0 + 0.21*pow(min(15.0/kPlus,1.0),2.0);
			
			if(kPlus<=5.5){
				epsw[faceI] = 2.0*nuw[faceI]*kr[faceCellI]/sqr(y[faceI]);
			}else{
				epsw[faceI] = epsMult*sqr((nuPsiw*magGradUw[faceI]))/nuPsiw;
			}
			
			epW = epsw[faceI];
			epWS = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			epWR = epsMult*nuPsiw*magGradUw[faceI];
			eMp = magPsiw[faceI];
			eMu = 0.23*nuPsiw*magGradUw[faceI]*pow(min(25.0/kPlus,1.0),2.0);
			kpo = kPlus;
			
		}
		
		Info << "epw: " << epW << " epwS: " << epWS << " epwR: " << epWR << " eMp: " << eMp << " eMu: " << eMu << " kpo: " << kpo << endl;
	}

	scalar kw=0;
	scalar ksw=0.0;
   

	if(epsType_ == "psi5"){
			
		forAll(nutw, faceI)
		{
			//scalar nuEffw = nuw[faceI]+nutw[faceI];
			scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
			label faceCellI = patch().faceCells()[faceI];		
			
			scalar utauw = sqrt(nuPsiw*magGradUw[faceI]);
			scalar kPlus = ks_*utauw/nuw[faceI];
			
			scalar epsMult = magPsiw[faceI]/nuPsiw + 0.23*bk25*ksr.boundaryField()[patchI][faceI]*pow(min(25.0/kPlus,1.0),2.0)/nuw[faceI];
			
			if(kPlus<=5.5){
				epsw[faceI] = 2.0*nuw[faceI]*kr[faceCellI]/sqr(y[faceI]);
			}else{
				epsw[faceI] = epsMult*bk25*ksr.boundaryField()[patchI][faceI];
			}
			
			epW = epsw[faceI];
			epWS = 2.0*nuw[faceI]*sqr(gradkSqrt[faceI]);
			epWR = epsMult*nuPsiw*magGradUw[faceI];
			eMp = magPsiw[faceI];
			eMu = 0.23*nuPsiw*magGradUw[faceI]*pow(min(25.0/kPlus,1.0),2.0);
			kw = kr.boundaryField()[patchI][faceI];
			ksw = ksr.boundaryField()[patchI][faceI];
			
		}
		
		Info << "epw: " << epW << " epwS: " << epWS << " kw: " << kw<< " magpsiw: " << eMp << " ksqrtw: " << ksw << endl;
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
	os.writeKeyword("bz") << bz_ << token::END_STATEMENT << nl;
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
