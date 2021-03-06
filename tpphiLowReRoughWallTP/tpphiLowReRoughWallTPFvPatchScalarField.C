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

#include "tpphiLowReRoughWallTPFvPatchScalarField.H"
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

void tpphiLowReRoughWallTPFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("tpphiLowReRoughWallTPFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar tpphiLowReRoughWallTPFvPatchScalarField::calcYPlusLam
(
    const scalar kappa,
    const scalar E
) const
{
    scalar ypl = 11.0;

    for (int i = 0; i < 10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}



void tpphiLowReRoughWallTPFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
	os.writeKeyword("ks") << ks_ << token::END_STATEMENT << nl;
	os.writeKeyword("pkC") << pkC_ << token::END_STATEMENT << nl;
	os.writeKeyword("bz") << bz_ << token::END_STATEMENT << nl;
	os.writeKeyword("rType") << rType_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tpphiLowReRoughWallTPFvPatchScalarField::tpphiLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
	ks_(1e-10),
	pkC_(0.1455),
	bz_(0),
    yPlusLam_(calcYPlusLam(kappa_, E_)),
	rType_("calculated")
{
    checkType();
}


tpphiLowReRoughWallTPFvPatchScalarField::tpphiLowReRoughWallTPFvPatchScalarField
(
    const tpphiLowReRoughWallTPFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
	ks_(ptf.ks_),
	pkC_(ptf.pkC_),
	bz_(ptf.bz_),
    yPlusLam_(ptf.yPlusLam_),
	rType_(ptf.rType_)
{
    checkType();
}


tpphiLowReRoughWallTPFvPatchScalarField::tpphiLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
	ks_(dict.lookupOrDefault<scalar>("ks", 1e-10)),
	pkC_(dict.lookupOrDefault<scalar>("pkC", 0.1333)),
	bz_(dict.lookupOrDefault<label>("bz", 0)),
    yPlusLam_(calcYPlusLam(kappa_, E_)),
	rType_(dict.lookupOrDefault<word>("rType","calculated"))
{
    checkType();
}


tpphiLowReRoughWallTPFvPatchScalarField::tpphiLowReRoughWallTPFvPatchScalarField
(
    const tpphiLowReRoughWallTPFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
	pkC_(wfpsf.pkC_),
	bz_(wfpsf.bz_),
    yPlusLam_(wfpsf.yPlusLam_),
	rType_(wfpsf.rType_)
{
    checkType();
}


tpphiLowReRoughWallTPFvPatchScalarField::tpphiLowReRoughWallTPFvPatchScalarField
(
    const tpphiLowReRoughWallTPFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
	pkC_(wfpsf.pkC_),
	bz_(wfpsf.bz_),
    yPlusLam_(wfpsf.yPlusLam_),
	rType_(wfpsf.rType_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tpphiLowReRoughWallTPFvPatchScalarField::updateCoeffs()
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
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
	
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
	
    
    tmp<scalarField> tphw(new scalarField(nutw.size()));
	scalarField& phw = tphw();
	
	tmp<scalarField> tkplus(new scalarField(nutw.size()));
	scalarField& kPlus = tkplus();
	
	scalar pkF = SMALL;
	scalar pkout = SMALL;

    forAll(nutw, faceI)
    {
		scalar nuEffw = nuw[faceI]+nutw[faceI];
		scalar nuPsiw = nuw[faceI]+psiDV[faceI];
		
        label faceCellI = patch().faceCells()[faceI];		
		
		//scalar utauw = pow(magPsiSqrw[faceI] + SMALL, 0.25);
		scalar utauw = sqrt(nuEffw*magGradUw[faceI] + SMALL);

        scalar kPlus = ks_*utauw/nuw[faceI];
        scalar kP = ks_*utauw/nuw[faceI];
        scalar dzero = 0.03*ks_*min(1,pow(kPlus/30.0,0.67))*min(1,pow(kPlus/45.0,0.25))*min(1,pow(kPlus/60.0,0.25));


		
		if(rType_ == "channel") {
			
			phw[faceI] = pow((1.0/5.5)*log(ks_/nuw[faceI]) - (3.0*kappa_/5.5),2.0)/kr.boundaryField()[patchI][faceI];
			//Pout << "phi W: " << pow((1.0/5.5)*log(ks_/nuw[faceI]) - (3.0*kappa_/5.5),2.0) << endl;
			//Pout << "k W: " << kr.boundaryField()[patchI][faceI] << endl;
		
		}else if(rType_ == "calculated"){
			
			if(kPlus<30.0){
			  pkF = 2.0*pow((kPlus/90.0),2.0);
			}else if(kPlus<90.0){
			  pkF = 0.18 + (0.088/70.0)*kPlus;
			}else{
			  pkF = 0.19 + (0.09/70.0)*100.0;
			}
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "calcratio"){
			
			if(kPlus<5.5){
			  pkF = 1e-10;
			}else if(kPlus<45.0){
			  pkF = 0.12*pow(kPlus/45.0,1.5) - 0.0014*(1.0 - (kPlus/45.0));
			}else if(kPlus<90.0){
			  pkF = 0.015*(kPlus - 45.0)/45.0 + 0.12;
			}else{
			  pkF = 0.135;
			}
			
			if(pkF < 0.0){
				Info << "Negative phiw" << endl;
			}
			
			phw[faceI] = pkF;

		}else if(rType_ == "powclean"){

			pkF = pkC_;
			phw[faceI] = pkF;
			
		}else if(rType_ == "pow"){
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else if(kPlus<=90.0){
			  pkF = pkC_*pow(min(1.0,(kPlus-2.25)/90.0),0.04);
			}else{
			  pkF = pkC_;
			}
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "sqrtpow"){
			
			if(kPlus<=2.25){
			  pkF = 1e-5;
			}else if(kPlus<=90.0){
			  pkF = pkC_*pow((kPlus-2.25)/90.0,0.02);
			}else{
			  pkF = pkC_;
			}
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "tanh"){
			
			if(kPlus<=5.5){
			  pkF = 1e-10;
			}else if(kPlus<=90.0){
			  pkF = pkC_*((0.04*kPlus/90.0) + (tanh((kPlus-5.5)/5.0)));
			}else{
			  pkF = pkC_;
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "aupoix"){
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else{
			  //pkF = max(1e-10, pkC_*(tanh((kPlus-2.25)/15.0) + 0.0*((kPlus-2.25)/1180.0)*(1.0-exp(-((kPlus-2.25)/100.0)))));
			  pkF = max(1e-10, pkC_*(tanh((kPlus-2.25)/15.0) + 0.0*tanh((kPlus-2.25)/800.0)));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "aupoixs"){
			
			if(kPlus<=2.25){
			  pkF = 1e-5;
			}else{
			  //pkF = sqrt(max(1e-10, pkC_*(tanh((kPlus-2.25)/15.0) + 0.0*((kPlus-2.25)/1180.0)*(1.0-exp(-((kPlus-2.25)/100.0))))));
			  pkF = sqrt(max(1e-10, pkC_*(tanh((kPlus-2.25)/20.0)*tanh((kPlus-2.25)/7.0) + 0.0*tanh((kPlus-2.25)/800.0))));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "aupoix2"){
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else{
			  //pkF = max(1e-10, pkC_*(tanh((kPlus-2.25)/15.0) + 0.0*((kPlus-2.25)/1180.0)*(1.0-exp(-((kPlus-2.25)/100.0)))));
			  pkF = max(1e-10, pkC_*pow(tanh((kPlus-2.25)/15.0),3.0));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "aupoixs2"){
			
			if(kPlus<=2.25){
			  pkF = 1e-5;
			}else{
			  //pkF = sqrt(max(1e-10, pkC_*(tanh((kPlus-2.25)/15.0) + 0.0*((kPlus-2.25)/1180.0)*(1.0-exp(-((kPlus-2.25)/100.0))))));
			  pkF = sqrt(max(1e-10, pkC_*pow(tanh((kPlus-2.25)/15.0),3.0 )));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "knopp"){
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else if(kPlus > 90.0){
			  pkF = max(1e-10, pkC_*(tanh((kPlus-2.25)/20.0) - 0.0*(0.8-exp(-(kPlus/300.0)))) + 0.0*(1.0-exp(-(kPlus-90.0)/650.0)));
			}else{
			  pkF = max(1e-10, pkC_*(tanh((kPlus-2.25)/20.0) - 0.0*(0.8-exp(-(kPlus/300.0)))));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "knopps"){  
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else if(kPlus > 90.0){
			  pkF = sqrt(max(1e-10, pkC_*(tanh((kPlus-2.25)/20.0) - 0.0*(0.8-exp(-(kPlus/300.0)))) + 0.0*(1.0-exp(-(kPlus-90.0)/650.0))));
			}else{
			  pkF = sqrt(max(1e-10, pkC_*(tanh((kPlus-2.25)/20.0) - 0.0*(0.8-exp(-(kPlus/300.0))))));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "knoppn"){
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else if(kPlus > 90.0){
			  pkF = max(1e-10, pkC_*pow(tanh((kPlus-2.25)/14.0),3.0));
			}else{
			  pkF = max(1e-10, pkC_*pow(tanh((kPlus-2.25)/14.0),3.0));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "knoppns"){  
			
			if(kPlus<=2.25){
			  pkF = 1e-10;
			}else if(kPlus > 90.0){
			  pkF = sqrt(max(1e-10, pkC_*pow(tanh((kPlus-2.25)/14.0),3.0)));
			}else{
			  pkF = sqrt(max(1e-10, pkC_*pow(tanh((kPlus-2.25)/14.0),3.0)));
			} 
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "develder"){ 
			
			if(kPlus<5.001){
			  pkF = 1e-10;
			}else if(kPlus<95.0){
			  pkF = min(0.12+0.015*pow((kPlus)/90.0,0.02),pkC_);
			}else{
			  pkF = pkC_;
			}
			
			phw[faceI] = pkF;
			
		}else if(rType_ == "fixedratio"){ 
		
			if(kPlus < 2.0){
			  pkF = (1e-10)*(1.0 - kPlus/2.0) + (kPlus/2.0)*pkC_;
			}else{
			  pkF = pkC_;
			}

			phw[faceI] = pkF;

		}else if(rType_ == "fixed"){ 
		
			phw[faceI] = pkC_/kr.boundaryField()[patchI][faceI];
			
		}else if(rType_ == "smooth"){ 
		
			phw[faceI] = SMALL;
			
		}else{ 
		
			phw[faceI] = 0.1455;
			
		}
    
		//if(patch().name() == "WALL_TOP"){
		//	Pout << "kPlus: " << kP << endl;
		//	Pout << "ks: " << ks_ << endl;
		//	Pout << "utauw: " << utauw << endl;
		//	Pout << "dzero: " << dzero << endl;
		//}

		pkout = phw[faceI];
	}

	//Pout<< "Phi/k w: " << pkout << endl;
	
	operator==(tphw);
 
    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<scalarField> tpphiLowReRoughWallTPFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];

    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();
    const scalarField kwc = k.boundaryField()[patchI].patchInternalField();
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    return pow(Cmu_, 0.25)*y*sqrt(kwc)/nuw;
}

void tpphiLowReRoughWallTPFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchScalarField::evaluate(commsType);
}

void tpphiLowReRoughWallTPFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, tpphiLowReRoughWallTPFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
