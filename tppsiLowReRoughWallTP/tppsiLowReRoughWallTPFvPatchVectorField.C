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

#include "tppsiLowReRoughWallTPFvPatchVectorField.H"
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

void tppsiLowReRoughWallTPFvPatchVectorField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("tppsiLowReRoughWallTPFvPatchVectorField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar tppsiLowReRoughWallTPFvPatchVectorField::calcYPlusLam
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



void tppsiLowReRoughWallTPFvPatchVectorField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
	os.writeKeyword("ks") << ks_ << token::END_STATEMENT << nl;
	os.writeKeyword("cr") << cr_ << token::END_STATEMENT << nl;
	os.writeKeyword("bz") << bz_ << token::END_STATEMENT << nl;
	os.writeKeyword("pswType") << pswType_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tppsiLowReRoughWallTPFvPatchVectorField::tppsiLowReRoughWallTPFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
	ks_(0.0),
	cr_(1.0),
	bz_(0),
	pswType_("nut"),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


tppsiLowReRoughWallTPFvPatchVectorField::tppsiLowReRoughWallTPFvPatchVectorField
(
    const tppsiLowReRoughWallTPFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
	ks_(ptf.ks_),
	cr_(ptf.cr_),
	bz_(ptf.bz_),
	pswType_(ptf.pswType_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


tppsiLowReRoughWallTPFvPatchVectorField::tppsiLowReRoughWallTPFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
	ks_(dict.lookupOrDefault<scalar>("ks", 0.0)),
	cr_(dict.lookupOrDefault<scalar>("cr", 1.0)),
	bz_(dict.lookupOrDefault<label>("bz", 0)),
	pswType_(dict.lookupOrDefault<word>("pswType", "nut")),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


tppsiLowReRoughWallTPFvPatchVectorField::tppsiLowReRoughWallTPFvPatchVectorField
(
    const tppsiLowReRoughWallTPFvPatchVectorField& wfpsf
)
:
    fixedValueFvPatchVectorField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
	cr_(wfpsf.cr_),
	bz_(wfpsf.bz_),
	pswType_(wfpsf.pswType_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


tppsiLowReRoughWallTPFvPatchVectorField::tppsiLowReRoughWallTPFvPatchVectorField
(
    const tppsiLowReRoughWallTPFvPatchVectorField& wfpsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
	cr_(wfpsf.cr_),
	bz_(wfpsf.bz_),
	pswType_(wfpsf.pswType_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tppsiLowReRoughWallTPFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get patch indices
    const label patchI = patch().index();
	const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Get patch normal vector
    const fvBoundaryMesh& boundary = mesh.boundary();
    vectorField nf = boundary[patchI].nf();
	
    
	// Load Turbulence Model object
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");	
    const scalarField& y = rasModel.y()[patchI];
	
    const volScalarField& kr = db().lookupObject<volScalarField>("k");
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
	const volScalarField& epsr = db().lookupObject<volScalarField>("epsilon");
	
	const volVectorField& vort = db().lookupObject<volVectorField>("vorticity");
    const scalarField magVort = mag(vort.boundaryField()[patchI]);
    const vectorField unitVort = vort.boundaryField()[patchI]/magVort;

    //const volVectorField& uField = mesh.lookupObject<volVectorField>("U");
    //const surfaceTensorField du = fvc::interpolate(fvc::grad(uField));
    //tensorField dui = du.boundaryField()[patchI];
    //const vectorField GradUw = dui & nf;
	//const scalarField magGradUw = magVort;

	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
	const fvPatchScalarField& nutw = lookupPatchField<volScalarField, scalar>("nut");
	
	const fvPatchVectorField& Uw = lookupPatchField<volVectorField, vector>("U");
	const vectorField GradUw = Uw.snGrad();
    const scalarField magGradUw = mag(GradUw);
	

    const volVectorField& tppsiw = db().lookupObject<volVectorField>("tppsi");
	const scalarField magPsiw = mag(tppsiw.boundaryField()[patchI]*kr.boundaryField()[patchI]);
    const scalarField psiDV = magPsiw/(magGradUw + SMALL);

    //const fvPatchVectorField& tppsiw = *this;
	//const scalarField magPsiw = mag(tppsiw.patchInternalField()*kr.boundaryField()[patchI]);
    //const scalarField psiDV = magPsiw/(magGradUw + SMALL);
	
	tmp<vectorField> tpsw(new vectorField(Uw.size()));
	vectorField& psw = tpsw();
	
	vector psiwout = vector(0,0,0);

    forAll(nutw, faceI)
    {
		scalar nuEffw = nuw[faceI]+nutw[faceI];
		scalar nuPsiw = nuw[faceI]+psiDV[faceI];
			
		label faceCellI = patch().faceCells()[faceI];		
			
		scalar utauw = sqrt(nuEffw*magGradUw[faceI]);
        
		scalar kPlus = ks_*utauw/nuw[faceI];


		if(pswType_ == "nut"){
			if(kPlus <= 2.25){
				psw[faceI] = vector(0,0,0);
			}else{	    
				psw[faceI] = nutw[faceI]*vort.boundaryField()[patchI][faceI]/(kr.boundaryField()[patchI][faceI] + SMALL);
			}
		}

        if(pswType_ == "nutf"){
            if(kPlus <= 5.5){
                psw[faceI] = vector(0,0,0);
            }else{      
                psw[faceI] = cr_*nutw[faceI]*vort[faceCellI]/(kr.boundaryField()[patchI][faceI] + SMALL);
            }
        }

        if(pswType_ == "nutm"){
            if(kPlus <= 5.5){
                psw[faceI] = vector(0,0,0);
            }else{      
                psw[faceI] = 0.028*(kr.boundaryField()[patchI][faceI]/epsr.boundaryField()[patchI][faceI] + SMALL)*vort.boundaryField()[patchI][faceI];
            }
        }
		
		if(pswType_ == "nutplus"){	
            if(kPlus <= 2.25){
                psw[faceI] = vector(0,0,0);
            }else{      
    			psw[faceI] = (1.0 + cr_*(pow(1.0/kPlus, 0.333)))*(0.21*tpr.boundaryField()[patchI][faceI]*(kr.boundaryField()[patchI][faceI])/(epsr.boundaryField()[patchI][faceI] + SMALL))*vort.boundaryField()[patchI][faceI];
    		} 
        } 

		if(pswType_ == "dev"){
			if(kPlus<=5.5){
				psw[faceI] = vector(0,0,0);
			}else{	    
				psw[faceI] = min(1.0,(kPlus-5.0)/90.0)*(epsr.boundaryField()[patchI][faceI]/sqr(kr.boundaryField()[patchI][faceI]))*(nuPsiw*vort[faceCellI])/(0.3*magVort[faceI]);
			}
		}
		
		if(pswType_ == "itime"){
			if(kPlus <= 5.5){
				psw[faceI] = vector(0,0,0);
			}else{	    
				psw[faceI] = epsr[faceCellI]*vort[faceCellI]/(kr[faceCellI]*sqr(magGradUw[faceI]));
			}
		}
		
		if(pswType_ == "zero"){
			if(kPlus <= 5.5){
				psw[faceI] = vector(0,0,0);
			}else{	    
				psw[faceI] = vector(0,0,0);
			}
		}

        if(pswType_ == "uts"){
            if(kPlus <= 5.5){
                psw[faceI] = vector(0,0,0);
            }else{      
                psw[faceI] = cr_*nuEffw*vort.boundaryField()[patchI][faceI]/(kr.boundaryField()[patchI][faceI] + SMALL);
            }
        }

        if(pswType_ == "con"){
            if(kPlus <= 5.5){
                psw[faceI] = vector(0,0,0);
            }else{      
                psw[faceI] = (1.0 + cr_)*0.3*min(1.0,sqr(kPlus/90.0))*unitVort[faceI];
            }  
        }

        if(pswType_ == "tkeplusa"){
            if(kPlus < 2.25){
                psw[faceI] = vector(0,0,0);
            }else{  
                scalar tkeplus = max(1.0e-10,(1.0/0.3)*tanh(((log(kPlus/30.0)/log(10.0))+1.0-tanh(kPlus/125.0))*tanh(kPlus/125.0)));      
                psw[faceI] = cr_*(nutw[faceI]/(nuw[faceI]+nutw[faceI]))*(-GradUw[faceI]/magGradUw[faceI])/tkeplus;
            }
        }

        if(pswType_ == "tkeplusk"){
            if(kPlus < 2.25){
                psw[faceI] = vector(0,0,0);
            }else{  
                scalar tkeplus = max(1.0e-10,(1.0/0.3)*min(1.0,kPlus/90.0));      
                psw[faceI] = cr_*(nutw[faceI]/(nuw[faceI]+nutw[faceI]))*(-GradUw[faceI]/magGradUw[faceI])/tkeplus;
            }
        }
		
		psiwout = psw[faceI];
		
                //if(patch().name() == "FOIL_TOP"){
        //Pout << "Psw: " << psw[faceI]  << endl;
                //}

		//if(patch().name() == "FOIL_TOP"){
		//   Pout << "grad Uw: " << GradUw[faceI] << "  Vorticity: " << vort.boundaryField()[patchI][faceI]  << endl;
		//}
    }

    //Pout<< "Psi w: " << psw << endl;
    //Pout<< "magVort w: " << magVort << endl;
    //Pout<< "Ugrad Mag w: " << magGradUw << endl;
    //Pout<< "nf w: " << nf << endl;

    //Info<< "Tppsi w: " << psiwout << endl;

	
	//if(patch().name() == "FOIL_TOP"){
	//	Pout << "min nut W: " << min(nutw) << endl;
	//}  

	vectorField::operator=(psw);
	
    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<scalarField> tppsiLowReRoughWallTPFvPatchVectorField::yPlus() const
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

void tppsiLowReRoughWallTPFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchVectorField::evaluate(commsType);
}

void tppsiLowReRoughWallTPFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, tppsiLowReRoughWallTPFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
