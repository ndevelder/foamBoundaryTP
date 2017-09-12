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

#include "fLowReRoughWallTPFvPatchScalarField.H"
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

void fLowReRoughWallTPFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("fLowReRoughWallTPFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar fLowReRoughWallTPFvPatchScalarField::calcYPlusLam
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



void fLowReRoughWallTPFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
	os.writeKeyword("ks") << ks_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fLowReRoughWallTPFvPatchScalarField::fLowReRoughWallTPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
	ks_(0.0),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


fLowReRoughWallTPFvPatchScalarField::fLowReRoughWallTPFvPatchScalarField
(
    const fLowReRoughWallTPFvPatchScalarField& ptf,
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
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


fLowReRoughWallTPFvPatchScalarField::fLowReRoughWallTPFvPatchScalarField
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
	ks_(dict.lookupOrDefault<scalar>("ks", 0.0)),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


fLowReRoughWallTPFvPatchScalarField::fLowReRoughWallTPFvPatchScalarField
(
    const fLowReRoughWallTPFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


fLowReRoughWallTPFvPatchScalarField::fLowReRoughWallTPFvPatchScalarField
(
    const fLowReRoughWallTPFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
	ks_(wfpsf.ks_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fLowReRoughWallTPFvPatchScalarField::updateCoeffs()
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
	dictionary tpCoeffDict(rasDictionary.subDict("turbulentPotentialCoeffs"));
	const scalar& sigmaKr = readScalar(tpCoeffDict.lookup("sigmaKInit")) ;
	
	const volScalarField& kr = mesh.lookupObject<volScalarField>("k");
	const volScalarField& tpr = mesh.lookupObject<volScalarField>("tpphi");
	const volScalarField& epsr = mesh.lookupObject<volScalarField>("epsilon");
	
	const volVectorField& vort = mesh.lookupObject<volVectorField>("vorticity");
	const scalarField magVort = mag(vort.boundaryField()[patchI].patchInternalField());
		
	const fvPatchScalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
	const fvPatchScalarField& nutw = lookupPatchField<volScalarField, scalar>("nut");
	
	const fvPatchVectorField& Uw = lookupPatchField<volVectorField, vector>("U");
	const scalarField magGradUw = mag(Uw.snGrad());
    
    const scalar Cmu25 = pow(Cmu_, 0.25);
	const scalar Cmu12 = pow(Cmu_, 0.5);
	
	tmp<scalarField> tfw(new scalarField(nutw.size()));
	scalarField& fw = tfw();
    
    forAll(nutw, faceI)
    {
		scalar nuEffw = nuw[faceI];
        label faceCellI = patch().faceCells()[faceI];		
		scalar utauw = sqrt(nuEffw*magGradUw[faceI]);
        scalar kPlus = ks_*utauw/nuEffw;
		
		//fw[faceI] = -20.0*sqr(nuw[faceI])*(kr[faceCellI])*tpr[faceCellI]/(epsr.boundaryField()[patchI][faceI]*sqr(sqr(y[faceI])));
		//fw[faceI] = -5.0*(epsr.boundaryField()[patchI][faceI]/kr[faceCellI])*tpr[faceCellI];
		//fw[faceI] = -2.0*nuw[faceI]*tpr[faceCellI]/sqr(y[faceI]);
		
		fw[faceI] = 0.0;
		
		//Info << fw[faceI] << endl;
		
    }
	
	operator==(fw);

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<scalarField> fLowReRoughWallTPFvPatchScalarField::yPlus() const
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

void fLowReRoughWallTPFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchScalarField::evaluate(commsType);
}


void fLowReRoughWallTPFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fLowReRoughWallTPFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
