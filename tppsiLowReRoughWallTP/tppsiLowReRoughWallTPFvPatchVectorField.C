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
#include "turbulentPotential.H"
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
    
	// Load Turbulence Model object
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");	
    const scalarField& y = rasModel.y()[patchI];
	
	const volScalarField& nuw = db().lookupObject<volScalarField>("nu");
	const volScalarField& nutw = db().lookupObject<volScalarField>("nut");
	const volScalarField& kr = db().lookupObject<volScalarField>("k");
	const volScalarField& tpr = db().lookupObject<volScalarField>("tpphi");
	
	const volVectorField& vort = db().lookupObject<volVectorField>("vorticity");
	const scalarField magVort = mag(vort.boundaryField()[patchI]);
	
	const fvPatchVectorField& Uw = lookupPatchField<volVectorField, vector>("U");
	const scalarField magGradUw = mag(Uw.snGrad());
	
    vectorField& tppsiw = *this;

    forAll(Uw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];		
		scalar utauw = sqrt(nuw.boundaryField()[patchI][faceI]*magGradUw[faceI]);
        scalar kPlus = ks_*utauw/nuw.boundaryField()[patchI][faceI];
        tppsiw[faceI] = nutw.boundaryField()[patchI][faceI]*vort.boundaryField()[patchI][faceI]/kr.boundaryField()[patchI][faceI];
    }
	
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


void tppsiLowReRoughWallTPFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
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
