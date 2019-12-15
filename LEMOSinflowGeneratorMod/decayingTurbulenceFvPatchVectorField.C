/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "decayingTurbulenceFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "OFstream.H"
#include "IFstream.H"
#include <ctime>
#include "collatedFileOperation.H"
#include "uncollatedFileOperation.H"

Foam::Random ranGen(Foam::label(time(0)));

const Foam::scalar OVERLAP = 0.3;

const Foam::label NUM_VORT = 10;

inline Foam::scalar lspot(Foam::scalar l)
{
    return 3*l;
}


Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    LField_(p.size()),
    refField_(p.size()),
    RField_(p.size()),
    curTimeIndex_(-1),
    vortons_(),
    Lund_(),
    R_(),
    ind_(2),
    direction_(1),
    perturb_(0),
    nVortField_(p.size()),
    inletShape_(1),
    Umean_(0,0,0),
    width_(0),
    midRadius_(0),
    center_(0,0,0),
    Radius_(0),
    n_(0)
{}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    LField_(mapper(ptf.LField_)),
    refField_(mapper(ptf.refField_)),
    RField_(mapper(ptf.RField_)),
    curTimeIndex_(-1),
    vortons_(ptf.vortons_),
    Lund_(mapper(ptf.Lund_)),
    R_(mapper(ptf.R_)),
    ind_(ptf.ind_),
    direction_(ptf.direction_),
    perturb_(ptf.perturb_),
    nVortField_(mapper(ptf.nVortField_)),
    inletShape_(ptf.inletShape_),
    Umean_(ptf.Umean_),
    width_(ptf.width_),
    midRadius_(ptf.midRadius_),
    center_(ptf.center_),
    Radius_(ptf.Radius_),
    n_(ptf.n_)
{}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    LField_("LField", dict, p.size()),
    refField_(p.size(),pTraits<vector>::zero),
    RField_("RField", dict, p.size()),
    curTimeIndex_(-1),
    Lund_(p.size(), pTraits<tensor>::zero),
    R_(RField_),
    ind_(2),
    direction_(dict.lookupOrDefault("direction", label(1))),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    nVortField_(p.size(),pTraits<scalar>::zero),
    inletShape_(dict.lookupOrDefault("inletShape", label(1))),
    Umean_(0,0,0),
    width_(0),
    midRadius_(0),
    center_(0,0,0),
    Radius_(0),
    n_(0)
{
    Info << "dict constructor for " << patch().name() << endl;
    /*
    if (Pstream::master())
    {
        fileName inputFile;

        const fileOperation& fp = Foam::fileHandler();
        if (isA<Foam::fileOperations::collatedFileOperation>(fp))
        {
            inputFile = db().time().rootPath()/db().time().globalCaseName()/"processors"+Foam::name(Pstream::nProcs())/db().time().timeName()/+"vortons_"+patch().name();
        }
        else if (isA<Foam::fileOperations::uncollatedFileOperation>(fp))
        {
            inputFile = db().time().path()/db().time().timeName()/+"vortons_"+patch().name();
        }

        IFstream os2(inputFile);

        if (os2.opened())
        {
            word name;
            os2.read(name);
            if (name=="vortons")
            {
                Info << "found vortons file: "<<name<< endl;
                os2 >> vortons_;
            }
        }

        Info << "vortons_.size()===" << vortons_.size() << endl;
    }
    */

    // 3. read in the boundary part of U
    if (dict.found("vortons"))
        vortons_ = SLList<decayingVorton>(dict.lookup("vortons"));
    
    Info << "vortons_.size()===" << vortons_.size() << endl;


    const vectorField& Cf = patch().Cf();

    if (inletShape_==1)
    {
        Umean_ = dict.lookup("Umean");
        width_ = readScalar(dict.lookup("width"));
        midRadius_ = readScalar(dict.lookup("midRadius"));
        center_ = dict.lookup("center");
        n_ = readScalar(dict.lookup("n"));
        scalar delta_ = width_/2;

        Info<<"width_======"<<width_<<endl;
        Info<<"delta_======"<<delta_<<endl;

        refField_ = 1.218*Umean_*pow( 1. - mag(mag(Cf - center_) - midRadius_)/1.01/delta_, n_);//for ring
    }
    else if (inletShape_==2)
    {
        Umean_ = dict.lookup("Umean");
        center_ = dict.lookup("center");
        Radius_ = readScalar(dict.lookup("Radius"));
        n_ = readScalar(dict.lookup("n"));
        vector Umax = Umean_/2./(7./8.-7./15.);

        refField_ = Umax*pow(1. - mag(Cf - center_)/Radius_, n_);
    }
    else
    {
        FatalErrorInFunction
            << "You can shoose 'ring', 'tube', instead of inletShape_!" << endl;
    }

    R_= RField_;
    if (dict.found("R"))
        R_ = Field<symmTensor>("R", dict, p.size());

    if (dict.found("value"))
        fixedValueFvPatchField<vector>::operator==(Field<vector>("value", dict, p.size()));
    else
        fixedValueFvPatchField<vector>::operator==(refField_);

    forAll(RField_, I)
    {
        if ((RField_[I].xx() < 0) || (RField_[I].yy() < 0) || (RField_[I].zz() < 0))
        FatalErrorIn("decayingTurbulenceFvPatchVectorField::")<< "Some RMS are negative in R"<<abort(FatalError);
    }

    bool bError = false;
    bool bRobustLund = dict.lookupOrDefault("robustLundTransformation",false);

    forAll(RField_,i)
    {
        if (RField_[i].component(symmTensor::XX)<=0.0)
        {
            Info << "Found a negative or 0 Re Stress in XX" << endl;
            bError=true;
            Lund_[i]=pTraits<tensor>::I;
        }
    }
    if (!bRobustLund and !bError)
    {
        Info << "Standard Lund transformation algorithm applied." << endl;
        Lund_.replace(tensor::XX, sqrt(RField_.component(symmTensor::XX)));
        Lund_.replace(tensor::YX, RField_.component(symmTensor::XY)/Lund_.component(tensor::XX));
        Lund_.replace(tensor::ZX, RField_.component(symmTensor::XZ)/Lund_.component(tensor::XX));
        Lund_.replace(tensor::YY, sqrt(RField_.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX))));
        Lund_.replace(tensor::ZY, (RField_.component(symmTensor::YZ) - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY));
        Lund_.replace(tensor::ZZ, sqrt(RField_.component(symmTensor::ZZ) - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY))));
    }
    else if (bRobustLund and !bError)
    {
        Info << "Robust Lund transformation algorithm applied." << endl;
        forAll(Lund_,i)
        {
            bool errorFound = false;
            scalar a11 = sqrt(RField_[i][0]);
            scalar a21 = RField_[i][1]/a11;
            scalar a31 = RField_[i][2]/a11;

            if (RField_[i][3]<=sqr(a21))  // test if a22 exists (risk of sqrt a negative number)
            {
                errorFound = true;
            };
            if (errorFound==false) // test if a32 exists (a22 must be different than zero)
            {
                scalar a22 = sqrt(RField_[i][3]-sqr(a21));
                if (a22==0.0)
                {
                    errorFound = true;
                };
            };
            if (errorFound==false) // test if a33 exists (risk of sqrt a negative number)
            {
                scalar a22 = sqrt(RField_[i][3]-sqr(a21));
                scalar a32 = (RField_[i][4]-a21*a31)/a22;
                if (RField_[i][5]<(sqr(a31)+sqr(a32)))
                {
                    errorFound = true;
                };
            };

            if (errorFound==true)
            {
                Lund_[i][0] = sqrt(RField_[i][0]);  //xx
                Lund_[i][4] = sqrt(RField_[i][3]);  //yy
                Lund_[i][8] = sqrt(RField_[i][5]);  //zz
            }
            else
            {
                Lund_[i][0] = sqrt(RField_[i][0]);  //xx
                Lund_[i][3] = RField_[i][1]/Lund_[i][0];  //yx
                Lund_[i][4] = sqrt(RField_[i][3]-sqr(Lund_[i][3]));  //yy
                Lund_[i][6] = RField_[i][2]/Lund_[i][0];  //zx
                Lund_[i][7] = (RField_[i][4]-Lund_[i][3]*Lund_[i][6])/Lund_[i][4];  //zy
                Lund_[i][8] = sqrt(RField_[i][5]-sqr(Lund_[i][6])-sqr(Lund_[i][7]));  //zz
            };

        };
    }else
    {
        Lund_=pTraits<tensor>::I;
    }

    if (dict.found("ind"))
        ind_ = readScalar(dict.lookup("ind"));

    const scalarField& sf = patch().magSf();

    int numsf = 0;

    forAll(sf, I)
    {
        if (sf[I] > LField_[I]*LField_[I])
        numsf++;
    }

    if (numsf > 0)
        Pout<<"Warning: there are "<<numsf<<" inlet cells (out of "<<patch().size()<<") where the integral length is smaller than the cell size."<<endl;
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    LField_(ptf.LField_),
    refField_(ptf.refField_),
    RField_(ptf.RField_),
    curTimeIndex_(-1),
    vortons_(ptf.vortons_),
    Lund_(ptf.Lund_),
    R_(ptf.R_),
    ind_(ptf.ind_),
    direction_(ptf.direction_),
    perturb_(ptf.perturb_),
    nVortField_(ptf.nVortField_),
    inletShape_(ptf.inletShape_),
    Umean_(ptf.Umean_),
    width_(ptf.width_),
    midRadius_(ptf.midRadius_),
    center_(ptf.center_),
    n_(ptf.n_)
{
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    LField_(ptf.LField_),
    refField_(ptf.refField_),
    RField_(ptf.RField_),
    curTimeIndex_(-1),
    vortons_(ptf.vortons_),
    Lund_(ptf.Lund_),
    R_(ptf.R_),
    ind_(ptf.ind_),
    direction_(ptf.direction_),
    perturb_(ptf.perturb_),
    nVortField_(ptf.nVortField_),
    inletShape_(ptf.inletShape_),
    Umean_(ptf.Umean_),
    width_(ptf.width_),
    midRadius_(ptf.midRadius_),
    center_(ptf.center_),
    n_(ptf.n_)
{
}

void Foam::decayingTurbulenceFvPatchVectorField::autoMap(const fvPatchFieldMapper& m)
{
    fixedValueFvPatchField<vector>::autoMap(m);
    m(LField_, LField_);
    m(refField_, refField_);
    m(RField_, RField_);
    m(R_, R_);
    m(nVortField_, nVortField_);
}

void Foam::decayingTurbulenceFvPatchVectorField::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const decayingTurbulenceFvPatchVectorField& tiptf = refCast<const decayingTurbulenceFvPatchVectorField>(ptf);

    LField_.rmap(tiptf.LField_, addr);
    refField_.rmap(tiptf.refField_, addr);
    RField_.rmap(tiptf.RField_, addr);
    R_.rmap(tiptf.R_, addr);
    nVortField_.rmap(tiptf.nVortField_, addr);
}


void Foam::decayingTurbulenceFvPatchVectorField::updateCoeffs()
{
    if (this->updated()) return;

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        doUpdate();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::decayingTurbulenceFvPatchVectorField::readInletData()
{
}

void Foam::decayingTurbulenceFvPatchVectorField::writeSourceFields()
{
}

void Foam::decayingTurbulenceFvPatchVectorField::doUpdate()
{
    Field<vector>& patchField = *this;

    Field<vector> tloc = patch().Cf();

    List<List<scalar> > l(Pstream::nProcs());
    List<List<vector> > cf(Pstream::nProcs());
    List<List<vector> > rf(Pstream::nProcs());

    l[Pstream::myProcNo()].setSize(LField_.size());
    cf[Pstream::myProcNo()].setSize(tloc.size());
    rf[Pstream::myProcNo()].setSize(refField_.size());

    label i = 0;
    forAll(tloc, I)
    {
        l[Pstream::myProcNo()][i]  = LField_[I];
        cf[Pstream::myProcNo()][i] = tloc[I];
        rf[Pstream::myProcNo()][i] = refField_[I];

        ++i;
    }

    Pstream::gatherList(l);
    Pstream::gatherList(cf);
    Pstream::gatherList(rf);

    List<scalar> L  = ListListOps::combine<List<scalar> >(l, accessOp<List<scalar> >());
    List<vector> CF = ListListOps::combine<List<vector> >(cf, accessOp<List<vector> >());
    List<vector> RF = ListListOps::combine<List<vector> >(rf, accessOp<List<vector> >());

    List<List<decayingVorton> > v(Pstream::nProcs());

    if (Pstream::master())
    {
        Pout<<"Starting inflgen, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
        forAll(CF, I)
        {
            if (L[I] > 0.0)
            {
                scalar x    = CF[I].x();
                scalar ls   = lspot(L[I]);
                scalar ymin = CF[I].y() - L[I];
                scalar zmin = CF[I].z() - L[I];

                for (label k = 0; k < NUM_VORT; k++)
                {
                    vector v((direction_ > 0) ? x-ls : x+ls, ymin+ranGen.scalar01()*2*L[I], zmin+ranGen.scalar01()*2*L[I]);
                    bool allowed = true;
                    for (SLList<decayingVorton>::const_iterator it = vortons_.begin(); it != vortons_.end(); ++it)
                    {
                        if (mag(v - it().location()) < OVERLAP*ls)
                        {
                            allowed = false;
                            break;
                        }
                    }
                    if (!allowed)
                        continue;
                    else
                    {
                        vortons_.insert(decayingVorton(L[I], v, RF[I], (direction_ > 0) ? x+ls : x-ls));
                        if(debug)
                        {
                            nVortField_[I]++;
                        }
                    }
                }
            }
        }

        Pout<<"The number of vortons is "<<vortons_.size()<<endl;

        v[Pstream::myProcNo()].setSize(vortons_.size());
        i = 0;
        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it != vortons_.end(); ++it)
            v[Pstream::myProcNo()][i++] = it();
    }

    Pstream::scatterList(v);

    List<decayingVorton> V = ListListOps::combine<List<decayingVorton> >(v, accessOp<List<decayingVorton> >());

    if (size() > 0)
    {
        Field<vector> turbulent(refField_.size(), pTraits<vector>::zero);

        forAll(patchField, I)
        {
            vector pt = tloc[I];

            forAll(V, vI)
            {
                const decayingVorton& vorton = V[vI];

                if (mag(vorton.location() - pt) < lspot(vorton.length()))
                    turbulent[I] += vorton.velocityAt(pt);
            }
        }

        R_ = ((ind_ - 1.0)/ind_)*R_ + (1.0/ind_)*sqr(turbulent);

        Field<tensor> C_(R_.size(), pTraits<tensor>::zero);
        forAll(C_, I)
        {
            C_[I].xx() = 1.0;
            C_[I].yy() = 1.0;
            C_[I].zz() = 1.0;
        }

        forAll(R_, I)
        {
            scalar D1 = R_[I].xx();
            if (D1 > 0)
                C_[I].xx() = 1.0/sqrt(D1);

            scalar D2 = R_[I].xx()*R_[I].yy() - R_[I].xy()*R_[I].xy();
            if (D1 > 0 && D2 > 0)
            {
                C_[I].yx() = -R_[I].xy()/sqrt(D1*D2);
                C_[I].yy() = sqrt(D1/D2);
            }

            scalar D3 = det(R_[I]);
            if (D2 > 0 && D3 > 0)
            {
                C_[I].zx() = (R_[I].xy()*R_[I].yz()-R_[I].yy()*R_[I].xz())/sqrt(D2*D3);
                C_[I].zy() = -(R_[I].xx()*R_[I].yz()-R_[I].xz()*R_[I].xy())/sqrt(D2*D3);
                C_[I].zz() = sqrt(D2/D3);
            }
        }

        ind_++;

        turbulent = C_&turbulent;
        turbulent = Lund_&turbulent;

        fixedValueFvPatchField<vector>::operator==(refField_+ turbulent);
    }

    if (Pstream::master())
    {
        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it != vortons_.end(); ++it)
        it().move(this->db().time().deltaT().value());

        bool modified;
        do
        {
            modified=false;
            for (SLList<decayingVorton>::iterator it = vortons_.begin(); it!=vortons_.end(); ++it)
            {
                if (direction_ > 0)
                    {
                        if (it().location().x() > it().xmax())
                        {
                            vortons_.remove(it);
                            modified=true;
                            break;
                        }
                    }
                    else
                    {
                        if (it().location().x() < it().xmax())
                        {
                            vortons_.remove(it);
                            modified=true;
                            break;
                        }
                    }
            }
        } while (modified);

        Pout<<"Finishing inflgen, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
    }
}

void Foam::decayingTurbulenceFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "LField", LField_);
    writeEntry(os, "refField", refField_);
    writeEntry(os, "RField", RField_);
    os.writeKeyword("direction")<<direction_<<token::END_STATEMENT<<nl;
    os.writeKeyword("inletShape")<<inletShape_<<token::END_STATEMENT<<nl;
    if (inletShape_==1)
    {
        os.writeKeyword("Umean")<<Umean_<<token::END_STATEMENT<<nl;
        os.writeKeyword("width")<<width_<<token::END_STATEMENT<<nl;
        os.writeKeyword("midRadius")<<midRadius_<<token::END_STATEMENT<<nl;
        os.writeKeyword("center")<<center_<<token::END_STATEMENT<<nl;
    }
    else if (inletShape_==2)
    {
        os.writeKeyword("Umean")<<Umean_<<token::END_STATEMENT<<nl;
        os.writeKeyword("center")<<center_<<token::END_STATEMENT<<nl;
        os.writeKeyword("Radius")<<Radius_<<token::END_STATEMENT<<nl;
    }
    else
    {
        FatalErrorInFunction
            << "You can shoose 'ring', 'tube', instead of inletShape_!" << endl;
    }
    writeEntry(os, "value", *this);
    writeEntry(os, "R", R_);
    os.writeKeyword("ind")<<ind_<<token::END_STATEMENT<<nl;

    /*if (Pstream::master())
    {
        fileName outputFile;
        const fileOperation& fp = Foam::fileHandler();
        if (isA<Foam::fileOperations::collatedFileOperation>(fp))
        {
            outputFile=db().time().rootPath()/db().time().globalCaseName()/"processors"+Foam::name(Pstream::nProcs())/db().time().timeName()/+"vortons_"+patch().name();
        }
        else if (isA<Foam::fileOperations::uncollatedFileOperation>(fp))
        {
            outputFile = db().time().path()/db().time().timeName()/+"vortons_"+patch().name();
        }
            OFstream os2(outputFile);
            os2.writeKeyword("vortons")<<vortons_<<token::END_STATEMENT<<nl;
    }*/

    // 3. write in the boundary part of U
    if (Pstream::master())
       os.writeKeyword("vortons")<<vortons_<<token::END_STATEMENT<<nl;
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        decayingTurbulenceFvPatchVectorField
    );
}

