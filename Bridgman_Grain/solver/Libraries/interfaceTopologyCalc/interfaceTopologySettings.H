/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

Class
    Foam::interfaceTopologyUtils::interfaceTopologySettings

Description
    Caches values for one subdict entry.
    Merges entry-specific values with defaults:
        1. Start with defaults
        2. Override with entry-specific values
    Entry-specific values take priority over defaults.

\*---------------------------------------------------------------------------*/

#ifndef interfaceTopologySettings_H
#define interfaceTopologySettings_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "dimensionSets.H"
#include "wordList.H"
#include "scalarList.H"
#include "word.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTopologyUtils
{

/*---------------------------------------------------------------------------*\
    interfaceTopologySettings — parameters for one field segmentation
\*---------------------------------------------------------------------------*/

class interfaceTopologySettings
{
    word              subdictName_;
    word              fieldName_;
    word              BPsiName_;
    scalar            Tm_;
    scalar            limT_;
    dimensionedScalar idimx_;
    scalar            ybottom_;
    scalar            yRef4_;
    scalar            interiorThreshold_;
    scalar            boundaryThreshold_;
    scalar            interfaceMinCells_;
    scalar            maskTolerance_;
    scalar            regionMaskEps_;

    scalar readScalar
    (
        const dictionary& entry,
        const dictionary& defaults,
        const word& key,
        const scalar fallback
    ) const
    {
        if (entry.found(key))   return entry.lookup<scalar>(key);
        if (defaults.found(key)) return defaults.lookup<scalar>(key);
        return fallback;
    }

    word readWord
    (
        const dictionary& entry,
        const dictionary& defaults,
        const word& key,
        const word& fallback
    ) const
    {
        if (entry.found(key))   return entry.lookup<word>(key);
        if (defaults.found(key)) return defaults.lookup<word>(key);
        return fallback;
    }

    void readFromDicts
    (
        const dictionary& entry,
        const dictionary& defaults
    )
    {
        // fieldName — entry key 'fieldName' first, then subdict name
        // supports both old style (subdict name = field name)
        // and new style (explicit fieldName or field keyword)
        if (entry.found("fieldName"))
            fieldName_ = entry.lookup<word>("fieldName");
        if (entry.found("FieldName"))
            fieldName_ = entry.lookup<word>("FieldName");    
        else if (entry.found("field"))
            fieldName_ = entry.lookup<word>("field");
        else if (entry.found("Field"))
            fieldName_ = entry.lookup<word>("Field");    
        else
            fieldName_ = subdictName_;  // fallback — subdict name is field name

        // BPsiName — three-level priority
        if (entry.found("BPsiName"))
        {
            BPsiName_ = entry.lookup<word>("BPsiName");
        }
        else if
        (
            defaults.isDict("BPsiNames")
         && defaults.subDict("BPsiNames").found(subdictName_)
        )
        {
            BPsiName_ =
                defaults.subDict("BPsiNames").lookup<word>(subdictName_);
        }
        else if (defaults.found("BPsiName"))
        {
            BPsiName_ = defaults.lookup<word>("BPsiName");
        }
        else
        {
            FatalErrorInFunction
                << "BPsiName not resolved for entry '" << subdictName_ << "'"
                << nl << "Set BPsiName in entry, default.BPsiNames, or default"
                << abort(FatalError);
        }

        Tm_      = readScalar(entry, defaults, "Tm",      0);
        limT_    = readScalar(entry, defaults, "limT",    0);
        ybottom_ = readScalar(entry, defaults, "ybottom", 0);
        yRef4_   = readScalar(entry, defaults, "yRef4",   0);

        idimx_ = dimensionedScalar
        (
            "idimx",
            dimless/dimLength,
            readScalar(entry, defaults, "idimx", 1)
        );

        interiorThreshold_ =
            readScalar(entry, defaults, "interiorThreshold", 0.96);
        boundaryThreshold_ =
            readScalar(entry, defaults, "boundaryThreshold", 0.05);
        interfaceMinCells_ =
            readScalar(entry, defaults, "interfaceMinCells", 5);
        maskTolerance_ =
            readScalar(entry, defaults, "maskTolerance", 2.0);
        regionMaskEps_ =
            readScalar(entry, defaults, "regionMaskEps", 0.001);
    }

public:

    interfaceTopologySettings
    (
        const word& subdictName,
        const dictionary& entry,
        const dictionary& defaults
    )
    :
        subdictName_(subdictName),
        fieldName_(subdictName),             //overwritten by readFromDicts
        BPsiName_("BPsi"),
        Tm_(0), limT_(0),
        idimx_("idimx", dimless/dimLength, 1),
        ybottom_(0), yRef4_(0),
        interiorThreshold_(0.96),
        boundaryThreshold_(0.05),
        interfaceMinCells_(5),
        maskTolerance_(2.0),
        regionMaskEps_(0.001)
    {
        readFromDicts(entry, defaults);
    }

    interfaceTopologySettings
    (
        const word& subdictName,
        const dictionary& entry
    )
    :
        interfaceTopologySettings(subdictName, entry, dictionary::null)
    {}

    void update(const dictionary& entry, const dictionary& defaults)
    {
        readFromDicts(entry, defaults);
        Info<< "interfaceTopologySettings: updated '" << subdictName_ << "'" << nl;
        print();
    }

    void update(const dictionary& entry)
    {
        update(entry, dictionary::null);
    }

    const word& subdictName()        const { return subdictName_;        }
    const word& fieldName()          const { return fieldName_;          }
    const word& BPsiName()           const { return BPsiName_;           }
    scalar      Tm()                 const { return Tm_;                 }
    scalar      limT()               const { return limT_;               }
    const dimensionedScalar& idimx() const { return idimx_;              }
    scalar      ybottom()            const { return ybottom_;            }
    scalar      yRef4()              const { return yRef4_;              }
    scalar      interiorThreshold()  const { return interiorThreshold_;  }
    scalar      boundaryThreshold()  const { return boundaryThreshold_;  }
    scalar      interfaceMinCells()  const { return interfaceMinCells_;  }
    scalar      maskTolerance()      const { return maskTolerance_;      }
    scalar      regionMaskEps()      const { return regionMaskEps_;      }

    void print() const
    {
        Info<< "\ninterfaceTopologySettings [" << subdictName_ << "]:" << nl
            << "    fieldName         : " << fieldName_         << nl
            << "    BPsiName          : " << BPsiName_          << nl
            << "    Tm                : " << Tm_                << nl
            << "    limT              : " << limT_              << nl
            << "    idimx             : " << idimx_.value()     << nl
            << "    ybottom           : " << ybottom_           << nl
            << "    yRef4             : " << yRef4_             << nl
            << "    interiorThreshold : " << interiorThreshold_ << nl
            << "    boundaryThreshold : " << boundaryThreshold_ << nl
            << "    interfaceMinCells : " << interfaceMinCells_ << nl
            << "    maskTolerance     : " << maskTolerance_     << nl
            << "    regionMaskEps     : " << regionMaskEps_     << endl;
    }
};


/*---------------------------------------------------------------------------*\
    interfaceComboSettings — parameters for a combined field entry
\*---------------------------------------------------------------------------*/

class interfaceComboSettings
{
public:

    enum class operation { multiply, add, min, max };

    enum class entryMode
    {
        raw,       // use field values directly
        binarized  // binarize each field before combining
    };

private:

    word            comboName_;
    wordList        factorEntries_;  // entry names e.g. "T", "CalcT[0]"
    operation       operation_;
    entryMode       entryMode_;
    scalarList      weights_;
    word            BPsiName_;

    static operation readOperation(const word& w)
    {
        if (w == "multiply") return operation::multiply;
        if (w == "add")      return operation::add;
        if (w == "min")      return operation::min;
        if (w == "max")      return operation::max;

        FatalErrorInFunction
            << "Unknown combo operation '" << w << "'" << nl
            << "Valid: multiply add min max"
            << abort(FatalError);
        return operation::multiply;
    }

public:

    interfaceComboSettings
    (
        const word& comboName,
        const dictionary& entry,
        const dictionary& defaults
    )
    :
        comboName_(comboName),
        factorEntries_(entry.lookup("factorEntries")),
        operation_(readOperation
        (
            entry.lookupOrDefault<word>("operation", "multiply")
        )),
        entryMode_
        (
            entry.lookupOrDefault<word>("entryMode", "binarized") == "raw"
          ? entryMode::raw
          : entryMode::binarized
        ),
        BPsiName_
        (
            entry.found("BPsiName")
          ? entry.lookup<word>("BPsiName")
          : defaults.lookupOrDefault<word>("BPsiName", "BPsi")
        )
    {
        if (factorEntries_.size() < 2)
        {
            FatalErrorInFunction
                << "combo '" << comboName_
                << "' must have at least 2 factorEntries"
                << abort(FatalError);
        }

        if (entry.found("weights"))
        {
            weights_ = scalarList(entry.lookup("weights"));

            if (weights_.size() != factorEntries_.size())
            {
                FatalErrorInFunction
                    << "weights size (" << weights_.size()
                    << ") != factorEntries size ("
                    << factorEntries_.size() << ")"
                    << " in combo '" << comboName_ << "'"
                    << abort(FatalError);
            }
        }
        else
        {
            weights_.setSize(factorEntries_.size());
            weights_ = scalar(1) / factorEntries_.size();
        }
    }

    const word&       comboName()     const { return comboName_;     }
    const wordList&   factorEntries() const { return factorEntries_; }
    operation         op()            const { return operation_;     }
    entryMode         mode()          const { return entryMode_;     }
    const scalarList& weights()       const { return weights_;       }
    const word&       BPsiName()      const { return BPsiName_;      }

    void print() const
    {
        Info<< "\ninterfaceComboSettings [" << comboName_ << "]:" << nl
            << "    factorEntries : " << factorEntries_ << nl
            << "    BPsiName      : " << BPsiName_      << endl;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTopologyUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
