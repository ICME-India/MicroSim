#ifndef readPhysicalDiftolmxe2ehfeolkfjj_H
#define readPhysicalDiftolmxe2ehfeolkfjj_H

#include "scalar.H"
#include "fvScalarMatrix.H"
#include "SquareMatrix.H"
#include "Time.H"

namespace Foam::mpmsc
{
    #define str(arg) std::to_string(arg)
    
    using FieldPtr  = autoPtr<volScalarField>;
    using Field3Ptr = autoPtr<volVectorField>;
    using FieldArr  = List<FieldPtr>;

    using tmpReal   = tmp<volScalarField>;
    using tmpVec3   = tmp<volVectorField>;

    const vector dir_x(1, 0, 0);
    const vector dir_y(0, 1, 0);
    const vector dir_z(0, 0, 1);

    const dimensionedScalar 
        dt      ("dt",  dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.0),
        idx     ("dx",  dimensionSet(0, 1, 0, 0, 0, 0, 0), 1.0),
        idx2    = idx * idx;


    inline auto fun_g_prime(const volScalarField& phi)
    {
        return 2.0 * phi * (1.0 - phi) * (1.0 - 2.0*phi);
    };


    inline auto fun_H_phi(const volScalarField& phi)
    {
        return 6.0 * phi * (1.0 - phi);
    };


    struct GrainFieldConstants_t
    {
        // autoPtr<IOdictionary> IODict;

        label   NumPhases;
        scalar  Epsilon,
                Mobility_sol,
                Mobility_liq,
                Gamma,
                Coeff_A1,
                Coeff_A2,
                Coeff_A3,
                Coeff_Asym,
                Coeff_H_phi,
                Coeff_AntiCap,
                Tequlibrium,
                Treference;

        vector  TGradient,
                pullingVelocity,
                TreferencePoint; // point at which reference temperature is given

        List<scalar>
                Angles_X,
                Angles_Y,
                Angles_Z;
                
        bool    UseTemperatureCorrector;
        
        // calculated conts
        List<tensor> Rotation;


        // functions
        template<class T>
        inline void read(const IOdictionary& dict, const word& key, T& value)
        {
            value = dict.lookup<T>(key);
        }

        GrainFieldConstants_t() = default;
        GrainFieldConstants_t(
            const IOdictionary& IODict, const Foam::fvMesh& mesh, const Foam::Time &runTime);

        inline void printConstants() const;
    };


    inline tensor getrotationmatrix_z(scalar angle) noexcept
    {
        return tensor
        (
            cos(angle), -sin(angle), 0,
            sin(angle),  cos(angle), 0,
            0, 0, 1
        );   
    }


    inline auto createVSF(const Foam::fvMesh &mesh, word name)
    {
        return new volScalarField(
            IOobject(
                name, 
                mesh.time().name(),
                mesh,
                IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
            ), mesh,
            dimensionedScalar("zero", dimless, 0.0) 
        );
    }


    inline auto createVSF_NO(const Foam::fvMesh &mesh, word name)
    {
        return new volScalarField(
            IOobject(
                name, mesh,
                IOobject::READ_IF_PRESENT, IOobject::NO_WRITE               //NO_READ
            ), mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0)
        );
    }


    inline auto createVVF(const Foam::fvMesh &mesh, word name)
    {
        return new volVectorField(
            IOobject(
                name, mesh,
                IOobject::NO_READ, IOobject::NO_WRITE
            ), mesh,
            Foam::dimensionedVector("zero", Foam::dimless, point(0,0,0))
        );
    }

    
    inline GrainFieldConstants_t::GrainFieldConstants_t(
        const IOdictionary& iod,
        const Foam::fvMesh& mesh, const Foam::Time &runTime)
    {

        read(iod, "NumPhases", NumPhases);

        read(iod, "EPSILON",          Epsilon);
        read(iod, "GG_Mobility_sol",  Mobility_sol);
        read(iod, "GG_Mobility_liq",  Mobility_liq);
        read(iod, "GG_Gamma",         Gamma);
        read(iod, "GG_Coeff_A1",      Coeff_A1);
        read(iod, "GG_Coeff_A2",      Coeff_A2);
        read(iod, "GG_Coeff_A3",      Coeff_A3);

        // read(iod, "GG_Coeff_Asym",    Coeff_Asym,   0.0);
        read(iod, "GG_Coeff_H_phi",   Coeff_H_phi);
        read(iod, "GG_Coeff_AntiCap",   Coeff_AntiCap);

        // read(iod, "GG_Angles_Z",      Angles_Z);
        // read(iod, "GG_Angles_Y",      Angles_Y);
        // read(iod, "GG_Angles_X",      Angles_X);

        read(iod, "Temperature_Equilibrium",     Tequlibrium);
        read(iod, "Temperature_Base",            Treference);
        read(iod, "Temperature_Gradient",        TGradient);
        read(iod, "PullingVelocity",             pullingVelocity);
        //read(iod, "Temperature_BasePoint",       TreferencePoint);

        // TreferencePoint += vector(0.5, 0.5, 0.5);
        
        read(iod, "UseTemperatureCorrector",  UseTemperatureCorrector);   

        auto anglesz = List<scalar>(iod.lookup("GG_Angles_Z"));

        Angles_Z.setSize(NumPhases);
        Rotation.setSize(NumPhases); 
        forAll(Angles_Z, i)
        {
            Angles_Z[i] = anglesz[i];
            Rotation[i] = getrotationmatrix_z(degToRad(-1.0*Angles_Z[i]));
        }
    }


    inline void GrainFieldConstants_t::printConstants() const
    {
        Info << nl
            << "========== Physical Constants ==========" << nl

            << "NumPhases                 : " << NumPhases << nl

            << "EPSILON                   : " << Epsilon << nl
            << "GG_Mobility_sol           : " << Mobility_sol << nl
            << "GG_Mobility_liq           : " << Mobility_liq << nl
            << "GG_Gamma                  : " << Gamma << nl
            << "GG_Coeff_A1               : " << Coeff_A1 << nl
            << "GG_Coeff_A2               : " << Coeff_A2 << nl
            << "GG_Coeff_A3               : " << Coeff_A3 << nl
            << "Coeff_H_phi               : " << Coeff_H_phi << nl
            << "Coeff_AntiCap               : " << Coeff_AntiCap << nl

            << "GG_Angles_Z               : " << Angles_Z << nl
            << "GG_Angles_Y               : " << Angles_Y << nl
            << "GG_Angles_X               : " << Angles_X << nl

            << "Temperature_Equilibrium   : " << Tequlibrium << nl
            << "Temperature_Base          : " << Treference << nl
            << "Temperature_Gradient      : " << TGradient << nl
            << "PullingVelocity           : " << pullingVelocity << nl
            //<< "Temperature_BasePoint     : " << TreferencePoint << nl

            << "UseTemperatureCorrector   : " << UseTemperatureCorrector << nl

            << "Rotation                  : " << Rotation << nl

            << "========================================" << endl;
    }


    scalar getdx(const volScalarField& field)
    {
        const auto& mesh = field.mesh();
        auto V = mesh.V()[0];
        return cbrt(V);
    }


    void fillLastPhase(FieldArr& sources)
    {
        auto&& lastPhase = (*sources.last());
        forAll(lastPhase, idx)
        {
            scalar sum = 0;
            for (label i = 0; i < sources.size()-1; i++)
            {
                sum += (*sources[i])[idx];
            }
            lastPhase[idx] = 1.0-sum;
        }
    }


    void fillGradient(volScalarField& field, 
        scalar start, scalar delta, int axis)
    {
        const auto& mesh = field.mesh();
        field = start + mesh.C().component(axis) * delta / idx;
    }


    void fillTemperatureGradient(volScalarField& T, 
        scalar Treference, vector TGradient, vector TreferencePoint)
    {
        const fvMesh& mesh = T.mesh();
        const vectorField& C = mesh.C();        
        scalarField& Tint = T.primitiveFieldRef();
        forAll(Tint, celli)
        {
            const vector dX = C[celli] - TreferencePoint;
            Tint[celli] = Treference + (dX & TGradient);
        }
        T.correctBoundaryConditions();
    }


    void fillSphere(volScalarField& field, 
        const point& centre, const scalar radius, const scalar value)
    {
        const auto& mesh = field.mesh();
        forAll(field, idx)
        {
            const point& c = mesh.C()[idx];
            if (mag(c - centre) <= radius)
            {
                field[idx] = value;
            }
        }
    }
     
        
    void fillCuboid
    (
    volScalarField& field,
    const point& minCorner, 
    const point& maxCorner, 
    const scalar value
    )
    {
        const auto& mesh = field.mesh();

        forAll(field, idx)
        {
            const point& c = mesh.C()[idx];
            if (c.x() >= minCorner.x() && c.x() <= maxCorner.x() &&
                c.y() >= minCorner.y() && c.y() <= maxCorner.y() &&
                c.z() >= minCorner.z() && c.z() <= maxCorner.z())
            {
                field[idx] = value;
            }
        }   
    }
    
    
}

#endif
