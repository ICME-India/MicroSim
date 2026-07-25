//! The interpolation for smoothening of the fields
surfaceScalarField("alpha", fvc::interpolate(alpha)),
surfaceScalarField("mu", fvc::interpolate(mu)),
geometricZeroField(),
simple
