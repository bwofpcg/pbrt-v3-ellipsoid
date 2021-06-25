//
//  ellipsoid.hpp
//  PBRT-V3
//
//  Created by bjw on 1/12/21.
//
// This file implements Ellipsoid BRDF lobes when encoded using an S matrix form (a symmetric 3x3 matrix)
// Note: this version does not include a fresnel term but does inlcude an RGB scale factor 
// It also implements a material that can include a lambertian lobe and up to two ellipsoid lobes

#ifndef ellipsoid_h
#define ellipsoid_h

#include "pbrt.h"
#include "material.h"
#include "reflection.h"

namespace pbrt {

// A single ellipsoid lobe BRDF parameterized by a symmetric matrix S
class EllipsoidSReflection : public BxDF {
  public:
    EllipsoidSReflection(const Spectrum &scaling,
                         const Vector3f &Sdiagonals, const Vector3f &Soffdiagonals)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          scaling(scaling),
          Sdiagonals(Sdiagonals),
          Soffdiagonals(Soffdiagonals) { computeSInverse(); }
    // isotropic GGX constructor
    EllipsoidSReflection(const Spectrum &scaling, Float alpha) : EllipsoidSReflection(scaling,Vector3f(1.0/(alpha*alpha),1.0/(alpha*alpha),1.0)*alpha*sqrt(Pi),Vector3f(0.0,0.0,0.0)) {}
    // isotropic GGX constructor
    EllipsoidSReflection(const Spectrum &scaling, Float alphaX, Float alphaY) : EllipsoidSReflection(scaling,Vector3f(1.0/(alphaX*alphaX),1.0/(alphaY*alphaY),1.0)*sqrt(Pi*alphaX*alphaY),Vector3f(0.0,0.0,0.0)) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // This implementation of an ellipsoid microfacet BRDF is parameterized by a symmetric matrix S
    // In terms of the A matrix from the original ellipsoid definition, S is defined as:
    // S = C * (A^T * A)^(-1)  where ^(-1) is the matrix inverse and C can be any arbitrary constant.
    // The corresponding NDF is D(h) = 1 / (PI * sqrt(det(S^-1)*abs(n^T*(S^-1)*n)) * (h^T * S * h)^2
    //
    // Note: by convention we often normalize the S matrix by selecting C such such
    //   that 1 == PI * sqrt(det(S^-1)*abs(n^T*(S^-1)*n) but this is not required and C can simply be set to one.
    const Spectrum scaling;  //scale factor for this BRDF
    //note: S is a symmetric 3x3 matrix which is specified by 6 elements (3 diagonals and 3 off-diagonals)
    //the diagonals are (Sxx,Syy,Szz) and the off-diagonals are (Sxy,Syz,Szx)
    const Vector3f Sdiagonals;  //diagonal elements of S in RGB channels
    const Vector3f Soffdiagonals;  //off-diagonals elements of S in RGB channels
    //we also need the inverse of the symmetric matrix S
    Vector3f invSdiagonals;  //diagonal elements of inverse_S in RGB channels
    Vector3f invSoffdiagonals;  //off-diagonals elements of inverse S in RGB channels
    Float detInvS;  //the determinant of the matrix S_inverse
    
    // compute the inverse of S from the value of S
    void computeSInverse();
    // compute the shadowing/masking factor (note: uses value of inverseS)
    Float computeShadMask(const Vector3f &wo, const Vector3f &wi) const;
};


// Need to wrap a material around ellipsoid lobe for testing
// ellipsoid + optional lambertian term as a simple way to test ellipsoid BRDF
class EllipsoidSMaterial : public Material {
  public:
    EllipsoidSMaterial(const std::shared_ptr<Texture<Spectrum>> &lambertianAlbedo,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid0_scale,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid0_S_diagonals,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid0_S_offdiagonals,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid1_scale,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid1_S_diagonals,
                    const std::shared_ptr<Texture<Spectrum>> &ellipsoid1_S_offdiagonals
                    )
        : lambAlbedo(lambertianAlbedo),
          e0scale(ellipsoid0_scale),
          e0Sdiag(ellipsoid0_S_diagonals), e0Soffdiag(ellipsoid0_S_offdiagonals),
          e1scale(ellipsoid1_scale),
          e1Sdiag(ellipsoid1_S_diagonals), e1Soffdiag(ellipsoid1_S_offdiagonals) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> lambAlbedo;
    //note these textures are really 3d vectors, not actually colors or spectra but pbrt only supports float or spectrum textures
    std::shared_ptr<Texture<Spectrum>> e0scale, e0Sdiag, e0Soffdiag;  
    std::shared_ptr<Texture<Spectrum>> e1scale, e1Sdiag, e1Soffdiag;  
};

EllipsoidSMaterial *CreateEllipsoidSMaterial(const TextureParams &mp);

} //namespace pbrt

#endif /* ellipsoid_h */
