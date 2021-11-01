//
//  ellipsoid.cpp
//  PBRT-V3
//
//  Created by bjw on 1/12/21.
//

#include "ellipsoid.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {
//some routines for working with symmetric 3x3 matrices encoded as two vector: 
//    one for the diagonal elements and one for the offdiagonals 

// compute the two sided product of a vector with a symmetric matrix (ie v*S*v)
static inline Float compute_vSv(const Vector3f &Sdiag, const Vector3f &Soffdiag, const Vector3f &v) {
    return Sdiag.x*v.x*v.x + Sdiag.y*v.y*v.y + Sdiag.z*v.z*v.z
    + 2*(Soffdiag.x*v.x*v.y + Soffdiag.y*v.y*v.z + Soffdiag.z*v.z*v.x);
}

//compute z*S*v where S is a symmetric matrix and z is the unit vector 0,0,1
static inline Float compute_zdirSv(const Vector3f &Sdiag, const Vector3f &Soffdiag, const Vector3f &v) {
    return v.x*Soffdiag.z + v.y*Soffdiag.y + v.z*Sdiag.z;
}

static inline Vector3f compute_Sv(const Vector3f &Sdiag, const Vector3f &Soffdiag, const Vector3f &v) {
    return Vector3f(v.x*Sdiag.x + v.y*Soffdiag.x + v.z*Soffdiag.z,
                    v.x*Soffdiag.x + v.y*Sdiag.y + v.z*Soffdiag.y,
                    v.x*Soffdiag.z + v.y*Soffdiag.y + v.z*Sdiag.z);
}

static inline Float compute_determinantS(const Vector3f &Sdiag, const Vector3f &Soffdiag) {
    return Sdiag.x*Sdiag.y*Sdiag.z + 2*(Soffdiag.x*Soffdiag.y+Soffdiag.z)
    - Sdiag.x*Soffdiag.y*Soffdiag.y - Sdiag.y*Soffdiag.z*Soffdiag.z - Sdiag.z*Soffdiag.x*Soffdiag.x;
}

// We also need the inverse matrix for S which can be computed here
void EllipsoidSReflection::computeSInverse() {
    //compute inverse of symmetric matrix S using cofactors.
    Float Vxx = Sdiagonals.z*Sdiagonals.y - Soffdiagonals.y*Soffdiagonals.y;
    Float Vxy = Soffdiagonals.y*Soffdiagonals.z - Sdiagonals.z*Soffdiagonals.x;
    Float Vzx = Soffdiagonals.y*Soffdiagonals.x - Sdiagonals.y*Soffdiagonals.z;
    Float invDet = 1.0 / (Sdiagonals.x*Vxx + Soffdiagonals.x*Vxy + Soffdiagonals.z*Vzx);
    Float Vyy = Sdiagonals.z*Sdiagonals.x - Soffdiagonals.z*Soffdiagonals.z;
    Float Vzz = Sdiagonals.y*Sdiagonals.x - Soffdiagonals.x*Soffdiagonals.x;
    Float Vyz = Soffdiagonals.x*Soffdiagonals.z - Soffdiagonals.y*Sdiagonals.x;
    invSdiagonals = Vector3f(invDet*Vxx, invDet*Vyy, invDet*Vzz);
    invSoffdiagonals = Vector3f(invDet*Vxy, invDet*Vyz, invDet*Vzx);
    detInvS = invDet;
//    printf("invS diag %f %f %f off %f %f %f\n",invSdiagonals.x,invSdiagonals.y,invSdiagonals.z,invSoffdiagonals.x,invSoffdiagonals.y,invSoffdiagonals.z);
}

// compute the shadowing/masking factor (note: uses value of inverseS)
Float EllipsoidSReflection::computeShadMask(const Vector3f &wo, const Vector3f &wi) const {
    const bool useBlackVertical = true;   //use shadow/masking variant with extra vertical black "microfacets" so that mean normal matches macro-surface normal
    const bool useHeightCorrelatedShadMask = false;  //use height correlated version of bidirectional shadow/masking
    //note that the surface normal in local coordinate is +Z (ie n = 0,0,1)
    //and that invS is proportional to (A^T)*A and that
    //since the fraction is invariant to scalings of A, we can use invS in place of (A^T)*T
    Float absCosView = AbsCosTheta(wo);
    Float absCosLight = AbsCosTheta(wi);
    Float nAAn = invSdiagonals.z;
    Float vAAv = compute_vSv(invSdiagonals, invSoffdiagonals, wo);
    Float lAAl = compute_vSv(invSdiagonals, invSoffdiagonals, wi);
    Float vAAn = compute_zdirSv(invSdiagonals, invSoffdiagonals, wo);
    Float lAAn = compute_zdirSv(invSdiagonals, invSoffdiagonals, wi);
    Float numeratorV = 2*nAAn*absCosView;
    Float numeratorL = 2*nAAn*absCosLight;
    if (CosTheta(wo) < 0) { vAAn = -vAAn; }  //Not totally sure handling of directions from backside is correct, needs testing
    if (CosTheta(wi) < 0) { lAAn = -lAAn; }
    Float Gv, Gl;
    if (useBlackVertical) {
        Gv = numeratorV / ( sqrt(vAAv*nAAn) + std::max(vAAn, numeratorV-vAAn) );
        Gl = numeratorL / ( sqrt(lAAl*nAAn) + std::max(lAAn, numeratorL-lAAn) );
    } else {
        Gv = numeratorV / (sqrt(vAAv*nAAn) + vAAn);
        Gl = numeratorL / (sqrt(lAAl*nAAn) + lAAn);
        if (Gv > 1.0) Gv = 1.0;
        if (Gl > 1.0) Gl = 1.0;
    }
    if (useHeightCorrelatedShadMask) {
        return Gv*Gl / (Gv + Gl - Gv*Gl);
    }
    return Gv*Gl;
}

//Note matrix S may follow one of three conventions:
//unnormalized - S may contain an arbirary scalefactor and a normalization factor is needed when evaluating
//normalized - S has been scaled so that its normalization factor equal one, and thus has no effect
//scaled-normalized - S has been normalized and then had an external scalefactor applied (eg to incorporate an albedo).  In this case the normalization factor should not be applied during evaluation
//Generally we normalize the S matrices beforehand, but they may become unormalized due to texture interpolation or other factors, so in this code we apply the normalization factor

Spectrum EllipsoidSReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    if (wi.z <= 0) return Spectrum(0.f);
    if (wo.z <= 0) return Spectrum(0.f);  //wrong side of surface
    Float absCosView = AbsCosTheta(wo);
    Float absCosLight = AbsCosTheta(wi);
    Vector3f h = wi + wo;   //unnormalized half vector
    Float hLen2 = h.LengthSquared();
    Float hSh = compute_vSv(Sdiagonals, Soffdiagonals, h);
    Float normalization = Pi * sqrt(detInvS * invSdiagonals.z);
    Float val = hLen2*hLen2 / (hSh*hSh * 4*absCosView*absCosLight * normalization);
    val *= computeShadMask(wo,wi);
    return scaling*val;
}

static inline std::string VecToString(const Vector3f &v) {
    return StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
}

std::string EllipsoidSReflection::ToString() const {
    return std::string("[ EllipsoidSReflection scaling: ") + scaling.ToString() +
           std::string(" Sdiagonal: ") + VecToString(Sdiagonals) +
           std::string(" Soffdiagonals: ") + VecToString(Soffdiagonals) + std::string(" ]");
}

//Implements Shirley et al.'s concentric mapping from a square to a circle
//This better preserves structure and stratification as compared to a simple radius,angle mapping
static Point2f ConcentricMapping_squareToCircle(const Point2f &inRandom) {
  Float x,y;
  Float dx = 2*inRandom.x - 1.0;
  Float dy = 2*inRandom.y - 1.0;
  Float adx = abs(dx);
  Float ady = abs(dy);
  if ((adx < 1e-8)&&(ady < 1e-8)) { //special case to avoid divide by zero problems
    x = y = 0.0;
  } else if (adx > ady) {
    Float angle = (Pi/4.0)*dy/dx;
    double r = dx;
    x = r*cos(angle);
    y = r*sin(angle);
  } else {
    double angle = (Pi/4.0)*dx/dy;
    double r = dy;
    x = r*sin(angle);
    y = r*cos(angle);
  }
  return Point2f(x,y);
}

Spectrum EllipsoidSReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                        const Point2f &u, Float *pdf,
                                        BxDFType *sampledType) const {
    // Sample microfacet orientation $\wh$ and reflected direction $\wi$
    if (wo.z == 0) return 0.;
    //(a) sample random point on the unit disk
    Point2f ptA;
    bool useConcentricMapping = true;
    if (useConcentricMapping) {
        ptA =  ConcentricMapping_squareToCircle(u);
    } else {
        Float r = sqrt(u.x);
        Float theta = 2*Pi*u.y;
        ptA = Point2f(r*cos(theta),r*sin(theta));
    }
    //(b) compress to crescent and project to sphere
    Float vAAv = compute_vSv(invSdiagonals, invSoffdiagonals, wo);
    Float vAAn = compute_zdirSv(invSdiagonals, invSoffdiagonals, wo);
    Float nAAn = invSdiagonals.z;
    Float s = 0.5*(1 + vAAn/sqrt(vAAv * nAAn));
    Point3f ptB;
    ptB.x = s*ptA.x + (1-s)*sqrt(1-ptA.y*ptA.y);
    ptB.y = ptA.y;
    ptB.z = sqrt(1 - ptB.x*ptB.x - ptB.y*ptB.y);
    //(c&d) compute rotated basis vectors.
    //We compute A^T * e directly so that we do not need to factor the S matrices into the A matrices
    Vector3f q = Vector3f(-wo.y,wo.x,0);  // q = l X n
    if ((wo.x==0)&&(wo.y==0)) {
        q = Vector3f(1,0,0);  //if wo==n then we can pick q to be any vector perpendicular to n
    }
    Float len_Av = sqrt(vAAv);
    Float len_ATq = sqrt(compute_vSv(Sdiagonals, Soffdiagonals, q));
    Vector3f ATe3 = compute_Sv(invSdiagonals, invSoffdiagonals, wo) / len_Av;
    Vector3f ATe2 = q / len_ATq;
    Vector3f ecross = Cross(ATe2, ATe3);
//    Float detS = compute_determinantS(Sdiagonals, Soffdiagonals);
    Vector3f ATe1 = compute_Sv(invSdiagonals, invSoffdiagonals, ecross)*sqrt(1/detInvS);
    // create the sampled microfacet normal direction
    Vector3f m = ptB.x*ATe1 + ptB.y*ATe2 + ptB.z*ATe3;
    m = Normalize(m);
    *wi = Reflect(wo, m);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

    // Compute PDF of _wi_ for microfacet reflection
    *pdf = Pdf(wo, *wi);
    Spectrum val = f(wo, *wi);
//    printf("weight = %g  cos %g  cos %g\n",val[0]*AbsCosTheta(*wi)/(*pdf),AbsCosTheta(wo),AbsCosTheta(*wi));
    return val;
}

//This interface is a bit unclear as to which direction was sampled wi or wo
//assuming it is wi that was sampled based on input wo, but may need to be updated for more general case
Float EllipsoidSReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;
    //compute BRDF value (which should always include normalization for this use)
    if (wi.z <= 0) return 0;
    if (wo.z <= 0) return 0;  //wrong side of surface
    //compute microfacect distribution value
    Vector3f h = wi + wo;   //unnormalized half vector
    Float hLen2 = h.LengthSquared();
    Float hSh = compute_vSv(Sdiagonals, Soffdiagonals, h);
    Float normalization = Pi * sqrt(detInvS * invSdiagonals.z);
    Float Dval = hLen2*hLen2 / (hSh*hSh * normalization);
    //also need mono-directional shadowing (but without the cosine factor and clamping)
    Float nAAn = invSdiagonals.z;
    Float vAAv = compute_vSv(invSdiagonals, invSoffdiagonals, wo);
    Float vAAn = compute_zdirSv(invSdiagonals, invSoffdiagonals, wo);
    Float luneRatio = 2*nAAn / (sqrt(vAAv*nAAn) + vAAn);
    return Dval*luneRatio/4;
}

// Material wrapper class that can include one lambertian lobe and up to two ellipsoid lobes
void EllipsoidSMaterial::ComputeScatteringFunctions(
      SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
      bool allowMultipleLobes) const {
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
    // Initialize diffuse component if present
    if (lambAlbedo) {
        Spectrum kd = lambAlbedo->Evaluate(*si).Clamp();
        if (!kd.IsBlack()) si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));
    }
    // Initialize ellipsoid component
    if (e0scale) {
        //TODO: figure out to disable interpolation in these lookups
        Spectrum scale = e0scale->Evaluate(*si).Clamp();
        Spectrum _Sdiag = e0Sdiag->Evaluate(*si);  //no clamping for matrix elements
        Spectrum _Soffdiag = e0Soffdiag->Evaluate(*si);
        //convert matrix elements to vectors instead of spectra
        Vector3f Sdiag(_Sdiag[0],_Sdiag[1],_Sdiag[2]);
        Vector3f Soffdiag(_Soffdiag[0],_Soffdiag[1],_Soffdiag[2]);
        si->bsdf->Add(ARENA_ALLOC(arena, EllipsoidSReflection)(scale,Sdiag,Soffdiag));
    }
    if (e1scale) {
        //TODO: figure out to disable interpolation in these lookups
        Spectrum scale = e1scale->Evaluate(*si).Clamp();
        Spectrum _Sdiag = e1Sdiag->Evaluate(*si);  //no clamping for matrix elements
        Spectrum _Soffdiag = e1Soffdiag->Evaluate(*si);
        //convert matrix elements to vectors instead of spectra
        Vector3f Sdiag(_Sdiag[0],_Sdiag[1],_Sdiag[2]);
        Vector3f Soffdiag(_Soffdiag[0],_Soffdiag[1],_Soffdiag[2]);
        si->bsdf->Add(ARENA_ALLOC(arena, EllipsoidSReflection)(scale,Sdiag,Soffdiag));
    }
    /*
    Spectrum ks = Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack()) {
        if (alphaX) {
            Float roughX = alphaX->Evaluate(*si);
            if (alphaY) {
                Float roughY = alphaY->Evaluate(*si);  //anisotropic GGX
                si->bsdf->Add(ARENA_ALLOC(arena, EllipsoidSReflection)(ks,roughX,roughY));
            } else { //isotropic GGX
                si->bsdf->Add(ARENA_ALLOC(arena, EllipsoidSReflection)(ks,roughX));
            }
        } else {
            printf("Error textures not yet implemented\n");
            throw std::runtime_error( "Ellipsoid textures not yet implemented" );
            //TODO: add loading of SMatrix from textures
        }
    }
     */
}

EllipsoidSMaterial *CreateEllipsoidSMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTextureOrNull("lambertian_albedo");
    std::shared_ptr<Texture<Spectrum>> e0scale =
        mp.GetSpectrumTextureOrNull("ellip0_scale");
    std::shared_ptr<Texture<Spectrum>> e0S_diag =
        mp.GetSpectrumTextureOrNull("ellip0_S_diagonals");
    std::shared_ptr<Texture<Spectrum>> e0S_offd =
        mp.GetSpectrumTextureOrNull("ellip0_S_offdiagonals");
    if (!e0scale) {
        e0S_diag = NULL;
        e0S_offd = NULL;
    } 
    std::shared_ptr<Texture<Spectrum>> e1scale =
        mp.GetSpectrumTextureOrNull("ellip1_scale");
    std::shared_ptr<Texture<Spectrum>> e1S_diag =
        mp.GetSpectrumTextureOrNull("ellip1_S_diagonals");
    std::shared_ptr<Texture<Spectrum>> e1S_offd =
        mp.GetSpectrumTextureOrNull("ellip1_S_offdiagonals");
    if (!e1scale) {
        e1S_diag = NULL;
        e1S_offd = NULL;
    }
     return new EllipsoidSMaterial(Kd, e0scale,e0S_diag,e0S_offd, e1scale,e1S_diag,e1S_offd);
}

} //namespace pbrt

