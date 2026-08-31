#include "GammaMaterial.h"
#include <FECore/writeplot.h>

BEGIN_FECORE_CLASS(GammaMaterial, FETransIsoMooneyRivlin)
END_FECORE_CLASS();

// Create GammaMaterialPoint instead of FEElasticMaterialPoint
FEMaterialPointData *GammaMaterial::CreateMaterialPointData() {
    	auto pt = new GammaMaterialPoint;
    	if (m_ac) pt->Append(m_ac->CreateMaterialPointData());
    	return pt;
}

// Recomputes the same fiber direction FETransIsoMooneyRivlin::DevStress() uses, then
// calls the active_contraction property directly so the active-only contribution can
// be inspected on its own, without the passive Mooney-Rivlin/fiber stress mixed in.
mat3ds GammaMaterial::GetActiveStress(FEMaterialPoint &mp) {
    	if (m_ac == nullptr) return mat3ds(0.0);

    	mat3d Q = GetLocalCS(mp);
    	vec3d fiber = m_fiber->unitVector(mp);
    	vec3d a0 = Q*fiber;

    	return m_ac->ActiveStress(mp, a0);
}

bool FEPlotActivePK2Stress::Save(FEDomain &dom, FEDataStream &a) {
    	GammaMaterial *mat = dynamic_cast<GammaMaterial*>(dom.GetMaterial());
    	if (mat == nullptr) return false;

    	writeAverageElementValue<mat3ds>(dom, a, [&](const FEMaterialPoint &mp) {
    		FEMaterialPoint &mmp = const_cast<FEMaterialPoint&>(mp);
    		const FEElasticMaterialPoint &ep = *mp.ExtractData<FEElasticMaterialPoint>();
    		mat3ds s = mat->GetActiveStress(mmp);
    		return ep.pull_back(s);
    	});

    	return true;
}

BEGIN_FECORE_CLASS(GammaContraction, FEActiveContractionMaterial)
	ADD_PARAMETER(m_pmax, "pmax");
	ADD_PARAMETER(m_lamOpt, "lam_opt");
	ADD_PARAMETER(m_enableForceLengthRelation, "enable_force_length_relation");
END_FECORE_CLASS();

// Uses Opengammas active stress calculation
mat3ds GammaContraction::ActiveStress(FEMaterialPoint &mp, const vec3d &a0) {
	GammaMaterialPoint &pt = *mp.ExtractData<GammaMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// Calculate active stress using Opengammas formula
	double lamRelative = lamd/m_lamOpt;
	double f = 1.0;
	if (m_enableForceLengthRelation
 	    && 0.6 <= lamRelative && lamRelative <= 1.4) {
		f = -25.0/4.0 * lamRelative*lamRelative + 25.0/2.0 * lamRelative - 5.25;
	}
	double saf = 1.0/lamd * m_pmax * f * pt.m_gamma;

	return AxA*saf;
}

// Do not use active stiffness
tens4ds GammaContraction::ActiveStiffness(FEMaterialPoint &mp, const vec3d &a0) {
	return tens4ds(0.0);
}
