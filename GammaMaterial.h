#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FEBioMech/FETransIsoMooneyRivlin.h>
#include <FEBioMech/FEActiveContractionMaterial.h>
#include <FECore/FEPlotData.h>

/*
 * Opengamma Material Point
 * Adds field m_gamma to FEElasticMaterialPoint
 * in order to couple with Opengammas FastMonodomainSolver
 */ 
class GammaMaterialPoint 
	: public FEElasticMaterialPoint
{
public:
    	GammaMaterialPoint(FEMaterialPointData *mp = nullptr)
    	        : FEElasticMaterialPoint(mp), m_gamma(0) {}

    	double m_gamma; // coupling vatiable
};

/*
 * Opengamma Material
 * Behaves the same as FETransIsoMooneyRivlin Material
 * but uses GammaMaterialPoint instead of FEElasticMaterialPoint
 */
class GammaMaterial 
	: public FETransIsoMooneyRivlin
{
public:

    	GammaMaterial(FEModel *pfem) 
    	        : FETransIsoMooneyRivlin(pfem) {}

    	virtual FEMaterialPointData *CreateMaterialPointData() override;

    	// Recompute just the active fiber-stress contribution (same Cauchy-frame value
    	// GammaContraction adds into DevStress()), independent of the passive stress.
    	mat3ds GetActiveStress(FEMaterialPoint &mp);

    	DECLARE_FECORE_CLASS();
};

/*
 * Opengamma Active Contraction
 * Implements Opengammas active contraction model,
 * see https://opendihu.readthedocs.io/en/latest/settings/muscle_contraction_solver.html
 * Use in combination with GammaMaterial
 */
class GammaContraction 
	: public FEActiveContractionMaterial
{
public:
    	GammaContraction(FEModel *pfem)
    	        : FEActiveContractionMaterial(pfem) {}

    	virtual mat3ds ActiveStress(FEMaterialPoint &mp, const vec3d &a0) override;
    	virtual tens4ds ActiveStiffness(FEMaterialPoint &mp, const vec3d &a0) override;

    	double m_pmax; 				// maximum PK2 active stress
	double m_lamOpt; 			// 1.2 constant in Opengamma
	bool m_enableForceLengthRelation; 	// whether f(lam/lam_opt) should be multiplied

    	DECLARE_FECORE_CLASS();
};

/*
 * Plots the active fiber-stress contribution alone, pulled back to the PK2 (reference)
 * frame the same way FEBio's built-in "PK2 stress" var does, so the two are directly
 * comparable, and so this can be compared against Opendihu's active tension (T) output.
 */
class FEPlotActivePK2Stress
	: public FEPlotDomainData
{
public:
    	FEPlotActivePK2Stress(FEModel *pfem)
    	        : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) {}

    	bool Save(FEDomain &dom, FEDataStream &a) override;
};
