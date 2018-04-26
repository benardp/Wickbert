#include "ParticleShaderCurvature.h"
#include "ParticleMaterial.h"

REGISTER_PARTICLESTUFF(ParticleShaderCurvature,"Shader:ParticleShaderCurvature");

ParticleShaderCurvature::ParticleShaderCurvature(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderCurvature"))
{
	imp_int = 0;
	material=0;
	kind = 0;
	scale = 0.001;
}

void ParticleShaderCurvature::attachAttributes()
{
	ParticleShader::attachAttributes();
	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
}

/// parameters
int ParticleShaderCurvature::qlen()
{
	return 2;
}

void ParticleShaderCurvature::getq(double *q)
{
	q[0] = kind;
	q[1] = scale;
}

void ParticleShaderCurvature::setq(double *q)
{
	kind = q[0];
	scale = q[1];
}

void ParticleShaderCurvature::qname(char **qn)
{
	qn[0] = "Curvature Kind";
	qn[1] = "Curvature Scale";
}

char *ParticleShaderCurvature::qtip(int i)
{
	switch (i) {
		default:
		case 0:
			return "Kind of curvature: 0=Mean, 1=Gaussian";
		case 1:
			return "Curvature to [0,1] color scaling factor";
	}
}

void ParticleShaderCurvature::drawShape(int i)
{
	if (!imp_int) return;
	if (!material) return;

	Implicit *imp = imp_int->getImplicit();

	if (!imp) return;

	double k;
	switch (kind) {
	case 0:
		k = imp->meanCurvature(position->getPosition(i));
		break;
	case 1:
		k = imp->gaussianCurvature(position->getPosition(i));
		break;
	}

	k *= scale;

	if (k > 1.0) k = 1.0;
	if (k < -1.0) k = -1.0;

	gmVector4 frontdiff = gmVector4(0.0,1.0,0.0,1.0);
	if (k < 0.0) {
		frontdiff[2] = -k;
	} else {
		frontdiff[0] = k;
	}
	material->setDiffuseFront(i,frontdiff);
	material->setDiffuseBack(i,0.5*frontdiff);
}
