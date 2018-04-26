#include "Rainbow.h"
#include "ParticleMaterial.h"

REGISTER_PARTICLESTUFF(Rainbow,"Shader:Rainbow");

Rainbow::Rainbow(Particles *ps)
	:ParticleShader(ps,std::string("Rainbow"))
{
	new PSParamString(this,&source,"ImplicitInterrogator:kg","source","source",
		"Source to index into the rainbow colormap");
	new PSParamDouble(this,&red,-1.0,"red","red value",
		"Value at the low end of the spectrum");
	new PSParamDouble(this,&violet,1.0,"violet","violet value",
		"Value at the hight end of the spectrum");
	new PSParamBool(this,&wrap,false,"wrap","Cycle color map",
		"clamps color index when false; cycles color map when true.");

	new Attached<ParticleMaterial>(this,&material);

	source_data = NULL;
}

void Rainbow::attachAttributes()
{
	ParticleShader::attachAttributes();

	PSParamPerParticle *source_param = NULL;
	ParticleAttribute *source_attr = NULL;

	std::string attr_name = source.substr(0,source.find(':'));
	std::string param_name = source.substr(source.find(':')+1);
	if (source_attr = ps->getAttribute<ParticleAttribute>(attr_name))
		source_param = source_attr->perparticle.findparam(param_name);
	if (source_param)
		source_data = (DoubleVector *)source_param->ref();
	else
		source_data = NULL;
}

void Rainbow::drawParticle(int i)
{
	if (!source_data) return;

	double f = source_data->at(i);

	f = f - red;

	double w = violet - red;
	if (w == 0.0) return;
	f *= 6.0/w;

	if (wrap) {
		f = fmod(f,6.0);
		if (f < 0.0) f += 6.0;
	} else
		f = f < 0.0 ? 0.0 : f > 6.0 ? 6.0 : f;

	gmVector4 rgba;

	if (f < 1.0)
		rgba = gmVector4(1.0,f,0.0,1.0);
	else if (f < 2.0)
		rgba = gmVector4(2.0-f,1.0,0.0,1.0);
	else if (f < 3.0)
		rgba = gmVector4(0.0,1.0,f-2.0,1.0);
	else if (f < 4.0)
		rgba = gmVector4(0.0,4.0-f,1.0,1.0);
	else if (f < 5.0)
		rgba = gmVector4(f-4.0,0.0,1.0,1.0);
	else
		rgba = gmVector4(1.0,0.0,6.0-f,1.0);

	material->diffuseFront[i] = rgba;
}


