/**
 * @file ParticleMaterial.cpp
 * @author Terry Fleury
 */

#include "ParticleMaterial.h"
#include "cleangl.h"

REGISTER_PARTICLESTUFF(ParticleMaterial,"Attribute:ParticleMaterial");

/**
 * Add particle material attribute to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleMaterial::ParticleMaterial(Particles *ps, const std::string& name)
  : ParticleAttribute(ps, name) 
{
	new PSParamgmVector4(this,&all_diffuseFront,gmVector4(1.0,0.0,0.0,1.0),
		"all_fdiff","global front diffuse","Diffuse reflection of all front faces");
	new PSParamgmVector4(this,&all_diffuseBack,gmVector4(0.0,0.0,1.0,1.0),
		"all_bdiff","global back diffuse","Diffuse reflection of all back faces");
	new PSParamgmVector4(this,&sel_diffuseFront,gmVector4(1.0,1.0,0.0,1.0),
		"sel_fdiff","selected front diffuse","Diffuse reflection of selected front faces");
	new PSParamgmVector4(this,&sel_diffuseBack,gmVector4(0.5,0.5,0.0,1.0),
		"sel_bdiff","selected back diffuse","Diffuse reflection of selected back faces");

	new PSParamgmVector4PerParticle(this,&diffuseFront,
		"fdiff","front diffuse","Diffuse reflection of front faces");
	new PSParamgmVector4PerParticle(this,&diffuseBack,
		"bdiff","back diffuse","Diffuse reflection of back faces");
}

void ParticleMaterial::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);

	if (ps)
	{
      diffuseFront.resize(ps->size());
      diffuseBack.resize(ps->size());
	  for (unsigned int i=0;i<ps->size();i++) {
		diffuseFront[i] = gmVector4(1.0,1.0,0.0,1.0);
		diffuseBack[i] = gmVector4(1.0,0.0,1.0,1.0);
	  }
	}
}

/**
 * Add a material diffuse property corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleMaterial::particleAdded() 
{
  // When adding a new particle, give it the default GL color
  diffuseFront.push_back(gmVector4(0.8,0.8,0.8,1.0));
  diffuseBack.push_back(gmVector4(0.8,0.8,0.8,1.0));
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleMaterial::particleRemoved(unsigned int i) 
{
  diffuseFront[i] = diffuseFront.back();
  diffuseBack[i] = diffuseBack.back();
  diffuseFront.pop_back();
  diffuseBack.pop_back();
}

/// Remove all particles from the vectors
void ParticleMaterial::clear()
{
  diffuseFront.clear();
  diffuseBack.clear();
}

void ParticleMaterial::setMaterial(unsigned int i)
{
	gmVector4 df_v4 = all_diffuseFront[3]*all_diffuseFront + (1.0 - all_diffuseFront[3])*diffuseFront[i];
	gmVector4 db_v4 = all_diffuseBack[3]*all_diffuseBack + (1.0 - all_diffuseBack[3])*diffuseBack[i];

	if (i == ps->selectedParticle) {
		df_v4 = sel_diffuseFront[3]*sel_diffuseFront + (1.0 - sel_diffuseFront[3])*diffuseFront[i];
		db_v4 = sel_diffuseBack[3]*sel_diffuseBack + (1.0 - sel_diffuseBack[3])*diffuseBack[i];
	}


	GLfloat df[4],db[4];

	df_v4.copyTo(df);
	db_v4.copyTo(db);

	glMaterialfv(GL_FRONT, GL_DIFFUSE, df);
	glMaterialfv(GL_BACK, GL_DIFFUSE, db);}

