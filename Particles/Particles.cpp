/**
 * Implementation of Particles.
 * @file Particles.cpp
 * @author John Hart, Ed Bachta, Wen Su
 */

#include "Particles.h"
#include "ParticleSystem.h"
#include "ParticleAttribute.h"
#include "ParticleBehavior.h"
#include "ParticleShader.h"

#include <time.h>
#include <iostream>
#include <strstream>

//Actually to correctly query the duration 
//of open gl calls, one needs to call glFinish()
//this makes the program slower.
//therefore I added this "option" here.
#define MEASURE_REAL_DRAW_TIME
#ifdef MEASURE_REAL_DRAW_TIME
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#endif
using namespace std;


Particles::Particles(ParticleSystem *ps, const std::string& n)
	:name(n)
{
	dt = 0.03;
	particleSystem=ps;
	if (ps) ps->particles.push_back(this);
	particleSize=0;
	draggable=true;
	play = true;
	show = true;
	selectedParticle = -1;
}

bool
Particles::remove(ParticleStuff *stuff)
{
	if (ParticleAttribute *attr = dynamic_cast<ParticleAttribute *>(stuff)) {
		ParticleAttributes::iterator attr_it = attributes.begin();
		while (attr_it->second != attr && attr_it != attributes.end())
			attr_it++;
		if (attr_it == attributes.end()) return false;
		delete attr_it->second;
		attributes.erase(attr_it);
		return true;
	} else if (ParticleBehavior *beha = dynamic_cast<ParticleBehavior *>(stuff)) {
		ParticleBehaviors::iterator beha_it = find(behaviors.begin(), behaviors.end(), beha);
		if (beha_it == behaviors.end()) return false;
		delete *beha_it;
		behaviors.erase(beha_it);
		return true;
	} else if (ParticleShader *shad = dynamic_cast<ParticleShader *>(stuff)) {
		ParticleShaders::iterator shad_it = find(shaders.begin(), shaders.end(), shad);
		if (shad_it == shaders.end()) return false;
		delete *shad_it;
		shaders.erase(shad_it);
		return true;
	}
	return false;
}



/** Add a new particle to the system.
 *  This routine is called any time the # of particles is increased.
 *  It first increases the number of particles returned by size(), then
 *  calls all of the particleAdded() routines in the attributes, behaviors and shaders.
 *  The particle is always added at the end of the vectors, so its index should be size()-1.
 *  @returns index of the new particle.
 */
unsigned int Particles::addParticle()
{
	int newparticle = particleSize++;

	// Call attribute callbacks for particle addition
	// position must be called first
	for(ParticleAttributes::iterator pa_it = attributes.begin();pa_it != attributes.end();++pa_it)
		pa_it->second->particleAdded();
	
	// Call behavior callbacks for particle addition
	for(ParticleBehaviors::iterator b_it = behaviors.begin();b_it != behaviors.end();++b_it)
		(*b_it)->particleAdded();

	// Call behavior callbacks for particle addition
	for(ParticleShaders::iterator s_it = shaders.begin();s_it != shaders.end();++s_it)
		(*s_it)->particleAdded();

	return newparticle;
}

/**
 * Remove a particle from the system. We remove particle i by copying 
 * particle n's data into particle i's elements and then removing particle n. 
 * This change of index must be observed in all other objects which maintain 
 * references by index.
 * @param i Index of the particle to remove.
 */
void Particles::removeParticle(unsigned int i)
{
	// decrease the size order is important
	--particleSize;

	// Call attribute callbacks for particle removal
	for(ParticleAttributes::iterator pa_it = attributes.begin();pa_it != attributes.end();++pa_it)
		pa_it->second->particleRemoved(i);

	// Call behavior callbacks for particle removal
	for(ParticleBehaviors::iterator b_it = behaviors.begin();b_it != behaviors.end();++b_it)
		(*b_it)->particleRemoved(i);

	// Call shaders callbacks for particle removal
	for(ParticleShaders::iterator s_it = shaders.begin();s_it != shaders.end();++s_it)
		(*s_it)->particleRemoved(i);

}

//! Destructor deletes all behaviors and attribute.
Particles::~Particles()
{
	unsigned int i;
	for(ParticleAttributes::iterator pa_it = attributes.begin();pa_it != attributes.end();++pa_it)
		delete pa_it->second;
	for(i=0;i<behaviors.size();i++)
		delete behaviors[i];
	for(i=0;i<shaders.size();i++)
		delete shaders[i];
}

/*! Searches for the behavior based on its name.
	\return	The index of the behavior in the vector.
*/
int Particles::findBehavior(const std::string &name)
{
	for (unsigned int i = 0; i < behaviors.size(); i++)
	{
		if (behaviors[i]->name == name)
			return i;
	}
	return -1;
}

/// Searches for the shader based on its name.
int Particles::findShader(const std::string &name)
{
	for (unsigned int i = 0; i < shaders.size(); i++)
	{
		if (shaders[i]->name == name)
			return i;
	}

	return -1;
}

void Particles::fullUpdate()
{
	if (!play)
		return;
	prepareAttributes();
	applyForces();
	applyConstraints();
	integrate();
	cleanup();
}

void Particles::prepareAttributes()
{
	time_t t;
	for (ParticleAttributes::iterator it = attributes.begin(); it != attributes.end(); ++it)
	{
		t = clock();
		it->second->prepare();
		t = clock() - t;
		it->second->total_time += (double) t / CLOCKS_PER_SEC;
		++(it->second->iterations);
	}
}

void Particles::applyForces()
{
	time_t t;
	for(unsigned int i=0;i<behaviors.size();++i) {
		t = clock();
		behaviors[i]->applyForce();
		t = clock() - t;
		behaviors[i]->total_time += (double) t / CLOCKS_PER_SEC;
		++(behaviors[i]->iterations);
	}
}

void Particles::applyConstraints()
{
	time_t t;
	for(unsigned int i = 0; i < behaviors.size(); i++) {
		t = clock();
		behaviors[i]->applyConstraint();
		t = clock() - t;
		behaviors[i]->total_time += (double) t / CLOCKS_PER_SEC;
		++(behaviors[i]->iterations);
	}
}

void Particles::integrate()
{
	time_t t;
	for(unsigned int i=0;i<behaviors.size();i++) {
		t = clock();
		behaviors[i]->integrate();
		t = clock() - t;
		behaviors[i]->total_time += (double) t / CLOCKS_PER_SEC;
		++(behaviors[i]->iterations);
	}
}

void Particles::cleanup()
{
	time_t t;
	for(unsigned int i=0;i<behaviors.size();i++) {
		t = clock();
		behaviors[i]->cleanup();
		t = clock() - t;
		behaviors[i]->total_time += (double) t / CLOCKS_PER_SEC;
		++(behaviors[i]->iterations);
	}
}

void Particles::attachAttributes() {
	ParticleAttributes::iterator attr_it;
	for (attr_it = attributes.begin(); attr_it != attributes.end(); attr_it++)
		attr_it->second->attachAttributes();
	ParticleBehaviors::iterator beha_it;
	for (beha_it = behaviors.begin(); beha_it != behaviors.end(); beha_it++)
		(*beha_it)->attachAttributes();
	ParticleShaders::iterator shad_it;
	for (shad_it = shaders.begin(); shad_it != shaders.end(); shad_it++)
		(*shad_it)->attachAttributes();
}

/**
 * Prints the contents of this particle system.
 */
void Particles::printContent() {
	std::cout << "-- Particle System Contents --" << std::endl;
	std::cout << "Attributes:" << std::endl;
	for(ParticleAttributes::iterator a_it = attributes.begin();
		a_it != attributes.end();
		a_it++) 
	a_it->second->printContent();
	std::cout << "Behaviors:" << std::endl;
	for(unsigned int j=0;j<behaviors.size();j++)
		behaviors[j]->printContent();
} 

void Particles::draw()
{
	if (!show)
		return;

	unsigned int s,i;
	time_t t;


#ifdef MEASURE_REAL_DRAW_TIME
		glFinish();
#endif
	for(s = 0; s < shaders.size(); s++) {
		glPushName(GLuint(-1));
		glPushMatrix();
		t = clock();
		shaders[s]->drawPre();
#ifdef MEASURE_REAL_DRAW_TIME
		glFinish();
#endif
		t = clock() - t;
		shaders[s]->pre_time += (double) t / CLOCKS_PER_SEC;
		shaders[s]->iterations++;
		glPopMatrix();
		glPopName();
	}

	for (i = 0; i < size(); i++) {
		glPushName(i);
		glPushMatrix();
		for(s = 0; s < shaders.size(); s++) {
			t = clock();
			shaders[s]->drawParticle(i);
#ifdef MEASURE_REAL_DRAW_TIME
			glFinish();
#endif
			t = clock() - t;
			shaders[s]->draw_time += (double) t / CLOCKS_PER_SEC;
			shaders[s]->iterations++;
		}
		glPopMatrix();
		glPopName();
	}

	for(s = 0; s < shaders.size(); s++) {
		glPushName(GLuint(-1));
		glPushMatrix();
		t = clock();
		shaders[s]->drawPost();
#ifdef MEASURE_REAL_DRAW_TIME
		glFinish();
#endif
		t = clock() - t;
		shaders[s]->post_time += (double) t / CLOCKS_PER_SEC;
		shaders[s]->iterations++;
		glPopMatrix();
		glPopName();
	}
}

void Particles::event(int e)
{
	for(unsigned int i=0;i<shaders.size();i++)
	{
		shaders[i]->event(e);
	}
}

void Particles::error(std::string e)
{
	particleSystem->error("[" + name + "] " + e);
}

void Particles::copyFrom(Particles *p)
{
#if 0
	/** Copy attributes.
	 */
	ParticleAttributes::iterator attr_it;
	for (attr_it = p->attributes.begin(); attr_it != p->attributes.end(); attr_it++) {
		ParticleAttribute *attr = (ParticleAttribute *)(NEW_PARTICLESTUFF(attr_it->second->getClass()).release());
		if (attr) {
			*attr = *(attr_it->second);
			attr->setParticleSystem(this);
		}
	}

	/** Copy behaviors.
	 */
	ParticleBehaviors::iterator beha_it;
	for (beha_it = p->behaviors.begin(); beha_it != p->behaviors.end(); beha_it++) {
		ParticleBehavior *beha = (ParticleBehavior *)(NEW_PARTICLESTUFF((*beha_it)->getClass()).release());
		if (beha) {
			*beha = **beha_it;
			beha->setParticleSystem(this);
		}
	}

	/** Copy shaders.
	 */
	ParticleShaders::iterator shad_it;
	for (shad_it = p->shaders.begin(); shad_it != p->shaders.end(); shad_it++) {
		ParticleShader *shad = (ParticleShader *)(NEW_PARTICLESTUFF((*shad_it)->getClass()).release());
		if (shad) {
			*shad = **shad_it;
			shad->setParticleSystem(this);
		}
	}
#else
	std::string preserve_name = name;
	std::strstream buf;
	buf << p;
	buf >> this;
	name = preserve_name;

	attachAttributes();
#endif
}

std::string Particles::profile()
{
	std::string report = "Particles: " + name + "\n";
	
	char buf[32];
	
	sprintf(buf,"%5d",this->size());
	report += "Number of particles: " + std::string(buf) + "\n";

	double totalTimeForClass=0;


	for (ParticleAttributes::iterator it = attributes.begin(); it != attributes.end(); ++it)
	{
		if (it->second->iterations>0)
		{
			double time=1000.0*(it->second->total_time)/(it->second->iterations);
			if (time>0)
			{ 
				totalTimeForClass+=time;
				sprintf(buf,"%6.3f",time);
				report += std::string(buf) + " ms in preparing attribute:" + it->first+ "\n";
				it->second->iterations=0;
				it->second->total_time=0;
			}
		}
	}
	sprintf(buf,"%6.3f",totalTimeForClass);
	report += std::string("________________________")+ "\n";
	report += std::string(buf) + " ms in attributes" + "\n";
	report += std::string("________________________")+ "\n";
	report += std::string("________________________")+ "\n";

	totalTimeForClass=0;
	for (unsigned int i = 0; i < behaviors.size(); i++) {
		if (behaviors[i]->iterations > 0) {
			double time=1000.0*behaviors[i]->total_time/behaviors[i]->iterations;
			totalTimeForClass+=time;
			sprintf(buf,"%6.3f",time);
			report += std::string(buf) + " ms in behavior " + behaviors[i]->name + "\n";
		}
		behaviors[i]->iterations=0;
		behaviors[i]->total_time=0;
	}
	
	sprintf(buf,"%6.3f",totalTimeForClass);
	report += std::string("________________________")+ "\n";
	report += std::string(buf) + " ms in behaviors " + "\n";
	report += std::string("________________________")+ "\n";
	report += std::string("________________________")+ "\n";

	totalTimeForClass=0;
	for (unsigned int i = 0; i < shaders.size(); i++) {
		if (shaders[i]->iterations > 0) {
			double time=1000.0*shaders[i]->pre_time/shaders[i]->iterations;
			totalTimeForClass+=time;
			sprintf(buf,"%6.3f",time);
			report += std::string(buf) + " ms in behavior (pre) " + shaders[i]->name + "\n";
			time=1000.0*shaders[i]->draw_time/shaders[i]->iterations;
			totalTimeForClass+=time;
			sprintf(buf,"%6.3f",time);
			report += std::string(buf) + " ms in behavior (draw) " + shaders[i]->name + "\n";
			time=1000.0*shaders[i]->post_time/shaders[i]->iterations;
			totalTimeForClass+=time;
			sprintf(buf,"%6.3f",time);
			report += std::string(buf) + " ms in behavior (post) " + shaders[i]->name + "\n";
		}
		shaders[i]->draw_time=0;
		shaders[i]->pre_time=0;
		shaders[i]->post_time=0;
		shaders[i]->iterations=0;	
	}
	
	report += std::string("________________________")+ "\n";
	sprintf(buf,"%6.3f",totalTimeForClass);
	report += std::string(buf) + " in shaders" +"\n";
	report += std::string("________________________")+ "\n";
	report += std::string("________________________")+ "\n";


	return report;
}

void *
Particles::getVariable(std::string ref)
{
	unsigned int colon = ref.find(':');
	if (colon == ref.length())
		return NULL;
	std::string attrname = ref.substr(0,colon);
	ref = ref.substr(colon+1);

	ParticleAttribute *attr = getAttributeGeneric(attrname);
	if (!attr) return NULL;

	return attr->getVariable(ref);
}

/** Outputs a Particles of the form "classname [name] { <parameters> } \n".
 */
std::ostream &operator<<(std::ostream &out, Particles *p) {
	unsigned int i;

	out << "\t" << '"' << p->name << '"' << " {" << std::endl;
	for (ParticleAttributes::iterator attr_it = p->attributes.begin(); attr_it != p->attributes.end(); attr_it++)
		out << attr_it->second;
	for (i = 0; i < p->behaviors.size(); i++)
		out << p->behaviors[i];
	for (i = 0; i < p->shaders.size(); i++)
		out << p->shaders[i];
	out << "\t}" << endl;

	return out;
}

/** Removes all of the attributes.
*/
void Particles::removeAttributes() {
	ParticleAttributes::iterator attr_it;
	for (attr_it = attributes.begin(); attr_it != attributes.end(); attr_it++)
		delete attr_it->second;
	attributes.clear();
}

/** Removes all of the behaviors.
*/
void Particles::removeBehaviors() {
	ParticleBehaviors::iterator beha_it;
	for (beha_it = behaviors.begin(); beha_it != behaviors.end(); beha_it++)
		delete *beha_it;
	behaviors.clear();
}

/** Removes all of the shaders.
*/
void Particles::removeShaders() {
	ParticleShaders::iterator shad_it;
	for (shad_it = shaders.begin(); shad_it != shaders.end(); shad_it++)
		delete *shad_it;
	shaders.clear();
}

void Particles::removeStuff()
{
	removeAttributes();
	removeBehaviors();
	removeShaders();
}

/** Reads in several particlestuff entries of the form:
 *    classname ["name"] { params }
 * but only reads the classname, the optional, quoted name, and the opening bracket
 * of the params list
 */
std::istream &operator>>(std::istream &in, Particles *p) {
	std::string token;
	char c;

	if ((in >> c) && c == '"') {
		// we have a quoted name for the particles, read it in
		token = "";
		in >> std::noskipws;	// read in spaces too
		while ((in >> c) && c != '"') token += c;
		p->name = token;
	}

	in >> std::skipws;	// skip over whitespace from now on
	if (!(in >> token) || token != "{") {
		p->error("Opening bracket missing in Particles " + p->name);
		return in;
	}

	while ((in >> token) && token != "}") {
		// read in particle stuff
		std::auto_ptr<ParticleStuff> auto_stuff = NEW_PARTICLESTUFF(token);
		ParticleStuff *stuff = auto_stuff.release();
		if (!stuff) {
			p->error("Cannot create ParticleStuff \"" + token + "\". ");
			// should really eat up the parameter list
			//Exactly ;)
			//so why not just do it???? 
			//the disappearence of polylines leads to a crash otherwise... - Elmar
			std::string s;
			std::getline(in,s);
			continue;
		}
		stuff->setParticleSystem(p);
		if (!(in >> c)) continue;
		if (c == '"') {
			// quoted stuff name
			token = "";
			in >> std::noskipws;
			while ((in >> c) && c != '"')
				token += c;
			in >> std::skipws;
			stuff->setName(token);
			if (!(in >> c)) continue;	// should be the "{" or "["
		}
		if (c == '[') {
			in >> stuff->attachedattributes;
			in >> c;	// eat the '{'
		}
		// now read parameters
		in >> stuff->params;
	}

	if (token != "}") {
		p->error("Missing } at end of Particles \"" + p->name + "\".");
	}

	return in;
}

// Grab the attributes to add them to factory
extern int
	grabAdaptiveRepulsionData,
	grabFaceOrientation,
	grabFacePosition,
	grabImplicitInterrogator,
	grabLightPosition,
	grabMeshInterrogator,
	grabParticleAcceleration,
	grabParticleAge,
	grabParticleBoundingBox,
	grabParticleDensity,
	grabParticleDesiredOrientation,
	grabParticleEigen,
	grabParticleLocality,
	grabParticleLocalityGrid,
	grabParticleMaterial,
	grabParticleMesh,
	grabParticleOrientation,
	grabParticleNormal,
	grabParticlePosition,
	grabParticleScalar,
	grabParticleVelocity,
	grabParticleVisibility,
	grabSurfaceInterrogator,
	grabVertexPosition,
	grabVertexOrientation,
	grabOBJPosition,
	grabViewDependence,
	grabVertexScalar,
	//Elmar
	grabFindFirstIntersectionWithImplicits,
	grabDuplicateParticlesAlongNormal,
	grabUniformDistribution,
	grabImplicitInterrogatorCached,
	grabShaderSurfaceInterrogation,
	grabParticleVector,
	grabContours;


// Grab the behaviors to add them to factory
extern int
	grabParticleCreation,
	grabParticleFate,
	grabParticleFateContour,
	grabParticleValueDeath,
	grabParticleRepulsion,
	grabAsymmetricRepulsion,
	grabSurfaceAdhesion,
	grabSurfaceDeformation,
	grabSurfacePropagation,
	grabKDPropagation,
	grabViscousDynamics,
	grabAccelerationIntegrator,
	grabFeatureDetector,
	grabSilhouetteDetector,
	grabShadowDetector,
	grabSuggestiveContourAdhesion,
	grabDeCarloSuggestiveContourAdhesion,
	grabSuggestiveContourStabilization,
	grabShadowAdhesion,
	grabSpecularAdhesion,
	grabSilhouetteAdhesion,
	grabCurveAdhesion,
//	grabSkinMeshShape,
    grabKeepInBounds,
	grabParabolicsAdhesion,
	grabSingularityAdhesion,
	grabParticleRepulsionContour,
	grabSingularityRepulsion;
#if 0	//BONUS to be fixed
	grabVelocityIntegrator,
	grabSPHConstraint,
	grabSilhouetteAdhesion,
	grabSilhouetteFate,
	grabTopologyInterrogator,
	grabMeshShape,
	grabCriticalPointInterrogator,
	grabClusterMeshInterrogator,
	grabVectorField,
	grabRBFPropagation,
	grabVisibleSurfaceAdhesion,
	grabNeighborAlignment,
	grabSPHInitial,
	grabSPHGravity,
	grabSPHPressure,
    grabSPHViscosity,
	grabSPHTension;
#endif

// Grab the shaders to add them to factory
extern int
	grabParticleShaderDisk,
	//grabParticleShaderDiskDevelop,
	grabParticleShaderCylinder,
	grabParticleShaderSphere,
	grabUseMaterial,
	grabUseMaterialSwitch,
	grabCopyParticle,
	grabParticleShaderTriangle,
	grabParticleShaderText,
	grabOrientParticle,
	grabParticleShaderLink,
	grabShowCell,
	grabRainbow,
	grabIntervalSubdivision,
	grabShaderBoundingBox,
#ifdef WB_USE_CG 
	grabParticleShaderVoronoi,
	grabParticleShaderVoronoiDevelop,
	grabParticleShaderVoronoiPhong,
#endif
	grabShaderSuggestiveContour,
	grabShaderDeCarloSuggestiveContour,
	grabShaderSilhouetteContour,	
	grabShaderSilhouetteSimple,
	grabShaderIntersectionContour,
	grabShaderShadowContour,
	grabShaderSpecularContour,
	grabShaderParabolicsContour,
	grabShaderAnimation,
	grabShaderProcessSurfaceEvent,
	grabShaderMeshRenderer;

#if 0	//BONUS to be fixed
	grabParticleShaderChain,
	grabParticleShaderStroke,
	grabParticleShaderGradientField,
	grabParticleShaderLink,
	grabParticleShaderCluster,
	grabShaderVectorField,
	grabParticleShaderCurvature,
	grabParticleShaderRadialCurvature,
	grabParticleShaderContourHider,
	grabParticleShaderViewpointStability,
	grabParticleShaderOriented,
#endif


// Dummy function to activate automatic factory code.
void grab_particlestuff()
{
	// attributes
	grabParticleBoundingBox++;
	grabAdaptiveRepulsionData++;
	grabImplicitInterrogator++;
	grabParticlePosition++;
	grabParticleVelocity++;
	grabParticleLocality++;
	grabParticleLocalityGrid++;
	grabParticleOrientation++;
	grabParticleNormal++;
	grabParticleMaterial++;
	grabParticleAcceleration++;
	grabParticleDensity++;
	grabParticleAge++;
	grabParticleMesh++;
	grabParticleDesiredOrientation++;
	grabParticleEigen++;
	grabParticleScalar++;
	grabSurfaceInterrogator++;
	grabViewDependence++;
	grabParticleVisibility++;
	grabMeshInterrogator++;
	grabVertexPosition++;
	grabVertexOrientation++;
	grabOBJPosition++;
	grabFacePosition++;
	grabFaceOrientation++;
	grabLightPosition++;
	grabVertexScalar++;
	//ELMAR
	grabFindFirstIntersectionWithImplicits++;
	grabDuplicateParticlesAlongNormal++;
	grabUniformDistribution++;
	grabShaderSurfaceInterrogation++;
	grabImplicitInterrogatorCached++;
	grabParticleVector++;
	grabContours++;
	
	// behaviors
	grabParticleFate++;
	grabParticleValueDeath++;
	grabParticleRepulsion++;
	grabAsymmetricRepulsion++;
	grabSurfaceAdhesion++;
	grabSurfaceDeformation++;
	grabSurfacePropagation++;
	grabKDPropagation++;
	grabParticleCreation++;
	grabViscousDynamics++;
	grabAccelerationIntegrator++;
	grabFeatureDetector++;
	grabSilhouetteDetector++;
	grabShadowDetector++;
	grabSilhouetteAdhesion++;
	grabShadowAdhesion++;
	grabSpecularAdhesion++;
	grabSuggestiveContourAdhesion++;
	grabDeCarloSuggestiveContourAdhesion++;
	grabSuggestiveContourStabilization++;
	grabSingularityRepulsion++;
	grabCurveAdhesion++;
//	grabSkinMeshShape++;
	grabKeepInBounds++;
	grabParabolicsAdhesion++;
	grabSingularityAdhesion++;
	grabParticleRepulsionContour++;

#if 0	//BONUS to be fixed
	grabVelocityIntegrator++;
	grabSPHConstraint++;
	grabSilhouetteFate++;
	grabTopologyInterrogator++;
	grabMeshShape++;
	grabCriticalPointInterrogator++;
	grabClusterMeshInterrogator++;
	grabRBFPropagation++;
	grabVectorField++;
	grabVisibleSurfaceAdhesion++;
	grabNeighborAlignment++;
	grabSPHInitial++;
	grabSPHGravity++;
	grabSPHPressure++;
	grabSPHViscosity++;
	grabSPHTension++;
#endif

	// shaders
	grabParticleShaderDisk++;
	//grabParticleShaderDiskDevelop++;
	grabParticleShaderCylinder++;
	grabParticleShaderSphere++;

	grabCopyParticle++;
	grabUseMaterial++;
	grabUseMaterialSwitch++;
	grabOrientParticle++;
	grabParticleShaderText++;
	grabParticleShaderTriangle++;
	grabShowCell++;
	grabParticleShaderLink++;
	grabShaderSilhouetteContour++;
	grabShaderSilhouetteSimple++;
	grabShaderIntersectionContour++;
	grabShaderShadowContour++;
	grabShaderSuggestiveContour++;
	grabShaderDeCarloSuggestiveContour++;
	grabShaderParabolicsContour++;
	grabShaderSpecularContour++;
	grabShaderAnimation++;
	grabShaderProcessSurfaceEvent++;
	grabRainbow++;
	grabIntervalSubdivision++;
	grabShaderBoundingBox++;
	grabShaderMeshRenderer++;
	
#ifdef WB_USE_CG
	grabParticleShaderVoronoi++;
	//grabParticleShaderVoronoiDevelop++;
	grabParticleShaderVoronoiPhong++;
#endif

#if 0	//BONUS to be fixed
	grabParticleShaderLink++;
	grabParticleShaderChain++;
	grabParticleShaderStroke++;
	grabParticleShaderGradientField++;
	grabParticleShaderCluster++;
	grabParticleShaderCurvature++;
	grabShaderVectorField++;
	grabParticleShaderRadialCurvature++;
	grabParticleShaderContourHider++;
	grabParticleShaderViewpointStability++;
	grabParticleShaderOriented++;
#endif
}
