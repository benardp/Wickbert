#include "ParticleBoundingBox.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(ParticleBoundingBox,"Attribute:ParticleBoundingBox");

ParticleBoundingBox::ParticleBoundingBox(Particles* ps, const std::string& name)
:ParticleAttribute(ps, name)
{
	new PSParamgmVector3(this,&min,gmVector3(-10,-10,-10),
		"min","min corner","Coordinates of the lower-left-back corner of the bounding box");
	new PSParamgmVector3(this,&max,gmVector3(10,10,10),
		"max","max corner","Coordinates of the upper-right-front corner of the bounding box");
	new PSParamButton(this,new ParticleBoundingBoxTighten(this),"tighten","tighten",
		"Resize the bounding box to tightly bound the current particle positions.");
	new Attached<ParticlePosition>(this,&position);
}

bool ParticleBoundingBox::inBounds(const gmVector3 &p) const
{
	return
		(min[0] <= p[0] && p[0] <= max[0]) &&
		(min[1] <= p[1] && p[1] <= max[1]) &&
		(min[2] <= p[2] && p[2] <= max[2]);
}

void ParticleBoundingBox::tighten()
{
	computeBoundingBox(min,max);
}

void ParticleBoundingBox::computeBoundingBox(gmVector3& _min, gmVector3& _max) const
{
	std::vector<double> X,Y,Z;
	gmVector3 posi;
	for(unsigned int i = 0; i<ps->size();  i++ ){
		posi = position->getPosition(i);
		X.push_back(posi[0]);
		Y.push_back(posi[1]);
		Z.push_back(posi[2]);
	}
	_min[0] = *min_element(X.begin(),X.end());
	_min[1] = *min_element(Y.begin(),Y.end());
	_min[2] = *min_element(Z.begin(),Z.end());
	_max[0] = *max_element(X.begin(),X.end());
	_max[1] = *max_element(Y.begin(),Y.end());
	_max[2] = *max_element(Z.begin(),Z.end());
}

void ParticleBoundingBox::computeIntersection(const int i, 
											  const gmVector3& myCameraPosition,
											  const double bboxOffset,
											  gmVector3& intersectBbox ) const
{
	gmVector3 l(min);
	gmVector3 h(max);
	gmVector3 posi(position->getPosition(i));
	gmVector3 viewi( myCameraPosition - posi);
	
	viewi.normalize();
	double t[3] = {gmGOOGOL, gmGOOGOL, gmGOOGOL};
	double tmin = gmGOOGOL;
	for(unsigned int j = 0; j < 3; j++){
		if(gmIsZero(viewi[j])) continue;
		t[j] = (l[j] - posi[j])/viewi[j];
		if(t[j] > 0) continue;
		t[j] = (h[j] - posi[j])/viewi[j];
	}
	for(unsigned int j = 0; j < 3; j++){
		if(tmin > t[j]) tmin = t[j];
	}
	for(unsigned int j = 0; j < 3; j++){
		intersectBbox[j] = posi[j] +  bboxOffset * tmin * viewi[j] ;
	}
}

void  ParticleBoundingBox::computeIntersection(const gmVector3& P, const gmVector3& myCameraPosition, const double bboxOffset, gmVector3& v1, gmVector3& v2) const
{
	gmVector3 l(min);
	gmVector3 h(max);

	gmVector3 view(P - myCameraPosition);	
	view.normalize();
	double t1, t2;
	
	double t;
	gmVector3 intersect;
	for(unsigned int j = 0; j < 3; j++){
		if(gmIsZero(view[j])) continue;
		t = (l[j] - myCameraPosition[j])/view[j];
		for(unsigned int j = 0; j < 3; j++){
			intersect[j] = myCameraPosition[j] +  t * view[j] ;
		}
		if(inBounds(intersect)){
			v1 = intersect;
			t1 = t;
			break;
		}
	}
	for(unsigned int j = 0; j < 3; j++){
		if(gmIsZero(view[j])) continue;
		t = (h[j] - myCameraPosition[j])/view[j];
		for(unsigned int j = 0; j < 3; j++){
			intersect[j] = myCameraPosition[j] +  t * view[j] ;
		}
		if(inBounds(intersect)){
			v2 = intersect;
			t2 = t;
			break;
		}
	}
	if(t2 < t1) std::swap(v1, v2);
	
}
