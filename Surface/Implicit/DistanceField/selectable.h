#ifndef SELECTABLE_H
#define SELECTABLE_H

#include "boundingvolume.h"
#include "renderable.h"

class Selectable : public Renderable {
 private:
  bool selected;
 public:

  Selectable() : selected(false) {}

  virtual bool Selected() const { return selected; }
  virtual void Select(int button = -1) { selected = true; }
  virtual void DeSelect(int button = -1) { selected = false; }
  virtual void Toggle(int button = -1) { selected = !selected; }
  virtual BoundingSphere GetBoundingSphere() = 0;
  virtual scalar DistanceTo(const Vector3d&) = 0;
  virtual void Move(const Vector3d&) {}
};

#endif
