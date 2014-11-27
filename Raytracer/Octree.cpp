#include "Octree.h"

OctreeNode::OctreeNode() { }
OctreeNode::OctreeNode(std::vector<TriangleD> triangles) : Triangles(triangles) { }

Octree::Octree() { }
Octree::Octree(std::vector<TriangleD> triangles) : root(OctreeNode(triangles)) { }