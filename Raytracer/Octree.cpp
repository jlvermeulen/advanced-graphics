#include "Octree.h"

OctreeNode::OctreeNode() { }
OctreeNode::OctreeNode(TriangleD* triangles, int nTriangles) : Triangles(triangles), NTriangles(nTriangles) { }

Octree::Octree() { }
Octree::Octree(TriangleD* triangles, int nTriangles)
{

}