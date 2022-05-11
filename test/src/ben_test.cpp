#include <gtest/gtest.h>
#include <iostream>

#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/transfer_functions.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "load_test_meshes.h"

class BenTestSuite : public MeshAssetSuite {};

using namespace geometrycentral;
using namespace geometrycentral::surface;

geometrycentral::surface::VertexData<geometrycentral::Vector2> getlogmap(SurfaceMesh& mesh,
                                                                         IntrinsicGeometryInterface& geometry) {
  VectorHeatMethodSolver solver(geometry);
  return solver.computeLogMap(mesh.vertex(0));
}

void setlogmap(geometrycentral::surface::VertexData<geometrycentral::Vector2> rawvalues,
               polyscope::SurfaceMesh* polymesh, const std::string& name = "log map") {
  std::vector<double> values;
  for (unsigned int i = 0; i < rawvalues.size(); i++) {
    values.push_back(rawvalues[i].norm());
  }
  polymesh->addVertexScalarQuantity(name, values);
}


TEST_F(BenTestSuite, Fox) {
  const auto& asset = getAsset("fox.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& origGeometry = *asset.geometry;

  SignpostIntrinsicTriangulation tri(mesh, origGeometry);

  tri.flipToDelaunay();
  tri.delaunayRefine(35, +INFINITY, 100000);

  // auto& subdiv = tri.getCommonSubdivision();
  // auto newmesh = *subdiv.buildSimpleMesh();

  // subdiv.constructMesh();
  // ManifoldSurfaceMesh& finemesh = *subdiv.mesh;
  // VertexData<Vector3> csPositions = subdiv.interpolateAcrossA(origGeometry.vertexPositions);
  // VertexPositionGeometry newGeometry(finemesh, csPositions);

  polyscope::init(); // initialize the gui

  /* */
  auto origpolymesh = polyscope::registerSurfaceMesh("my mesh", origGeometry.vertexPositions, mesh.getFaceVertexList());
  // auto newpolymesh =
  //     polyscope::registerSurfaceMesh("my fine mesh", newGeometry.vertexPositions, finemesh.getFaceVertexList());

  setlogmap(getlogmap(mesh, origGeometry), origpolymesh);
  setlogmap(tri.restrictToInput(getlogmap(tri.mesh, tri)), origpolymesh, "intrinsic log map");
  // setlogmap(getlogmap(finemesh, newGeometry), newpolymesh, "extrinsic log map");

  // exact
  auto mything = exactGeodesicDistance(mesh, origGeometry, mesh.vertex(0));
  origpolymesh->addVertexScalarQuantity("exact log map", mything);
  /* */

  polyscope::show(); // pass control to the gui until the user exits

  // EXPECT_TRUE(tri.isDelaunay());
}