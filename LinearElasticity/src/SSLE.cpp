//
// LE created by chenghau
//
#include "util.h"
#include "LEInputData.h"
#include "LENodeData.h"
#include "LEBCSetup.h"
#include "SSLEEquation.h"
#include <petscvec.h>
#include <IMGALoop.h>
#include "GetSurfaceGp.h"
#include "GetTrueSurfaceGP.h"
#include "CalcStress.h"


#pragma mark OptSug
#include <CalcError.h>
//#include "CheckSurface.h"
#include "DACoarse.h"
#include "DARefine.h"

using namespace PETSc;

int main(int argc, char *argv[])
{

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  LEInputData inputData;

  if (!(inputData.ReadFromFile()))
  {
    if (!rank)
    {
      throw TALYFEMLIB::TALYException() << "Can't read the config file \n";
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  inputData.PrintInputData();
  ///------------------------------ Command line option to restart from a checkpoint -------------------------------////
  bool resume_from_checkpoint = false;
  {
    PetscBool resume = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-resume_from_checkpoint", &resume, nullptr);
    resume_from_checkpoint = (resume == PETSC_TRUE);
  }
  Checkpointer checkpointer(inputData.CheckpointNumbackup, "CheckPoint");
  bool restart_sum = false;
  {
    PetscBool restart_sum_ = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-restart_sum", &restart_sum_, nullptr);
    restart_sum = (restart_sum_ == PETSC_TRUE);
  }
  /// --------------------------------------------------------------------------------------------------------------////

#ifndef PROFILING
  /// this is for quick testing through different refine lvl
  if (argc == 2 and inputData.BaselvlFromArgument)
  {
    inputData.mesh_def.refine_lvl_base = std::atoi(argv[1]);
    TALYFEMLIB::PrintStatus("---------------------------------------");
    TALYFEMLIB::PrintStatus("Refine level base from argument = ", inputData.mesh_def.refine_lvl_base);
    TALYFEMLIB::PrintStatus("---------------------------------------");
  }
  /// this is for quick testing through different refine lvl and lambda
  if (argc == 3 and inputData.BaselvlFromArgument)
  {
    inputData.mesh_def.refine_lvl_base = std::atoi(argv[1]);
    inputData.RatioGPSBM = (double)std::atoi(argv[2]) / 100;
    TALYFEMLIB::PrintStatus("---------------------------------------");
    TALYFEMLIB::PrintStatus("Refine level base from argument = ", inputData.mesh_def.refine_lvl_base);
    TALYFEMLIB::PrintStatus("Lambda from argument = ", inputData.RatioGPSBM);
    TALYFEMLIB::PrintStatus("---------------------------------------");
  }
#else
  TALYFEMLIB::PrintStatus("---------------------------------------");
  TALYFEMLIB::PrintStatus("please put -log_view after the run to print out the time");
  TALYFEMLIB::PrintStatus("---------------------------------------");
#endif

  const DENDRITE_UINT eleOrder = inputData.elemOrder;
  const DENDRITE_UINT levelBase = inputData.mesh_def.refine_lvl_base;
  const bool mfree = inputData.ifMatrixFree;

// Linear elasticity
#if (DIM == 2)
  static const char *varname[]{"UX", "UY"};
#endif
#if (DIM == 3)
  static const char *varname[]{"UX", "UY", "UZ"};
#endif


    ///------------------------------------------Specific thing to do with some cases ----------------------------------------------///
    if (inputData.bccaseType == CSV_FORCE) {

        // create k-d tree for traction points reading from CSV
         util_funcs::readTXT(inputData.CsvForce.filename,inputData.traction_position_,inputData.traction_vector_);
         util_funcs::removeDuplicates(inputData.traction_position_,inputData.traction_vector_);
        inputData.kdTree_ = new my_kd_tree_t(3 /*dim*/, inputData.traction_position_, 10/* max leaf */);

    }

  ///------------------------------------------Creation/loading of mesh----------------------------------------------///
  DomainExtents domainExtents(inputData.mesh_def.fullDADomain, inputData.mesh_def.physDomain);
  DA *octDA = nullptr;
  DistTREE dTree;
  SubDomain subDomain(domainExtents, resume_from_checkpoint);

  std::vector<GEOMETRY::STL *> stls; // 3D
  std::vector<GEOMETRY::MSH *> mshs; // 2D

  /// IBM setup
  const IBM_METHOD ibmMethod = IBM_METHOD::SBM; // very important to set here!!!!!!!!!! -> inside TalyVec and TalyMat, it will go through different btw NITSCHE and SBM
  IMGA *imga;
  imga = new IMGA(domainExtents, ibmMethod); // [fix]
  std::vector<GEOMETRY::Geometry *> ibm_geoms;
  std::vector<GeomRefinement> ibm_refinements;
  std::vector<std::vector<double>> node_geom;

  for (const auto &geom_def : inputData.ibm_geom_def)
  {
#if (DIM == 2)
    if (geom_def.type == IBMGeomDef::Type::MESHOBJECT_2D)
    {
      std::vector<double> node;
      util_funcs::read_msh(geom_def.mesh_path, node);
      node_geom.emplace_back(node);

      mshs.push_back(new GEOMETRY::MSH(geom_def.mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D));
      std::array<DENDRITE_REAL, DIM> shift;

      shift[0] = geom_def.InitialDisplacement[0];
      shift[1] = geom_def.InitialDisplacement[1];

        if(inputData.bccaseType==CSV_FORCE) {
            for (int dim = 0; dim<DIM;dim++) {
                inputData.shift_(dim) = shift[dim];
            }
        }

      auto geom_retain_side = RetainSide::OUT;
      if (geom_def.outer_boundary)
      {
        geom_retain_side = RetainSide::IN;
      }

      ibm_geoms.emplace_back(new GEOMETRY::Geometry(mshs.back(), Point<DIM>(shift), geom_retain_side));
      ibm_refinements.emplace_back(geom_def.geomRefine);
      GeomRefinement geomRefine;
    }
#endif

#if (DIM == 3)
    if (geom_def.type == IBMGeomDef::Type::MESHOBJECT)
    {
      stls.push_back(new GEOMETRY::STL(geom_def.mesh_path, GEOMETRY::InOutTest::RAY_TRACING));
      std::array<DENDRITE_REAL, DIM> shift;
      shift[0] = geom_def.InitialDisplacement[0];
      shift[1] = geom_def.InitialDisplacement[1];
      shift[2] = geom_def.InitialDisplacement[2];
      if(inputData.bccaseType==CSV_FORCE) {
        for (int dim = 0; dim<DIM;dim++) {
          inputData.shift_(dim) = shift[dim];
        }
      }

      auto geom_retain_side = RetainSide::OUT;
      if (geom_def.outer_boundary)
      {
        geom_retain_side = RetainSide::IN;
      }
      ibm_geoms.emplace_back(new GEOMETRY::Geometry(stls.back(), Point<DIM>(shift), geom_retain_side));
      ibm_refinements.emplace_back(geom_def.geomRefine);
    }
#endif
  }

//  if (inputData.SbmGeo == LEInputData::SBMGeo::RING)
//  {
//// make sure we go here pleaseee
//
////////////////////////////// SBM RING ////////////////////////////
//std::cout<<"SBM RING"<<std::endl;
//
    const double coordsCenter[DIM]{1.25, 1.25};
    subDomain.addObject(VOXEL::Circle(coordsCenter, 1.0 - 2.5 / pow(2, inputData.mesh_def.refine_lvl_base), RetainSide::IN));
    subDomain.addObject(VOXEL::Circle(coordsCenter, 0.25 + 2.5 / pow(2, inputData.mesh_def.refine_lvl_base), RetainSide::OUT));
//  }
//  else
//  {
    for (const auto &c : ibm_geoms)
    {
      subDomain.addObject(c);
    }
//  }

  /// carving subda
  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *physCoords, double physSize)
  {
    return (subDomain.functionToRetain(physCoords, physSize));
  };

  for (int i = 0; i < ibm_geoms.size(); i++)
  {
    imga->addGeometry(ibm_geoms.at(i), ibm_refinements.at(i));
  }

#if (DIM ==3)
    if (inputData.SbmGeo == LEInputData::SBMGeo::cantilever) {

        std::vector<std::pair<double, double>> minmax_cantilever(3, {std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()});

        for (int geoID = 0; geoID < imga->getGeometries().size(); geoID++)
        {
            std::vector<GEOMETRY::Triangles> m_triangles = imga->getGeometries()[geoID]->getSTL()[0].getTriangles();

            for (const auto& triangle : m_triangles) {
                for (int j = 0; j < 3; ++j) {
                    for (int d = 0; d< DIM;d++) {
                        minmax_cantilever[d].first = std::min(minmax_cantilever[d].first, triangle.triangleCoord[j][d] +
                                                                    inputData.ibm_geom_def[geoID].InitialDisplacement(d));
                        minmax_cantilever[d].second = std::max(minmax_cantilever[d].second, triangle.triangleCoord[j][d] +
                                                                                          inputData.ibm_geom_def[geoID].InitialDisplacement(d));

                    }
                }
            }
        }

        inputData.minmax_cantilever = minmax_cantilever;
    }
#endif

    /// Pre Distance Calculation
    PointCloud<double> CenterPts;

#if (DIM == 3)

    if (inputData.DistCalcType == LEInputData::KD_TREE)
    {
        util_funcs::GetTriangleCenter(stls, CenterPts);
    }
#endif


    my_kd_tree_t kd_tree(3 /*dim*/, CenterPts, {10 /* max leaf */});

    octDA = createSubDA(dTree, functionToRetain, levelBase, eleOrder);
  subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
//  int no_refine = util_funcs::performRefinementSubDA(octDA, domainExtents, dTree, inputData, &subDomain);

//  PrintStatus("Number of refinement = ", no_refine);

#pragma mark OptSug

//    //////////////// Preventing losing elements /////////////////
//  if (inputData.ibm_geom_def.size() != 0)
//  {
//
//    int BoundaryMaxRefinelvl = -1;
//
//    if (no_refine > 0)
//    {
//      for (int geomID = 0; geomID < inputData.ibm_geom_def.size(); geomID++)
//      {
//        BoundaryMaxRefinelvl = std::max(BoundaryMaxRefinelvl,
//                                        static_cast<int>(inputData.ibm_geom_def.at(geomID).refine_lvl));
//      }
//    }
//
//      while (true) {
//          DARefine refine(octDA, dTree.getTreePartFiltered(), domainExtents,inputData.mesh_def.refine_lvl_base, false);
//
//          DA *newDA = refine.getRefineSubDA(dTree, 0.03,
//                                            ot::RemeshPartition::SurrogateOutByIn);
//          if (newDA == nullptr) {
//              break;
//          }
//          std::swap(newDA, octDA);
//          delete newDA;
//
//          subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
//      }
//
//      while (true) {
//          DARefine refine(octDA, dTree.getTreePartFiltered(), domainExtents, inputData.mesh_def.refine_lvl_base + 1, true);
//
//          DA *newDA = refine.getRefineSubDA(dTree, 0.03,
//                                            ot::RemeshPartition::SurrogateOutByIn);
//          if (newDA == nullptr) {
//              break;
//          }
//          std::swap(newDA, octDA);
//          delete newDA;
//
//          subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
//      }
//
//      {
//          while (true) {
//              DACoarse refine(octDA, dTree.getTreePartFiltered(), domainExtents, inputData.mesh_def.refine_lvl_base, false);
//              DA *newDA = refine.getRefineSubDA(dTree, 0.03,
//                                                ot::RemeshPartition::SurrogateInByOut);
//              if (newDA == nullptr) {
//                  break;
//              }
//              std::swap(newDA, octDA);
//              delete newDA;
//
//              subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
//          }
//      }  }

  TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());

  /// --------------------------------------------------------------------------------------------------------------////
  const auto &treePartition = dTree.getTreePartFiltered();
  IO::writeBoundaryElements(octDA, treePartition, "boundary", "subDA", domainExtents);

#pragma mark OptSug
//  Vec nodalFalseElement;
//  std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> elementMarkers, elementMarkerVisualize;
//  util_funcs::generateNewMarkers(octDA, dTree, domainExtents, subDomain, elementMarkers, elementMarkerVisualize, nodalFalseElement, &inputData);

    imga->initIMGAComputation(octDA, treePartition);
    Marker elementMarker(octDA, treePartition, domainExtents, imga, MarkerType::GAUSS_POINT);

  SubDomainBoundary boundary(&subDomain, octDA, domainExtents);

  LEBCSetup LEBC(&boundary, &inputData);

  /// Gridfield setup
  TalyMesh<LENodeData> talyMesh(octDA->getElementOrder());

  // domain boundary
  SubDomainBoundary *subDomainBoundary = nullptr;
  subDomainBoundary = &boundary;

  // ndof
  static const DENDRITE_UINT ndof = LENodeData::LE_DOF;


  TimeInfo ti(0.0, inputData.dt, inputData.totalT);

  // for SC

  auto leEq = new TalyEquation<SSLEEquation, LENodeData>(octDA, dTree.getTreePartFiltered(), domainExtents, ndof,
                                                         &ti, true, subDomainBoundary, &inputData, ibmMethod, imga);

  leEq->equation()->setSbmCalc(&kd_tree);

  std::vector<PetscInt> dirichletNodes;
    LinearSolver *leSolver = setLinearSolver(leEq, octDA, ndof, mfree); // for IBM
//  LinearSolver *leSolver = setLinearSolver(leEq, octDA, dTree, ndof, mfree); // for da-opt

  leSolver->setIBMDirichletNodes(dirichletNodes); // it does not set anything (dirichletNodes do not have anything), but it is here to "prevent segmentation fault" !

  if (inputData.bccaseType == BOTTOM_FORCE)
  {
    GetSurfaceGP getSurfaceGp1(octDA, dTree.getTreePartFiltered(),
                               {VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)},
                               domainExtents, subDomainBoundary, inputData, imga);
    getSurfaceGp1.GetRegion(octDA->getCommActive(), inputData.DomainMax);
  } else if (inputData.bccaseType == CSV_FORCE) {
      //
      GetSurfaceGP getSurfaceGp1(octDA, dTree.getTreePartFiltered(),
                                 {VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)},
                                 domainExtents, subDomainBoundary, inputData, imga);
      getSurfaceGp1.GetMinMaxRegion(octDA->getCommActive(), inputData.DomainMin,inputData.DomainMax);

  }

  inputData.solverOptionsLE.apply_to_petsc_options("-le_");
  {
    KSP SSLEEquation_ksp = leSolver->ksp();
    KSPSetOptionsPrefix(SSLEEquation_ksp, "le_");
    KSPSetFromOptions(SSLEEquation_ksp);
  }

  PrintStatus("Setting BC for LE");
  leSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary
                                          {
                                            Boundary b;

                                            LEBC.setBoundaryConditions(b, pos);
                                            return b; });

#pragma mark IBM-3D-moving
    leEq->assignIBMConstructs(imga, elementMarker.getMarkers().data());
    leEq->setVectors({VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)},
                       SYNC_TYPE::ALL);

  /// Calculate Cmatrix here
  inputData.Cmatrix.resize(3 * (DIM - 1));
  util_funcs::CalcCmatrix(&inputData, inputData.Cmatrix);

  // TODO: put back optimal surrogate boundary and merge with thermoelasticity
#pragma mark OptSug
//  leEq->setVectors({VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX),
//                    VecInfo(nodalFalseElement, 1, LENodeData::NODE_ID)},
//                   SYNC_TYPE::ALL);
//
//  leEq->setVectors({VecInfo(nodalFalseElement, 1, LENodeData::NODE_ID)});
//  leEq->assignIBMConstructs(imga, elementMarkers.data(), LENodeData::NODE_ID);
//
//  CheckSurface checkSurface(octDA, dTree.getTreePartFiltered(), {VecInfo(nodalFalseElement, 1, LENodeData::NODE_ID)}, domainExtents, elementMarkers, &subDomain, subDomainBoundary);
//  // This will correct the element markers
//  checkSurface.correctCycles();
//
//  CheckSurface checkSurfaceVis(octDA, dTree.getTreePartFiltered(), {VecInfo(nodalFalseElement, 1, LENodeData::NODE_ID)}, domainExtents, elementMarkerVisualize, &subDomain, subDomainBoundary);
//  checkSurfaceVis.correctCycles();
//  const char *varnameMark[] = {"marker"};
//  std::vector<double> printMarker(elementMarkerVisualize.size());
//  for (int i = 0; i < elementMarkerVisualize.size(); i++)
//  {
//    printMarker[i] = (double)(elementMarkerVisualize[i].to_ulong());
//  }
//  IO::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), printMarker.data(), "Marker", "marker", varnameMark, domainExtents, true);

  /// Solve
  TimerGroup<MPITimer> timers;
  std::vector<std::string> timer_labels = {"Solving"};
  std::map<std::string, int> timer_tags;
  for (int i = 0; i < timer_labels.size(); i++)
  {
    timer_tags.insert(std::pair<std::string, int>(timer_labels[i], i));
    timers.AddTimer(timer_labels[i]);
  }
  PrintStatus("before solving!");
  timers.Start(timer_tags["Solving"]);
  leSolver->solve();
  timers.Stop(timer_tags["Solving"]);
  PrintStatus("solved!");
  timers.PrintTotalTimeSeconds();

  util_funcs::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), leSolver->getCurrentSolution(), "results",
                             "le", varname,
                             subDomain.domainExtents(), false, false, ndof);

  // TODO: fix the ordering in node data
  CalcStress calcStress(octDA, dTree.getTreePartFiltered(), {VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)}, domainExtents,
                        &subDomain, &inputData, ti);
#if (DIM ==3)
  static const char *stress_varname[]{"strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_xz", "strain_yz",
                                      "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_xz", "stress_yz", "vonMises"};
#endif

#if (DIM ==2)
    static const char *stress_varname[]{"strain_xx", "strain_yy", "strain_xy",
                                      "stress_xx", "stress_yy", "stress_xy", "vonMises"};
#endif

  std::vector<double> StressVectorPerElement;
  calcStress.getElementalstress(StressVectorPerElement);
  IO::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), StressVectorPerElement.data(), "Stress",
                     "Stress", stress_varname,
                     domainExtents, true, false, 6 * (DIM - 1) + 1);

  MPI_Barrier(MPI_COMM_WORLD); // here to make sure we at least have Stress output

  /// L2 error
  const auto analytic_sol = [](TALYFEMLIB::ZEROPTV pos, int dof, double time)
  {
    double r = sqrt(pow(pos.x() - 1, 2) + pow(pos.y() - 1, 2));
    return -r * log(r) / 2 / log(2);
  };

  Vec U_le = leSolver->getCurrentSolution();

  static const char *varname2[]{"Displacment"};

#if (DIM == 2)
  IS x_is, y_is;
  Vec U_mag = util_funcs::GetMag(U_le, x_is, y_is);

  util_funcs::save_timestep(octDA, treePartition, U_mag, 1, ti, subDomain, "leMag", varname2);
  VecInfo v(U_mag, 1, 0);
  Analytic LEAnalytic(octDA, treePartition, v, analytic_sol, subDomain.domainExtents());
  LEAnalytic.getL2error();
#endif

#if (DIM == 3)

  IS x_is;
  IS y_is;
  IS z_is;
  Vec U_mag = util_funcs::GetMag3D(U_le, x_is, y_is, z_is);

  util_funcs::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), leSolver->getCurrentSolution(), "results",
                             "U_mag", varname2,
                             subDomain.domainExtents(), false, false, ndof);

#endif

  MPI_Barrier(MPI_COMM_WORLD); // [fix the bug: missing some parts important

#if (DIM == 2)
  Vec U_x, U_y;
  util_funcs::GetVec(U_le, x_is, y_is, U_x, U_y);
  VecInfo vx(U_x, 1, 0);
  VecInfo vy(U_y, 1, 0);
#endif
#if (DIM == 3)
  Vec U_x, U_y, U_z;
  util_funcs::GetVec(U_le, x_is, y_is, z_is, U_x, U_y, U_z);
  VecInfo vx(U_x, 1, 0);
  VecInfo vy(U_y, 1, 0);
  VecInfo vz(U_z, 1, 0);
#endif
  VecInfo vle(U_le, LENodeData::LE_DOF, 0);

//  /// calculate distance square sum
//  if (imga->getGeometries().size() != 0 and inputData.SbmGeo != LEInputData::NONE)
//  {
//    GetSurfaceGP getSurfaceGp(octDA, treePartition, vle, domainExtents, subDomainBoundary, inputData, imga);
//    double EuclideanDistance, MaxDistance, RMSDistance, SpecialDistance;
//    getSurfaceGp.getEuclideanDistance(EuclideanDistance);
//    PrintStatus("EuclideanDistance = ", EuclideanDistance);
//    getSurfaceGp.getMaxDistance(MaxDistance);
//    PrintStatus("MaxDistance = ", MaxDistance);
//    getSurfaceGp.getRMSDistance(RMSDistance);
//    PrintStatus("RMSDistance = ", RMSDistance);
//  }

  /// Domain error
#pragma mark OptSug

#if (DIM ==2)
    double DomainError[2];
    Marker *elementMarkerBaseOnNode = new Marker(octDA, dTree.getTreePartFiltered(), domainExtents, imga, MarkerType::ELEMENT_NODES);
    inputData.CalcUxError = true;
    CalcError calcErrorx(octDA, dTree.getTreePartFiltered(), vx, domainExtents, &subDomain, &inputData,elementMarkerBaseOnNode->getMarkers());
    calcErrorx.getL2error(DomainError);
    TALYFEMLIB::PrintStatus("[Domain Error] L2, Lf (x-dir) = ", DomainError[0], " ", DomainError[1]);
    const auto &elemErrorX = calcErrorx.getElementalError();
    static const char *varnameErrorX[]{"ErrorX"};
    IO::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), elemErrorX.data(), "Error", "ElemErrorX", varnameErrorX,
                       domainExtents, true);


    inputData.CalcUxError = false;
    CalcError calcErrory(octDA, dTree.getTreePartFiltered(), vy, domainExtents, &subDomain, &inputData,elementMarkerBaseOnNode->getMarkers());
    calcErrory.getL2error(DomainError);
    TALYFEMLIB::PrintStatus("[Domain Error] L2, Lf (y-dir) = ", DomainError[0], " ", DomainError[1]);
    const auto &elemErrorY = calcErrory.getElementalError();
    static const char *varnameErrorY[]{"ErrorY"};
    IO::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), elemErrorY.data(), "Error", "ElemErrorY", varnameErrorY,
                       domainExtents, true);

#endif

  delete leEq;
  delete leSolver;

  dendrite_finalize(octDA);
}
