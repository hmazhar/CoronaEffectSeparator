// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Hammad Mazhar, Ida Critelli
// =============================================================================
// Simulation code for Corona Effect Separator
// =============================================================================

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/collision/ChCCollisionModelParallel.h"
#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_opengl/ChOpenGLWindow.h"
#include "input_output.h"

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::utils;
using namespace chrono::opengl;

real timestep = .0001;
real seconds_to_simulate = 10;

int max_iter = 10;
int num_steps = seconds_to_simulate / timestep;
std::string data_folder = "data/ces";

// mm to cm;
real ws = 100;
real kgm3_kgcm3 = 1.0 / 1000000.0;
Vector gravity(0, -9.80665 * ws, 0);

real drumspeed_rpm = 44.8;
real drumspeed_radss = drumspeed_rpm * ((2.0 * CH_C_PI) / 60.0);

utils::Generator* gen;

ChCoordsys<> CS_Spazzola, CS_Splitter2, CS_Splitter1, CS_Drum;
ChSharedBodyPtr Drum;

real ang = 0;

void RunTimeStep(ChSystemParallelDVI* mSys, const int frame) {
  if (frame % 400 == 0) {
    double r_g = 0.01 * ws;
    double dist = 0.99 * r_g;
    gen->createObjectsBox(utils::POISSON_DISK, dist, Vector(-0.1, 0.0796399719369236 + 0.56, -0.0817847646755438) * ws,
                          Vector(0.14, 0, 0.14) * ws, Vector(0, 0, 0));
  }
  ang -= drumspeed_radss * timestep;
  if (ang <= 0) {
    ang = 2 * CH_C_PI;
  }
  Quaternion q1;
  q1.Q_from_AngZ(ang);
  Drum->SetPos(Vector(0.178496411, 0.42749398, -0.081784765) * ws);
  Drum->SetPos_dt(Vector(0, 0, 0));
  Drum->SetRot(q1);

  Drum->SetWvel_loc(Vector(0, 0, -drumspeed_radss));
}

int main(int argc, char* argv[]) {
  GetLog() << "Initializing Simulator \n";

  ChSystemParallelDVI* system_parallel = new ChSystemParallelDVI;
  system_parallel->SetIntegrationType(ChSystem::INT_ANITESCU);
  system_parallel->GetSettings()->solver.solver_mode = SPINNING;
  system_parallel->GetSettings()->solver.max_iteration_normal = (30);
  system_parallel->GetSettings()->solver.max_iteration_sliding = 30;  //(max_iter * 2);
  system_parallel->GetSettings()->solver.max_iteration_spinning = (30);
  system_parallel->GetSettings()->solver.max_iteration_bilateral = 0;
  system_parallel->GetSettings()->solver.tolerance = (0);
  system_parallel->GetSettings()->solver.alpha = (0);
  system_parallel->GetSettings()->solver.contact_recovery_speed = (10);
  system_parallel->ChangeSolverType(APGD);
  system_parallel->GetSettings()->collision.use_aabb_active = true;
  system_parallel->GetSettings()->collision.aabb_min = R3(-1, -1, -1) * ws;
  system_parallel->GetSettings()->collision.aabb_max = R3(1, 1, 1) * ws;
  system_parallel->GetSettings()->min_threads = 2;
  system_parallel->GetSettings()->collision.collision_envelope = (0.001 * .2 * ws);
  system_parallel->GetSettings()->collision.bins_per_axis = I3(25, 50, 25);
  system_parallel->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_GJK;
  system_parallel->Set_G_acc(gravity);
  system_parallel->SetStep(timestep);

  GetLog() << "Building Model \n";
  ChSharedBodyPtr Spazzola = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Splitter2 = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Splitter = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  Drum = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Truss = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Other = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Conveyer = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

  ChSharedPtr<ChMaterialSurface> mat = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat->SetFriction(.5);
  mat->SetRollingFriction(1.0);

  ChSharedPtr<ChMaterialSurface> mat_conveyor = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat_conveyor->SetFriction(0);

  Vector pos = Vector(0);
  Vector dims = Vector(0);
  real mass = 0;
  {
    mass = 2.604131256;
    pos = Vector(-0.033926415, 0.433220509, -0.081784765) * ws;
    utils::InitializeObject(Spazzola, mass, mat, pos, QUNIT, true, true, 8, 8);
    dims = Vector(0.053, 0.15, 0) * ws;
    Spazzola->SetInertia(utils::CalcCylinderGyration(dims.x, dims.y) * mass);
    utils::AddCylinderGeometry(Spazzola.get_ptr(), dims.x, dims.y, Vector(0), Q_from_AngX(CH_C_PI / 2.0));
    utils::FinalizeObject(Spazzola, system_parallel);
    CS_Spazzola = ChCoordsys<>(pos);
  }
  {
    mass = 0.169782297;
    pos = Vector(0.416156012, 0.286549383, -0.081784765) * ws;
    dims = Vector(.002, 0.075, 0.15) * ws;
    Quaternion rot(0.97236992, 0, 0, 0.233445364);
    utils::InitializeObject(Splitter2, mass, mat, pos, rot, true, true, 8, 8);
    Splitter2->SetInertia(utils::CalcBoxGyration(dims) * mass);
    utils::AddBoxGeometry(Splitter2.get_ptr(), dims, Vector(0.002, 0, 0) * ws, QUNIT);
    utils::FinalizeObject(Splitter2, system_parallel);
  }
  {
    mass = 0.136861816;
    pos = Vector(0.324996411, 0.27399398, -0.081784765) * ws;
    dims = Vector(.002, 0.0535, 0.1515) * ws;
    utils::InitializeObject(Splitter, mass, mat, pos, QUNIT, true, true, 8, 8);
    Splitter->SetInertia(utils::CalcBoxGyration(dims) * mass);
    utils::AddBoxGeometry(Splitter.get_ptr(), dims, Vector(0.002, 0, 0) * ws, QUNIT);
    utils::FinalizeObject(Splitter, system_parallel);
  }
  {
    mass = 24.01511276;
    pos = Vector(0.178496411, 0.42749398, -0.081784765) * ws;
    dims = Vector(0.16, 0.15, 0) * ws;
    utils::InitializeObject(Drum, mass, mat, pos, QUNIT, true, false, 8, 8);
    Drum->SetInertia(utils::CalcCylinderGyration(dims.x, dims.y) * mass);
    utils::AddCylinderGeometry(Drum.get_ptr(), dims.x, dims.y, Vector(0), Q_from_AngX(CH_C_PI / 2.0));
    CS_Drum = ChCoordsys<>(pos);
    utils::FinalizeObject(Drum, system_parallel);
  }
  {
    utils::InitializeObject(Truss, 1, mat, Vector(0, 0, 0), QUNIT, true, true, 8, 8);
    utils::AddBoxGeometry(Truss.get_ptr(), Vector(.5, .6, .002) * ws,
                          Vector(0.1784964, 0.427494, 0.15 - 0.081784765) * ws, QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), Vector(.5, .6, .002) * ws,
                          Vector(0.1784964, 0.427494, -0.15 - 0.081784765) * ws, QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), Vector(.002, .6, .2) * ws, Vector(-0.25, 0.427494, -0.081784765) * ws,
                          QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), Vector(.002, .6, .2) * ws, Vector(.67, 0.427494, -0.081784765) * ws, QUNIT);
    utils::FinalizeObject(Truss, system_parallel);
  }
  {
    mass = 0.575578251;
    pos = Vector(0.181496411, 0.42339398, -0.081784765) * ws;
    utils::InitializeObject(Other, mass, mat, pos, Quaternion(1, 0, 0, 0), true, true, 8, 8);
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.22121571, 0.013748141, 0.448934435) * ws,
                          Vector(0.034434623, -.7, 0.030630547) * ws, Quaternion(0.707107, -0.0, -0.707107, -0.0));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.045869262 * ws, 0.16 * ws, Vector(0, -0.517959057661396, 0) * ws,
                               Quaternion(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.008766574 * ws, 0.16 * ws,
                               Vector(-0.381566155791125, -0.515717402690516, 0) * ws,
                               Quaternion(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.012493547 * ws, 0.16 * ws,
                               Vector(0.492612868005681, -0.513858743768626, 0) * ws, Quaternion(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.024619437 * ws, 0.16 * ws,
                               Vector(0.244352690753348, -0.517959057661396, 0) * ws, Quaternion(0.5, -0.5, -0.5, 0.5));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.121733706, 0.008766574, 0.165) * ws,
                          Vector(-0.324444835447566, -0.623217402690516, 0) * ws,
                          Quaternion(0.515154, 0.0, 0.0, 0.857097));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.121733706, 0.009012905, 0.165) * ws,
                          Vector(-0.0943376359362894, -0.62333298893951, 0) * ws,
                          Quaternion(0.515154, 0.0, 0.0, -0.857097));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.102920917, 0.008807300, 0.165) * ws,
                          Vector(0.0461378510072952, -0.614333388045196, 0) * ws,
                          Quaternion(0.739227, 0.0, 0.0, -0.673456));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.110302982, 0.009073557, 0.165) * ws,
                          Vector(0.201202979597005, -0.620296344571127, 0) * ws,
                          Quaternion(0.611299, 0.0, 0.0, -0.791400));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.111244909, 0.009917229, 0.165) * ws,
                          Vector(0.287088436539427, -0.62133744838765, 0) * ws,
                          Quaternion(0.792869, 0.0, 0.0, -0.609392));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.111244909, 0.010000000, 0.165) * ws,
                          Vector(0.463991547662121, -0.621358743768626, 0) * ws,
                          Quaternion(0.609392, 0.0, 0.0, -0.792869));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.115826155, 0.014764599, 0.447) * ws,
                          Vector(0.0545, -0.624455203991318, 0.106974099367162) * ws,
                          Quaternion(0.585725, -0.396140, 0.585725, -0.396140));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.115826155, 0.014764599, 0.447) * ws,
                          Vector(0.0545, -0.624455203991318, -0.106974099367162) * ws,
                          Quaternion(0.585725, 0.396140, -0.585725, -0.396140));
    utils::AddBoxGeometry(Other.get_ptr(), Vector(0.111244909, 0.018292686, 0.16500005) * ws,
                          Vector(0.455978023535446, -0.619225184316436, 5.E-08) * ws,
                          Quaternion(0.0, 0.792869, 0.609392, 0.0));
    utils::FinalizeObject(Other, system_parallel);
  }
  {
    mass = 3.832724939;
    pos = Vector(-0.0627155647639022, 0.0534995314229745, -0.0817847646755438) * ws;
    ChMatrix33<> mr;

    mr(0, 0) = -0.995090127679177;
    mr(1, 0) = 0.0989729144535987;
    mr(2, 0) = 0;
    mr(0, 1) = 0.0989729144535989;
    mr(1, 1) = 0.995090127679177;
    mr(2, 1) = 0;
    mr(0, 2) = 0;
    mr(1, 2) = 0;
    mr(2, 2) = -1;
    utils::InitializeObject(Conveyer, mass, mat_conveyor, pos, mr.Get_A_quaternion(), true, true, 8, 8);
    utils::AddBoxGeometry(Conveyer.get_ptr(), Vector(0.193724352473189, 0.005, 0.2) * ws,
                          Vector(0.0785377577161679, 0.554163002941846, 0) * ws, QUNIT);
    system_parallel->AddBody(Conveyer);
  }
  ChSharedPtr<ChMaterialSurface> mat_p = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat_p->SetFriction(.2);  // This is so that it is averaged properly
  mat_p->SetRollingFriction(1.0);
  gen = new utils::Generator(system_parallel);
  int Id_g = 100;

  utils::MixtureIngredientPtr& m1 = gen->AddMixtureIngredient(utils::SPHERE, 0.333);
  m1->setDefaultSize(0.0005 * ws);
  m1->setDefaultDensity(946 * kgm3_kgcm3);
  m1->setDistributionSize(0.0005 * ws, 1, 0.0005 * ws, 0.001 * ws);
  m1->setDefaultMaterialDVI(mat_p);

  utils::MixtureIngredientPtr& m2 = gen->AddMixtureIngredient(utils::BOX, 0.333);
  m2->setDefaultSize(0.0005);
  m2->setDefaultDensity(946 * kgm3_kgcm3);
  m2->setDistributionSize(0.0005 * ws, 1, 0.0005 * ws, 0.001 * ws);
  m2->setDefaultMaterialDVI(mat_p);

  utils::MixtureIngredientPtr& m3 = gen->AddMixtureIngredient(utils::CYLINDER, 0.333);
  m3->setDefaultSize(0.0005 * ws);
  m3->setDefaultDensity(946 * kgm3_kgcm3);
  m3->setDistributionSize(0.0005 * ws, 1, 0.0005 * ws, 0.001 * ws);
  m3->setDefaultMaterialDVI(mat_p);

  //  utils::MixtureIngredientPtr& m4 = gen->AddMixtureIngredient(utils::ELLIPSOID, 0.25);
  //  m4->setDefaultSize(0.0005 * ws);
  //  m4->setDefaultDensity(946 * kgm3_kgcm3);
  //  m4->setDistributionSize(0.0005 * ws, 1, 0.0005 * ws, 0.001 * ws);
  //  m4->setDefaultMaterialDVI(mat_p);
  //  gen->setBodyIdentifier(Id_g);

  //  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  //  gl_window.Initialize(1280, 720, "CES", system_parallel);
  //  gl_window.SetCamera(Vector(0, 0, -2) * ws, Vector(0, 0, 0), Vector(0, 1, 0), .05 * ws, 0.01 * ws, 100 * ws);
  //  gl_window.viewer->cloud.SetPointSize(.001 * ws);
  //  gl_window.viewer->contact_renderer.SetPointSize(.001 * ws);
  //  // gl_window.Pause();
  //  for (int i = 0; i < num_steps; i++) {
  //    if (gl_window.Active()) {
  //      if (gl_window.DoStepDynamics(timestep)) {
  //        RunTimeStep(system_parallel, i);
  //        TimingOutput(system_parallel);
  //        if (i == 0) {
  //          system_parallel->GetSettings()->collision.aabb_min =
  //              system_parallel->data_manager->measures.collision.min_bounding_point;
  //          system_parallel->GetSettings()->collision.aabb_max =
  //              system_parallel->data_manager->measures.collision.max_bounding_point;
  //        }
  //      }
  //      gl_window.Render();
  //    }
  //  }
  //  exit(0);

  int file = 0;
  int save_every = 1.0 / timestep / 600.0;
  for (int i = 0; i < num_steps; i++) {
    system_parallel->DoStepDynamics(timestep);
    RunTimeStep(system_parallel, i);
    TimingOutput(system_parallel);

    if (i == 0) {
      system_parallel->GetSettings()->collision.aabb_min =
          system_parallel->data_manager->measures.collision.min_bounding_point;
      system_parallel->GetSettings()->collision.aabb_max =
          system_parallel->data_manager->measures.collision.max_bounding_point;
    }

    if (i % save_every == 0) {
      std::stringstream ss;
      std::cout << "Frame: " << file << std::endl;
      ss << data_folder << "/" << file << ".txt";
      DumpAllObjectsWithGeometryPovray(system_parallel, ss.str(), true);
      file++;
    }
  }
  return 0;
}
