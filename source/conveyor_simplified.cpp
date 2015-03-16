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
//#include "ElectricParticleProperty.h"
//#include "ElectricForcesCES.h"

// Use the namespace of Chrono
using namespace chrono;
using namespace chrono::collision;
using namespace chrono::utils;
using namespace chrono::opengl;

real timestep = .0001;
real seconds_to_simulate = 10;

int max_iter = 10;
int num_steps = seconds_to_simulate / timestep;
std::string data_folder = "data/ces";
// real average_shape_size = .01;

ChVector<> gravity(0, -9.80665, 0);

real drumspeed_rpm = 44.8;
real drumspeed_radss = drumspeed_rpm * ((2.0 * CH_C_PI) / 60.0);

utils::Generator* gen;
// ElectricForcesCES ces_forces;

ChCoordsys<> CS_Spazzola, CS_Splitter2, CS_Splitter1, CS_Drum;

// Optional: define a callback to be exectuted at each creation of a metal particle:
class MyCreator_metal : public utils::CallbackGenerator {
  // Here do custom stuff on the just-created particle:
 public:
  virtual void PostCreation(ChSharedPtr<ChBody> mbody) {

    //    ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty);
    //    electric_asset->fraction = ElectricParticleProperty::e_fraction_sphere;
    //    electric_asset->material_type = ElectricParticleProperty::e_mat_metal;
    //    electric_asset->conductivity = 58000000;
    //    electric_asset->birthdate = this->systemreference->GetChTime();
    //    ChVector<> Cradii;    // use equivalent-inertia ellipsoid to get characteristic size:
    //    ChVector<> Ine = mbody->GetInertiaXX();
    //    Cradii.x = sqrt((5. / (2. * mbody->GetMass())) * (Ine.y + Ine.z - Ine.x));
    //    Cradii.y = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.z - Ine.y));
    //    Cradii.z = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.y - Ine.z));
    //    electric_asset->Cdim = Cradii * 2.;
    //    mbody->AddAsset(electric_asset);
  }
  // here put custom data that might be needed by the callback:
  ChSystem* systemreference;
};

// Optional: define a callback to be exectuted at each creation of a metal particle:
class MyCreator_plastic : public utils::CallbackGenerator {
  // Here do custom stuff on the just-created particle:
 public:
  virtual void PostCreation(ChSharedPtr<ChBody> mbody) {

    //    ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty);
    //    electric_asset->fraction = ElectricParticleProperty::e_fraction_sphere;
    //    electric_asset->material_type = ElectricParticleProperty::e_mat_plastic;
    //    electric_asset->conductivity = 0;
    //    electric_asset->birthdate = this->systemreference->GetChTime();
    //    ChVector<> Cradii;    // use equivalent-inertia ellipsoid to get characteristic size:
    //    ChVector<> Ine = mbody->GetInertiaXX();
    //    Cradii.x = sqrt((5. / (2. * mbody->GetMass())) * (Ine.y + Ine.z - Ine.x));
    //    Cradii.y = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.z - Ine.y));
    //    Cradii.z = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.y - Ine.z));
    //    electric_asset->Cdim = Cradii * 2.;
    //    mbody->AddAsset(electric_asset);
  }
  // here put custom data that might be needed by the callback:
  ChSystem* systemreference;
};

void RunTimeStep(ChSystemParallelDVI* mSys, const int frame) {
  //.5mm - 1mm (dia)
  // 1ms timestep
  if (frame % 250 == 0) {
    double r_g = 0.01;
    double dist =  0.99 * r_g;
    gen->createObjectsBox(utils::POISSON_DISK, dist, ChVector<>(-0.1, 0.0796399719369236 + 0.56, -0.0817847646755438), ChVector<>(0.14, 0, 0.14), ChVector<>(0, 0, 0));
  }
  // Larger radius will not create forces

  //  ces_forces.apply_forces(mSys,               // contains all bodies
  //                          CS_Drum,            // pos and rotation of axis of drum (not rotating reference!)
  //                          drumspeed_radss,    // speed of drum
  //                          frame);
}

int main(int argc, char* argv[]) {
  GetLog() << "Initializing Simulator \n";

  ChSystemParallelDVI* system_parallel = new ChSystemParallelDVI;
  system_parallel->SetIntegrationType(ChSystem::INT_ANITESCU);
  system_parallel->GetSettings()->solver.max_iteration_normal = (30);
  system_parallel->GetSettings()->solver.max_iteration_sliding = 30;    //(max_iter * 2);
  system_parallel->GetSettings()->solver.max_iteration_spinning = (0);
  system_parallel->GetSettings()->solver.max_iteration_bilateral = 30;
  system_parallel->GetSettings()->solver.tolerance = (.00);
  system_parallel->GetSettings()->solver.alpha = (0);
  system_parallel->GetSettings()->solver.contact_recovery_speed = (.2);
  system_parallel->ChangeSolverType(APGD);
  system_parallel->GetSettings()->min_threads = 2;
  system_parallel->GetSettings()->collision.collision_envelope = (0.001 * .05);
  system_parallel->GetSettings()->collision.bins_per_axis = I3(25, 50, 25);
  system_parallel->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
  system_parallel->Set_G_acc(gravity);
  system_parallel->SetStep(timestep);

  GetLog() << "Building Model \n";
  ChSharedBodyPtr Spazzola = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Splitter2 = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Splitter = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Drum = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Truss = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Other = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
  ChSharedBodyPtr Conveyer = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

  ChSharedPtr<ChMaterialSurface> mat = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat->SetFriction(0.5);

  {
    real mass = 2.604131256;
    ChVector<> pos(-0.033926415, 0.433220509, -0.081784765);
    utils::InitializeObject(Spazzola, mass, mat, pos, QUNIT, true, false, 8, 8);
    Spazzola->SetInertiaXX(ChVector<>(0.021381139, 0.021381139, 0.00358415));
    utils::AddCylinderGeometry(Spazzola.get_ptr(), 0.053, 0.15, ChVector<>(0), Q_from_AngX(CH_C_PI / 2.0));
    utils::FinalizeObject(Spazzola, system_parallel);
    CS_Spazzola = ChCoordsys<>(pos);
  }

  {
    real mass = 0.169782297;
    ChVector<> pos(0.416156012, 0.286549383, -0.081784765);
    ChQuaternion<> rot(0.97236992, 0, 0, 0.233445364);
    utils::InitializeObject(Splitter2, mass, mat, pos, rot, true, false, 8, 8);
    Splitter2->SetInertiaXX(ChVector<>(0.001599931, 0.001358747, 0.000411944));
    utils::AddBoxGeometry(Splitter2.get_ptr(), ChVector<>(.002, 0.075, 0.15), ChVector<>(0.002, 0, 0), QUNIT);
    Splitter2->SetBodyFixed(true);
    utils::FinalizeObject(Splitter2, system_parallel);
  }

  {
    real mass = 0.136861816;
    ChVector<> pos(0.324996411, 0.27399398, -0.081784765);
    utils::InitializeObject(Splitter, mass, mat, pos, QUNIT, true, false, 8, 8);
    Splitter->SetInertiaXX(ChVector<>(0.001202514, 0.001027289, 0.000176832));
    utils::AddBoxGeometry(Splitter.get_ptr(), ChVector<>(.002, 0.0535, 0.1515), ChVector<>(0.002, 0, 0), QUNIT);
    Splitter->SetBodyFixed(true);
    utils::FinalizeObject(Splitter, system_parallel);
  }

  {
    real mass = 24.01511276;
    ChVector<> pos(0.178496411, 0.42749398, -0.081784765);
    utils::InitializeObject(Drum, mass, mat, pos, QUNIT, true, false, 8, 8);
    Drum->SetInertiaXX(ChVector<>(0.332285546, 0.332285546, 0.308075523));
    utils::AddCylinderGeometry(Drum.get_ptr(), 0.16, 0.15, ChVector<>(0), Q_from_AngX(CH_C_PI / 2.0));
    CS_Drum = ChCoordsys<>(pos);
    utils::FinalizeObject(Drum, system_parallel);
  }
  {
    utils::InitializeObject(Truss, 1, mat, ChVector<>(0, 0, 0), QUNIT, true, true, 8, 8);
    utils::AddBoxGeometry(Truss.get_ptr(), ChVector<>(.5, .6, .002), ChVector<>(0.178496411, 0.42749398, 0.15 - 0.081784765), QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), ChVector<>(.5, .6, .002), ChVector<>(0.178496411, 0.42749398, -0.15 - 0.081784765), QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), ChVector<>(.002, .6, .2), ChVector<>(-0.25, 0.42749398, -0.0817847646755438), QUNIT);
    utils::AddBoxGeometry(Truss.get_ptr(), ChVector<>(.002, .6, .2), ChVector<>(.67, 0.42749398, -0.0817847646755438), QUNIT);
    utils::FinalizeObject(Truss, system_parallel);
  }
  {
    real mass = 0.575578251;
    ChVector<> pos(0.181496411, 0.42339398, -0.081784765);
    utils::InitializeObject(Other, mass, mat, pos, ChQuaternion<>(1, 0, 0, 0), true, false, 8, 8);
    Other->SetInertiaXX(ChVector<>(0.229283207, 0.138637511, 0.320113118));

//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.22121571, 0.013748141, 0.448934435), ChVector<>(0.034434623, -1.123553005, 0.030630547), ChQuaternion<>(0.707107, -0.0, -0.707107, -0.0));
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.010173215, 0.250351436, 0.448934435), ChVector<>(0.034434661, -1.22725133, 0.241673042), ChQuaternion<>(0.707107, -0.0, -0.707107, -0.0));
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.011042497, 0.250103825, 0.448934435), ChVector<>(0.034434585, -1.227003564, -0.179542666), ChQuaternion<>(0.0, -0.707107, 0.0, 0.707107));
//
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.013251174, 0.254046191, 0.2035), ChVector<>(-0.401248638, -1.230946291, 0.031500059), ChQuaternion<>(0.0, 0.0, 0.0, 1.0));
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.029255479, 0.249805368, 0.2035), ChVector<>(-0.061741847, -1.226705282, 0.031499996), ChQuaternion<>(-0.0, 0.0, 0.0, 1.0));
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.019763692, 0.250242141, 0.2035), ChVector<>(0.209277317, -1.2271419, 0.031499947), ChQuaternion<>(-0.0, 0.0, 0.0, 1.0));
//    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.015414096, 0.250130195, 0.2035), ChVector<>(0.467954961, -1.227029814, 0.031499901), ChQuaternion<>(1.0, -0.0, 0.0, -0.0));

    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.22121571, 0.013748141, 0.448934435), ChVector<>(0.034434623, -.7, 0.030630547), ChQuaternion<>(0.707107, -0.0, -0.707107, -0.0));

    utils::AddCylinderGeometry(Other.get_ptr(), 0.045869262, 0.16, ChVector<>(0, -0.517959057661396, 0), ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.008766574, 0.16, ChVector<>(-0.381566155791125, -0.515717402690516, 0), ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.012493547, 0.16, ChVector<>(0.492612868005681, -0.513858743768626, 0), ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(), 0.024619437, 0.16, ChVector<>(0.244352690753348, -0.517959057661396, 0), ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.121733706251767, 0.00876657438055399, 0.165), ChVector<>(-0.324444835447566, -0.623217402690516, 0),
                          ChQuaternion<>(0.515154, 0.0, 0.0, 0.857097));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.121733706251767, 0.00901290521459452, 0.165), ChVector<>(-0.0943376359362894, -0.62333298893951, 0),
                          ChQuaternion<>(0.515154, 0.0, 0.0, -0.857097));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.102920917614647, 0.00880730074763987, 0.165), ChVector<>(0.0461378510072952, -0.614333388045196, 0),
                          ChQuaternion<>(0.739227, 0.0, 0.0, -0.673456));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.110302982718799, 0.009073557274052, 0.165), ChVector<>(0.201202979597005, -0.620296344571127, 0),
                          ChQuaternion<>(0.611299, 0.0, 0.0, -0.791400));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.111244909897975, 0.00991722943912151, 0.165), ChVector<>(0.287088436539427, -0.62133744838765, 0),
                          ChQuaternion<>(0.792869, 0.0, 0.0, -0.609392));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.111244909897976, 0.01, 0.165), ChVector<>(0.463991547662121, -0.621358743768626, 0), ChQuaternion<>(0.609392, 0.0, 0.0, -0.792869));
    // utils::AddCylinderGeometry(Other.get_ptr(),0.0147645993377839,0.447,ChVector<>(0.0545,-0.516955203991318,0.150095419710721),ChQuaternion<>(0.5,0.0,0.0,0.5));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.115826155371625, 0.0147645993377839, 0.447), ChVector<>(0.0545, -0.624455203991318, 0.106974099367162),
                          ChQuaternion<>(0.585725, -0.396140, 0.585725, -0.396140));
    // utils::AddCylinderGeometry(Other.get_ptr(),0.0147645993377839,0.447,ChVector<>(0.0545,-0.516955203991318,-0.150095419710721),ChQuaternion<>(0.5,0.0,0.0,0.5));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.115826155371625, 0.0147645993377839, 0.447), ChVector<>(0.0545, -0.624455203991318, -0.106974099367162),
                          ChQuaternion<>(0.585725, 0.396140, -0.585725, -0.396140));
    utils::AddBoxGeometry(Other.get_ptr(), ChVector<>(0.111244909897976, 0.018292686227323, 0.16500005), ChVector<>(0.455978023535446, -0.619225184316436, 5.00000000014378E-08),
                          ChQuaternion<>(0.0, 0.792869, 0.609392, 0.0));
    Other->SetBodyFixed(true);
    utils::FinalizeObject(Other, system_parallel);
  }

  {
    real mass = 3.832724939;
    ChVector<> pos(-0.0627155647639022, 0.0534995314229745, -0.0817847646755438);
    ChQuaternion<> rot(1, 0, 0, 0);
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
    rot = mr.Get_A_quaternion();
    utils::InitializeObject(Conveyer, mass, mat, pos, rot, true, false, 8, 8);
    // Conveyer->SetInertiaXX(ChVector<>(0.029314394817887, 0.0796399719369236, 0.0514634926681141));
    utils::AddBoxGeometry(Conveyer.get_ptr(), ChVector<>(0.193724352473189, 0.005, 0.2), ChVector<>(0.0785377577161679, 0.554163002941846, 0), QUNIT);
    Conveyer->SetBodyFixed(true);
    system_parallel->AddBody(Conveyer);
  }
  //
  ChSharedPtr<ChLinkEngine> engine_drum = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
  ChSharedPtr<ChLinkEngine> engine_spazzola = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);

  engine_drum->Initialize(Drum, Truss, CS_Drum);
  engine_drum->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
  if (ChSharedPtr<ChFunction_Const> mfun = (engine_drum->Get_spe_funct().DynamicCastTo<ChFunction_Const>())) {
    mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]
  }
  system_parallel->Add(engine_drum);

  engine_spazzola->Initialize(Spazzola, Truss, CS_Spazzola);
  engine_spazzola->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
  if (ChSharedPtr<ChFunction_Const> mfun = (engine_spazzola->Get_spe_funct().DynamicCastTo<ChFunction_Const>())) {
    mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]
  }
  system_parallel->Add(engine_spazzola);

  ChSharedPtr<ChMaterialSurface> mat_p = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat_p->SetFriction(.2);    // This is so that it is averaged properly

  gen = new utils::Generator(system_parallel);
  int Id_g = 100;

  MyCreator_metal* callback_metal = new MyCreator_metal;
  callback_metal->systemreference = system_parallel;

  MyCreator_plastic* callback_plastic = new MyCreator_plastic;
  callback_plastic->systemreference = system_parallel;

  utils::MixtureIngredientPtr& m1 = gen->AddMixtureIngredient(utils::SPHERE, 0.25);
  m1->setDefaultSize(0.0005);
  m1->setDefaultDensity(8400);
  m1->setDistributionSize(0.0005, 1, 0.0005, 0.001);
  m1->setDefaultMaterialDVI(mat_p);
  m1->SetCallbackPostCreation(callback_metal);

  utils::MixtureIngredientPtr& m2 = gen->AddMixtureIngredient(utils::BOX, 0.25);
  m2->setDefaultSize(0.0005);
  m2->setDefaultDensity(946);
  m2->setDistributionSize(0.0005, 1, 0.0005, 0.001);
  m2->setDefaultMaterialDVI(mat_p);
  m2->SetCallbackPostCreation(callback_plastic);

  utils::MixtureIngredientPtr& m3 = gen->AddMixtureIngredient(utils::CYLINDER, 0.25);
  m3->setDefaultSize(0.0005);
  m3->setDefaultDensity(946);
  m3->setDistributionSize(0.0005, 1, 0.0005, 0.001);
  m3->setDefaultMaterialDVI(mat_p);
  m3->SetCallbackPostCreation(callback_plastic);

  utils::MixtureIngredientPtr& m4 = gen->AddMixtureIngredient(utils::ELLIPSOID, 0.25);
  m4->setDefaultSize(0.0005);
  m4->setDefaultDensity(946);
  m4->setDistributionSize(0.0005, 1, 0.0005, 0.001);
  m4->setDefaultMaterialDVI(mat_p);
  m4->SetCallbackPostCreation(callback_plastic);
  gen->setBodyIdentifier(Id_g);

//  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
//  gl_window.Initialize(1280, 720, "CES", system_parallel);
//  gl_window.SetCamera(ChVector<>(0, 0, -2), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), .05f, 0.01, 100);
//  gl_window.viewer->cloud.SetPointSize(.001);
//  gl_window.viewer->contact_renderer.SetPointSize(.001);
//  // gl_window.Pause();
//  for (int i = 0; i < num_steps; i++) {
//
//    if (gl_window.Active()) {
//      if (gl_window.DoStepDynamics(timestep)) {
//        RunTimeStep(system_parallel, i);
//        TimingOutput(system_parallel);
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
