///////////////////////////////////////////////////
//
//   Corona Electrostatic Separator
//
//   This program is based on the following
//   libraries:
//   - ChronoEngine
//   - ChronoEngine Parallel
//
// ------------------------------------------------
///////////////////////////////////////////////////
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/collision/ChCCollisionModelParallel.h"
#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGeometry.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_opengl/ChOpenGLWindow.h"

#include "input_output.h"
#include "ElectricParticleProperty.h"
#include "ElectricForcesCES.h"

// Use the namespace of Chrono
using namespace chrono;
using namespace chrono::collision;
using namespace chrono::utils;
using namespace chrono::opengl;

real timestep = .0005;
real seconds_to_simulate = 10;

int max_iter = 10;
int num_steps = seconds_to_simulate / timestep;

real average_shape_size = .003;

ChVector<> gravity(0, -9.80665, 0);

real drumspeed_rpm = 44.8;
real drumspeed_radss = drumspeed_rpm * ((2.0 * CH_C_PI) / 60.0);

utils::Generator* gen;

// Optional: define a callback to be exectuted at each creation of a metal particle:
class MyCreator_metal : public utils::CallbackGenerator {
  // Here do custom stuff on the just-created particle:
 public:
  virtual void PostCreation(ChSharedPtr<ChBody> mbody) {

    ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty);
    electric_asset->fraction = ElectricParticleProperty::e_fraction_sphere;
    electric_asset->material_type = ElectricParticleProperty::e_mat_metal;
    electric_asset->conductivity = 58000000;
    electric_asset->birthdate = this->systemreference->GetChTime();
    ChVector<> Cradii;    // use equivalent-inertia ellipsoid to get characteristic size:
    ChVector<> Ine = mbody->GetInertiaXX();
    Cradii.x = sqrt((5. / (2. * mbody->GetMass())) * (Ine.y + Ine.z - Ine.x));
    Cradii.y = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.z - Ine.y));
    Cradii.z = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.y - Ine.z));
    electric_asset->Cdim = Cradii * 2.;
    mbody->AddAsset(electric_asset);
  }
  // here put custom data that might be needed by the callback:
  ChSystem* systemreference;
};

// Optional: define a callback to be exectuted at each creation of a metal particle:
class MyCreator_plastic : public utils::CallbackGenerator {
  // Here do custom stuff on the just-created particle:
 public:
  virtual void PostCreation(ChSharedPtr<ChBody> mbody) {

    ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty);
    electric_asset->fraction = ElectricParticleProperty::e_fraction_sphere;
    electric_asset->material_type = ElectricParticleProperty::e_mat_plastic;
    electric_asset->conductivity = 0;
    electric_asset->birthdate = this->systemreference->GetChTime();
    ChVector<> Cradii;    // use equivalent-inertia ellipsoid to get characteristic size:
    ChVector<> Ine = mbody->GetInertiaXX();
    Cradii.x = sqrt((5. / (2. * mbody->GetMass())) * (Ine.y + Ine.z - Ine.x));
    Cradii.y = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.z - Ine.y));
    Cradii.z = sqrt((5. / (2. * mbody->GetMass())) * (Ine.x + Ine.y - Ine.z));
    electric_asset->Cdim = Cradii * 2.;
    mbody->AddAsset(electric_asset);
  }
  // here put custom data that might be needed by the callback:
  ChSystem* systemreference;
};

void RunTimeStep(ChSystemParallelDVI* mSys, const int frame) {


	if(frame%100==0){
		double r_g = 0.01;
		double dist = 2 * 0.99 * r_g;
		gen->createObjectsBox(utils::REGULAR_GRID,
				dist,
				ChVector<>(0.029314394817887 + 0.15, 0.0796399719369236 + 0.6, -0.0817847646755438),
				ChVector<>(0.14, 0,0.14),
				ChVector<>(0, 0, 0));
	}

}

int main(int argc, char* argv[]) {
  GetLog() << "Initializing Simulator \n";

  ChSystemParallelDVI* system_parallel = new ChSystemParallelDVI;
  system_parallel->SetIntegrationType(ChSystem::INT_ANITESCU);
  system_parallel->GetSettings()->solver.max_iteration_normal = (max_iter);
  system_parallel->GetSettings()->solver.max_iteration_sliding = (max_iter * 2);
  system_parallel->GetSettings()->solver.max_iteration_spinning = (0);
  system_parallel->GetSettings()->solver.tolerance = (.01);
  system_parallel->GetSettings()->solver.alpha = (0);
  system_parallel->GetSettings()->solver.contact_recovery_speed = (.5);
  system_parallel->ChangeSolverType(APGD);
  system_parallel->GetSettings()->min_threads = 8;
  system_parallel->GetSettings()->collision.collision_envelope = (average_shape_size * .01);
  system_parallel->GetSettings()->collision.max_body_per_bin = 100;
  system_parallel->GetSettings()->collision.min_body_per_bin = 50;
  system_parallel->GetSettings()->collision.bins_per_axis = I3(100, 150, 100);

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

  ChCoordsys<> CS_Spazzola, CS_Splitter2, CS_Splitter1, CS_Drum;
  {
    real mass = 2.60413125593896;
    ChVector<> pos(-0.0339264149957167, 0.433220509447178, -0.0817847646755438);
    utils::InitializeObject(Spazzola, mass, mat, pos, QUNIT, true, false, 8, 8);
    Spazzola->SetInertiaXX(ChVector<>(0.0213811392305318, 0.0213811392305318, 0.00358414990389192));

    //    {
    //      geometry::ChTriangleMeshConnected trimesh;
    //      trimesh.LoadWavefrontMesh("body_1_1.obj", false, false);
    //      ChSharedPtr<ChTriangleMeshShape> trimesh_shape(new ChTriangleMeshShape);
    //      trimesh_shape->SetMesh(trimesh);
    //      trimesh_shape->SetName("body_1_1.obj");
    //      trimesh_shape->Pos = ChVector<>(0, 0, 0);
    //      trimesh_shape->Rot = ChQuaternion<>(1, 0, 0, 0);
    //      Spazzola->GetAssets().push_back(trimesh_shape);
    //    }

    utils::AddCylinderGeometry(
        Spazzola.get_ptr(), 0.053, 0.15, ChVector<>(0), Q_from_AngX(CH_C_PI / 2.0));

    // utils::AddTriangleMeshGeometry(Spazzola.get_ptr(), "body_1_1.obj", "body_1_1.obj");
    system_parallel->AddBody(Spazzola);

    CS_Spazzola = ChCoordsys<>(pos);
  }

  {
    real mass = 0.169782296961291;
    ChVector<> pos(0.416156012458144, 0.286549383384347, -0.0817847646755439);
    ChQuaternion<> rot(0.972369920397677, 0, 0, 0.233445363855904);
    utils::InitializeObject(Splitter2, mass, mat, pos, rot, true, false, 8, 8);
    Splitter2->SetInertiaXX(
        ChVector<>(0.00159993103963307, 0.00135874722183981, 0.000411943807053511));
    // utils::AddTriangleMeshGeometry(Splitter2.get_ptr(), "body_2_1.obj", "body_2_1.obj");
    utils::AddBoxGeometry(
        Splitter2.get_ptr(), ChVector<>(.002, 0.075, 0.15), ChVector<>(0.002, 0, 0), QUNIT);
    Splitter2->SetBodyFixed(true);
    system_parallel->AddBody(Splitter2);
  }

  {
    real mass = 0.136861816137753;
    ChVector<> pos(0.324996410674918, 0.273993980472249, -0.0817847646755438);
    ChQuaternion<> rot(1, 0, 0, 0);
    utils::InitializeObject(Splitter, mass, mat, pos, rot, true, false, 8, 8);
    Splitter->SetInertiaXX(
        ChVector<>(0.00120251356875357, 0.00102728888618058, 0.000176832203169258));
    // utils::AddTriangleMeshGeometry(Splitter.get_ptr(), "body_3_1.obj", "body_3_1.obj");
    utils::AddBoxGeometry(
        Splitter.get_ptr(), ChVector<>(.002, 0.0535, 0.1515), ChVector<>(0.002, 0, 0), QUNIT);
    Splitter->SetBodyFixed(true);
    system_parallel->AddBody(Splitter);
  }

  {
    real mass = 24.0151127638079;
    ChVector<> pos(0.178496410674923, 0.427493980472233, -0.0817847646755438);
    ChQuaternion<> rot(1, 0, 0, 0);
    utils::InitializeObject(Drum, mass, mat, pos, rot, true, false, 8, 8);
    Drum->SetInertiaXX(ChVector<>(0.332285545611457, 0.332285545611457, 0.308075522985403));
    //    {
    //      geometry::ChTriangleMeshConnected trimesh;
    //      trimesh.LoadWavefrontMesh("body_4_1.obj", false, false);
    //      ChSharedPtr<ChTriangleMeshShape> trimesh_shape(new ChTriangleMeshShape);
    //      trimesh_shape->SetMesh(trimesh);
    //      trimesh_shape->SetName("body_4_1.obj");
    //      trimesh_shape->Pos = ChVector<>(0, 0, 0);
    //      trimesh_shape->Rot = ChQuaternion<>(1, 0, 0, 0);
    //      Drum->GetAssets().push_back(trimesh_shape);
    //    }
    utils::AddCylinderGeometry(
        Drum.get_ptr(), 0.16, 0.15, ChVector<>(0), Q_from_AngX(CH_C_PI / 2.0));
    // utils::AddTriangleMeshGeometry(Drum.get_ptr(), "body_4_1.obj", "body_4_1.obj");
    system_parallel->AddBody(Drum);
    CS_Drum = ChCoordsys<>(pos);
  }
  {
    utils::InitializeObject(Truss, 1, mat, ChVector<>(0, 0, 0), QUNIT, true, true, 8, 8);

    utils::AddBoxGeometry(
        Truss.get_ptr(),
        ChVector<>(.5, .6, .002),
        ChVector<>(0.178496410674923, 0.427493980472233, 0.15 - 0.0817847646755438),
        ChQuaternion<>(1, 0, 0, 0));
    utils::AddBoxGeometry(
        Truss.get_ptr(),
        ChVector<>(.5, .6, .002),
        ChVector<>(0.178496410674923, 0.427493980472233, -0.15 - 0.0817847646755438),
        ChQuaternion<>(1, 0, 0, 0));

    system_parallel->AddBody(Truss);
  }
  {
    real mass = 0.575578251261412;
    ChVector<> pos(0.181496410674923, 0.423393980472242, -0.0817847646755503);
    ChQuaternion<> rot(1, 0, 0, 0);
    utils::InitializeObject(Other, mass, mat, pos, rot, true, false, 8, 8);
    Other->SetInertiaXX(ChVector<>(0.229283207371073, 0.138637511031281, 0.320113118266443));

    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.221215709664422, 0.0137481406638926, 0.448934434860076),
                          ChVector<>(0.0344346226988081, -1.12355300540107, 0.0306305470658227),
                          ChQuaternion<>(0.707107, -0.0, -0.707107, -0.0));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0101732150271489, 0.250351436400839, 0.448934434860075),
                          ChVector<>(0.0344346606165717, -1.22725133002599, 0.241673041703092),
                          ChQuaternion<>(0.707107, -0.0, -0.707107, -0.0));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0110424966507172, 0.250103825371303, 0.448934434860075),
                          ChVector<>(0.0344345849372048, -1.22700356432596, -0.179542665947878),
                          ChQuaternion<>(0.0, -0.707107, 0.0, 0.707107));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0132511736464386, 0.254046190843515, 0.2035),
                          ChVector<>(-0.401248638358623, -1.23094629097478, 0.0315000593787049),
                          ChQuaternion<>(0.0, 0.0, 0.0, 1.0));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0292554789734579, 0.249805368165055, 0.2035),
                          ChVector<>(-0.0617418466516264, -1.2267052816943, 0.0314999958019597),
                          ChQuaternion<>(-0.0, 0.0, 0.0, 1.0));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0197636915985527, 0.250242141495163, 0.2035),
                          ChVector<>(0.209277316591803, -1.22714190001486, 0.031499947373644),
                          ChQuaternion<>(-0.0, 0.0, 0.0, 1.0));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.0154140963254538, 0.25013019465113, 0.2035),
                          ChVector<>(0.467954961389621, -1.22702981403608, 0.0314999008291813),
                          ChQuaternion<>(1.0, -0.0, 0.0, -0.0));
    utils::AddCylinderGeometry(Other.get_ptr(),
                               0.0458692620536463,
                               0.16,
                               ChVector<>(0, -0.517959057661396, 0),
                               ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(),
                               0.00876657438055399,
                               0.16,
                               ChVector<>(-0.381566155791125, -0.515717402690516, 0),
                               ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(),
                               0.0124935472670776,
                               0.16,
                               ChVector<>(0.492612868005681, -0.513858743768626, 0),
                               ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddCylinderGeometry(Other.get_ptr(),
                               0.0246194369069635,
                               0.16,
                               ChVector<>(0.244352690753348, -0.517959057661396, 0),
                               ChQuaternion<>(0.5, -0.5, -0.5, 0.5));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.121733706251767, 0.00876657438055399, 0.165),
                          ChVector<>(-0.324444835447566, -0.623217402690516, 0),
                          ChQuaternion<>(0.515154, 0.0, 0.0, 0.857097));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.121733706251767, 0.00901290521459452, 0.165),
                          ChVector<>(-0.0943376359362894, -0.62333298893951, 0),
                          ChQuaternion<>(0.515154, 0.0, 0.0, -0.857097));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.102920917614647, 0.00880730074763987, 0.165),
                          ChVector<>(0.0461378510072952, -0.614333388045196, 0),
                          ChQuaternion<>(0.739227, 0.0, 0.0, -0.673456));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.110302982718799, 0.009073557274052, 0.165),
                          ChVector<>(0.201202979597005, -0.620296344571127, 0),
                          ChQuaternion<>(0.611299, 0.0, 0.0, -0.791400));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.111244909897975, 0.00991722943912151, 0.165),
                          ChVector<>(0.287088436539427, -0.62133744838765, 0),
                          ChQuaternion<>(0.792869, 0.0, 0.0, -0.609392));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.111244909897976, 0.01, 0.165),
                          ChVector<>(0.463991547662121, -0.621358743768626, 0),
                          ChQuaternion<>(0.609392, 0.0, 0.0, -0.792869));
    // utils::AddCylinderGeometry(Other.get_ptr(),0.0147645993377839,0.447,ChVector<>(0.0545,-0.516955203991318,0.150095419710721),ChQuaternion<>(0.5,0.0,0.0,0.5));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.115826155371625, 0.0147645993377839, 0.447),
                          ChVector<>(0.0545, -0.624455203991318, 0.106974099367162),
                          ChQuaternion<>(0.585725, -0.396140, 0.585725, -0.396140));
    // utils::AddCylinderGeometry(Other.get_ptr(),0.0147645993377839,0.447,ChVector<>(0.0545,-0.516955203991318,-0.150095419710721),ChQuaternion<>(0.5,0.0,0.0,0.5));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.115826155371625, 0.0147645993377839, 0.447),
                          ChVector<>(0.0545, -0.624455203991318, -0.106974099367162),
                          ChQuaternion<>(0.585725, 0.396140, -0.585725, -0.396140));
    utils::AddBoxGeometry(Other.get_ptr(),
                          ChVector<>(0.111244909897976, 0.018292686227323, 0.16500005),
                          ChVector<>(0.455978023535446, -0.619225184316436, 5.00000000014378E-08),
                          ChQuaternion<>(0.0, 0.792869, 0.609392, 0.0));

    // utils::AddTriangleMeshGeometry(Other.get_ptr(), "body_3_1.obj", "body_3_1.obj");
    Other->SetBodyFixed(true);
    system_parallel->AddBody(Other);
  }

  {
    real mass = 3.8327249391131;
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
    Conveyer->SetInertiaXX(ChVector<>(0.029314394817887, 0.0796399719369236, 0.0514634926681141));
    // utils::AddTriangleMeshGeometry(Splitter.get_ptr(), "body_3_1.obj", "body_3_1.obj");
    utils::AddBoxGeometry(Conveyer.get_ptr(),
                          ChVector<>(0.193724352473189, 0.001, 0.15),
                          ChVector<>(0.0785377577161679, 0.554163002941846, 0),
                          QUNIT);
    Conveyer->SetBodyFixed(true);
    system_parallel->AddBody(Conveyer);
  }

  ChSharedPtr<ChLinkEngine> engine_drum = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
  ChSharedPtr<ChLinkEngine> engine_spazzola = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);

  engine_drum->Initialize(Drum, Truss, CS_Drum);
  engine_drum->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
  if (ChSharedPtr<ChFunction_Const> mfun =
          (engine_drum->Get_spe_funct().DynamicCastTo<ChFunction_Const>()))
    mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]
  system_parallel->Add(engine_drum);

  engine_spazzola->Initialize(Spazzola, Truss, CS_Spazzola);
  engine_spazzola->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
  if (ChSharedPtr<ChFunction_Const> mfun =
          (engine_spazzola->Get_spe_funct().DynamicCastTo<ChFunction_Const>()))
    mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]
  system_parallel->Add(engine_spazzola);

  ChSharedPtr<ChMaterialSurface> mat_p = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  mat_p->SetFriction(.0001);    // This is so that it is averaged properly

  gen = new utils::Generator(system_parallel);
  int Id_g = 100;

  MyCreator_metal* callback_metal = new MyCreator_metal;
  callback_metal->systemreference = system_parallel;

  MyCreator_plastic* callback_plastic = new MyCreator_plastic;
  callback_plastic->systemreference = system_parallel;

  utils::MixtureIngredientPtr& m1 = gen->AddMixtureIngredient(utils::SPHERE, 0.5);
  m1->setDefaultSize(0.003);
  m1->setDefaultDensity(8400);
  m1->setDistributionSize(0.003, 1, .002, .005);
  m1->setDefaultMaterialDVI(mat_p);
  m1->SetCallbackPostCreation(callback_metal);

  utils::MixtureIngredientPtr& m2 = gen->AddMixtureIngredient(utils::ELLIPSOID, 0.5);
  m2->setDefaultSize(0.003);
  m2->setDefaultDensity(946);
  m2->setDistributionSize(0.003, 1, .002, .005);
  m2->setDefaultMaterialDVI(mat_p);
  m2->SetCallbackPostCreation(callback_plastic);

  gen->setBodyIdentifier(Id_g);

  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "CES", system_parallel);
  gl_window.SetCamera(ChVector<>(0, 0, -10), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), .1f);

  gl_window.Pause();
  for (int i = 0; i < num_steps; i++) {

    if (gl_window.Active()) {
      if (gl_window.DoStepDynamics(timestep)) {
        RunTimeStep(system_parallel, i);
        TimingOutput(system_parallel);
      }
      gl_window.Render();
    }
  }
  exit(0);

  // gl_window.StartDrawLoop(timestep);

  int file = 0;
  for (int i = 0; i < num_steps; i++) {
    system_parallel->DoStepDynamics(timestep);

    int save_every = 1.0 / timestep / 60.0;    // save data every n steps
    if (i % save_every == 0) {
      std::stringstream ss;
      std::cout << "Frame: " << file << std::endl;
      if (argc == 2) {
        ss << "data/ces/"
           << "/" << file << ".txt";
      } else {
        ss << "data/ces/"
           << "/" << file << ".txt";
      }
      DumpAllObjectsWithGeometryPovray(system_parallel, ss.str(), true);
      file++;
    }
    RunTimeStep(system_parallel, i);
  }

  return 0;
}
