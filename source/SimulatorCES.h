#ifndef SIMULATORCES_H
#define SIMULATORCES_H

///
/// HEADER for the simulator of the Corona Electrostatic Separator
/// for the WEEE project
/// Authors I. Critelli, A. Tasora, 2014
///

#include "chrono_parallel/physics/ChSystemParallel.h"
//#include "physics/ChConveyor.h"
#include "physics/ChBodyAuxRef.h"
#include "core/ChFileutils.h"

#include "core/ChRealtimeStep.h"
#include "core/ChMath.h"
#include "core/ChDistribution.h"
#include "collision/ChCCollisionSystemBullet.h"
#include "particlefactory/ChParticleEmitter.h"
#include "particlefactory/ChParticleRemover.h"
#include "particlefactory/ChParticleProcessor.h"

#include <fstream>
#include <sstream>

#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"

#include "ElectricParticleProperty.h"
#include "ElectricForcesCES.h"
#include "ParserEmitter.h"
#include "ParserElectricForcesCES.h"
#include "ProcessFlow.h"

// Use the namespace of Chrono

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::particlefactory;

using namespace std;

class SimulatorCES {

 public:
  //
  // DATA of the simulator
  //

  // The ChronoENGINE physical system
  ChSystemParallelDVI mphysicalSystem;

  ElectricForcesCES ces_forces;    // this contains data for computing the CES electric forces

  ChParticleEmitter emitter;
  ChSharedPtr<ChRandomParticlePositionRectangleOutlet> emitter_positions;
  ChSharedPtr<ChRandomParticleAlignment> emitter_rotations;

  double drumspeed_rpm;      // [rpm]
  double drumspeed_radss;    //[rad/s]

  // material surfaces
  float surface_drum_friction;
  float surface_drum_rolling_friction;
  float surface_drum_spinning_friction;
  float surface_drum_restitution;
  float surface_plate_friction;
  float surface_plate_rolling_friction;
  float surface_plate_spinning_friction;
  float surface_plate_restitution;
  ChSharedPtr<ChMaterialSurface> surface_particles;

  double max_particle_age;

  double xnozzlesize;
  double znozzlesize;

  double flowmeter_xmin;
  double flowmeter_xmax;
  double flowmeter_width;
  double flowmeter_y;
  int flowmeter_bins;

  // Coordinate systems with position and rotation of important items in the
  // simulator. These are initializad with constant values, but if loading the
  // SolidWorks model, they will be changed accordingly to what is found in the CAD
  // file (see later, where the SolidWorks model is parsed).

  ChCoordsys<> conveyor_csys;
  ChCoordsys<> drum_csys;
  ChCoordsys<> nozzle_csys;
  ChCoordsys<> Splitter1_csys;
  ChCoordsys<> Splitter2_csys;
  ChCoordsys<> Spazzola_csys;

  // set as true for saving log files each n frames
  bool save_dataset;
  bool save_irrlicht_screenshots;
  bool save_POV_screenshots;
  int saveEachNframes;

  bool irr_cast_shadows;
  int totframes;
  bool init_particle_speed;
  double particle_magnification;    // for larger visualization of particle
  std::string solidworks_py_modelfile;
  std::string results_file;
  double timestep;
  double Tmax;
  bool splitters_collide;

  ///
  /// Create the SimulatorCES
  /// and initialize member data
  ///
  SimulatorCES() {
    // Initialize member data:

    solidworks_py_modelfile =
        "../CAD_conveyor/conveyor_Ida";    // note! do not add ".py" after the filename

    results_file = "output/results.dat";

    drumspeed_rpm = 44.8;
    drumspeed_radss = drumspeed_rpm * ((2.0 * CH_C_PI) / 60.0);    //[rad/s]

    // sphrad = 0.38e-3;
    // sphrad2 = 0.25e-3;
    // sphrad3 = 0.794e-3;

    surface_drum_friction = 0.5f;
    surface_drum_rolling_friction = 0;
    surface_drum_spinning_friction = 0;
    surface_drum_restitution = 0;
    surface_plate_friction = 0.2f;
    surface_plate_rolling_friction = 0;
    surface_plate_spinning_friction = 0;
    surface_plate_restitution = 0;

    surface_particles = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
    surface_particles->SetFriction(0.2f);
    surface_particles->SetRollingFriction(0);
    surface_particles->SetSpinningFriction(0);
    surface_particles->SetRestitution(0);

    max_particle_age = 2;

    // nozzle sizes
    znozzlesize = 0.182;    //**from CAD, nozzle z width
    xnozzlesize = 0.1;      //**from CAD, nozzle x width

    flowmeter_xmin = 0.28;
    flowmeter_xmax = flowmeter_xmin + 0.3;
    flowmeter_y = -0.1;
    flowmeter_width = 0.2;
    flowmeter_bins = 25;

    // Init coordinate systems with position and rotation of important items in the
    // simulator. These are initializad with constant values, but if loading the
    // SolidWorks model, they will be changed accordingly to what is found in the CAD
    // file (see later, where the SolidWorks model is parsed).
    /*
    //***ALEX disabled because initialized by SolidWorks file, anyway
    double conv_thick = 0.01;
    double conveyor_length = 0.6;
    conveyor_csys	= ChCoordsys<>( ChVector<>(0, -conv_thick, 0) ) ; // default position
    drum_csys		= ChCoordsys<>( ChVector<>(conveyor_length/2,
    -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
    nozzle_csys		= ChCoordsys<>( ChVector<>(0, 0.01, 0) ); // default position
    Splitter1_csys	= ChCoordsys<>( ChVector<>(conveyor_length/2+0.2,
    -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
    Splitter2_csys	= ChCoordsys<>( ChVector<>(conveyor_length/2+0.4,
    -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
    Spazzola_csys	= ChCoordsys<>( ChVector<>(conveyor_length/2-0.10,
    -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
    */

    // set as true for saving log files each n frames
    save_dataset = false;
    save_irrlicht_screenshots = false;
    save_POV_screenshots = false;
    saveEachNframes = 3;

    irr_cast_shadows = true;
    totframes = 0;
    init_particle_speed = true;
    particle_magnification = 3;    // for larger visualization of particle
    timestep = 0.001;
    Tmax = 0.1;
    splitters_collide = true;

    // Set small collision envelopes for objects that will be created from now on..
    ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);    // 0.002
    ChCollisionModel::SetDefaultSuggestedMargin(0.0005);     // 0.0008
    // Set contact breaking/merging tolerance of Bullet:
    ChCollisionSystemBullet::SetContactBreakingThreshold(0.001);

    // Important! dt is small, and particles are small, so it's better to keep this small...
    // not needed in INT_TASORA, only for INT_ANITESCU
    mphysicalSystem.SetMaxPenetrationRecoverySpeed(0.15);
    mphysicalSystem.SetMinBounceSpeed(0.1);

    // In the following there is a default initialization of the
    // particle creation system, based on ChParticleEmitter.
    // This is a default configuration, that is _overridden_ if you
    // call ParseSettings() and load a settings.ces file that contain different
    // configurations for the emitter.

    // ---Initialize the randomizer for positions
    emitter_positions = ChSharedPtr<ChRandomParticlePositionRectangleOutlet>(
        new ChRandomParticlePositionRectangleOutlet);
    emitter_positions->OutletWidth() = 0.1;       // default x outlet size, from CAD;
    emitter_positions->OutletHeight() = 0.182;    // default y outlet size, from CAD;
    emitter.SetParticlePositioner(emitter_positions);

    // ---Initialize the randomizer for alignments
    emitter_rotations =
        ChSharedPtr<ChRandomParticleAlignmentUniform>(new ChRandomParticleAlignmentUniform);
    emitter.SetParticleAligner(emitter_rotations);

    // ---Initialize the randomizer for creations, with statistical distribution

    // Create a ChRandomShapeCreator object (ex. here for metal particles)
    ChSharedPtr<ChRandomShapeCreatorSpheres> mcreator_metal(new ChRandomShapeCreatorSpheres);
    mcreator_metal->SetDiameterDistribution(
        ChSmartPtr<ChMinMaxDistribution>(new ::ChMinMaxDistribution(0.002, 0.003)));

    // Optional: define a callback to be exectuted at each creation of a metal particle:
    class MyCreator_metal : public ChCallbackPostCreation {
      // Here do custom stuff on the just-created particle:
     public:
      virtual void PostCreation(ChSharedPtr<ChBody> mbody,
                                ChCoordsys<> mcoords,
                                ChRandomShapeCreator& mcreator) {
        // Attach some optional visualization stuff
        // ChSharedPtr<ChTexture> mtexture(new ChTexture);
        // mtexture->SetTextureFilename("../objects/pinkwhite.png");
        // mbody->AddAsset(mtexture);
        ChSharedPtr<ChColorAsset> mvisual(new ChColorAsset);
        mvisual->SetColor(ChColor(0.9f, 0.4f, 0.2f));
        mbody->AddAsset(mvisual);
        // Attach a custom asset. It will hold electrical properties
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
    MyCreator_metal* callback_metal = new MyCreator_metal;
    callback_metal->systemreference = &this->mphysicalSystem;
    mcreator_metal->SetCallbackPostCreation(callback_metal);

    // Create a ChRandomShapeCreator object (ex. here for metal particles)
    ChSharedPtr<ChRandomShapeCreatorSpheres> mcreator_plastic(new ChRandomShapeCreatorSpheres);
    mcreator_plastic->SetDiameterDistribution(
        ChSmartPtr<ChMinMaxDistribution>(new ::ChMinMaxDistribution(0.002, 0.002)));

    // Optional: define a callback to be exectuted at each creation of a plastic particle:
    class MyCreator_plastic : public ChCallbackPostCreation {
      // Here do custom stuff on the just-created particle:
     public:
      virtual void PostCreation(ChSharedPtr<ChBody> mbody,
                                ChCoordsys<> mcoords,
                                ChRandomShapeCreator& mcreator) {
        // Attach some optional visualization stuff
        // ChSharedPtr<ChTexture> mtexture(new ChTexture);
        // mtexture->SetTextureFilename("../objects/bluwhite.png");
        // mbody->AddAsset(mtexture);
        ChSharedPtr<ChColorAsset> mvisual(new ChColorAsset);
        mvisual->SetColor(ChColor(0.3f, 0.6f, 0.7f));
        mbody->AddAsset(mvisual);
        // Attach a custom asset. It will hold electrical properties
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

        ++particlecounter;
        mbody->SetIdentifier(particlecounter);
      }
      // here put custom data of the callback
      ChSystem* systemreference;
      int particlecounter;
    };
    MyCreator_plastic* callback_plastic = new MyCreator_plastic;
    callback_plastic->systemreference = &this->mphysicalSystem;
    callback_plastic->particlecounter = 0;
    mcreator_plastic->SetCallbackPostCreation(callback_plastic);

    // Create a parent ChRandomShapeCreator that 'mixes' the two generators above,
    // mixing them with a given percentual:
    ChSharedPtr<ChRandomShapeCreatorFromFamilies> mcreatorTot(new ChRandomShapeCreatorFromFamilies);
    mcreatorTot->AddFamily(mcreator_metal, 0.4);      // 1st creator family, with percentual
    mcreatorTot->AddFamily(mcreator_plastic, 0.4);    // 2nd creator family, with percentual
    mcreatorTot->Setup();

    // Finally, tell to the emitter that it must use the 'mixer' above:
    emitter.SetParticleCreator(mcreatorTot);

    // ---Initialize the randomizer for velocities, with statistical distribution

    ChSharedPtr<ChRandomParticleVelocityConstantDirection> mvelo(
        new ChRandomParticleVelocityConstantDirection);
    mvelo->SetDirection(-VECT_Y);
    mvelo->SetModulusDistribution(0.0);

    emitter.SetParticleVelocity(mvelo);
  }

  ///
  /// Parser
  /// - load settings from a .ces input file, with simulator settings in JSON format

  bool ParseSettings(const char* filename) {
    try {
      // Prepare  input stream and copy to char* buffer
      ChStreamInAsciiFile settingfile(filename);
      std::stringstream buffer;
      buffer << settingfile.GetFstream().rdbuf();
      std::string mstr = buffer.str();
      const char* stringbuffer = mstr.c_str();

      rapidjson::Document document;
      document.Parse<0>(stringbuffer);
      if (document.HasParseError()) {
        std::string errstrA((const char*)(&stringbuffer[ChMax(document.GetErrorOffset() - 10, 0)]));
        errstrA.resize(10);
        std::string errstrB((const char*)(&stringbuffer[document.GetErrorOffset()]));
        errstrB.resize(20);
        throw(ChException("the file has bad JSON syntax," + std::string(document.GetParseError()) +
                          " \n\n[...]" + errstrA + " <--- " + errstrB + "[...]\n"));
      }
      if (!document.IsObject())
        throw(ChException("the file is not a valid JSON document"));

      char* token;

      token = (char*)"solidworks_exported_model";
      if (document.HasMember(token)) {
        if (!document[token].IsString()) {
          throw(ChException("Invalid filename string after '" + std::string(token) + "'"));
        }
        this->solidworks_py_modelfile = document[token].GetString();
      }
      token = (char*)"results_file";
      if (document.HasMember(token)) {
        if (!document[token].IsString()) {
          throw(ChException("Invalid filename string after '" + std::string(token) + "'"));
        }
        this->results_file = document[token].GetString();
      }
      token = (char*)"drum_rpm";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->drumspeed_rpm = document[token].GetDouble();
        this->drumspeed_radss = drumspeed_rpm * ((2.0 * CH_C_PI) / 60.0);    //[rad/s]
      }
      token = (char*)"save_each_Nsteps";
      if (document.HasMember(token)) {
        if (!document[token].IsInt()) {
          throw(ChException("Invalid integer number after '" + std::string(token) + "'"));
        }
        this->saveEachNframes = document[token].GetInt();
      }
      token = (char*)"save_dataset";
      if (document.HasMember(token)) {
        if (!document[token].IsBool()) {
          throw(ChException("Invalid true/false flag after '" + std::string(token) + "'"));
        }
        this->save_dataset = document[token].GetBool();
      }
      token = (char*)"save_irrlicht_screenshots";
      if (document.HasMember(token)) {
        if (!document[token].IsBool()) {
          throw(ChException("Invalid true/false flag after '" + std::string(token) + "'"));
        }
        this->save_irrlicht_screenshots = document[token].GetBool();
      }
      token = (char*)"save_POV_screenshots";
      if (document.HasMember(token)) {
        if (!document[token].IsBool()) {
          throw(ChException("Invalid true/false flag after '" + std::string(token) + "'"));
        }
        this->save_POV_screenshots = document[token].GetBool();
      }
      token = (char*)"timestep";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->timestep = document[token].GetDouble();
      }
      token = (char*)"Tmax";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->Tmax = document[token].GetDouble();
      }
      token = (char*)"surface_drum_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_drum_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_drum_rolling_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_drum_rolling_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_drum_spinning_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_drum_spinning_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_drum_restitution";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_drum_restitution = (float)document[token].GetDouble();
      }
      token = (char*)"surface_plate_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_plate_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_plate_rolling_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_plate_rolling_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_plate_spinning_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_plate_spinning_friction = (float)document[token].GetDouble();
      }
      token = (char*)"surface_plate_restitution";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->surface_plate_restitution = (float)document[token].GetDouble();
      }
      token = (char*)"surface_particles_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        // this->surface_particles_friction = (float)document[token].GetDouble();
        this->surface_particles->SetFriction((float)document[token].GetDouble());
      }
      token = (char*)"surface_particles_rolling_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        // this->surface_particles_rolling_friction = (float)document[token].GetDouble();
        this->surface_particles->SetRollingFriction((float)document[token].GetDouble());
      }
      token = (char*)"surface_particles_spinning_friction";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        // this->surface_particles_spinning_friction = (float)document[token].GetDouble();
        this->surface_particles->SetSpinningFriction((float)document[token].GetDouble());
      }
      token = (char*)"surface_particles_restitution";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        // this->surface_particles_restitution = (float)document[token].GetDouble();
        this->surface_particles->SetRestitution((float)document[token].GetDouble());
      }
      token = (char*)"max_particle_age";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->max_particle_age = (float)document[token].GetDouble();
      }
      token = (char*)"default_collision_envelope";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        ChCollisionModel::SetDefaultSuggestedEnvelope(document[token].GetDouble());
      }
      token = (char*)"default_collision_margin";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        ChCollisionModel::SetDefaultSuggestedMargin(document[token].GetDouble());
      }
      token = (char*)"default_contact_breaking";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        ChCollisionSystemBullet::SetContactBreakingThreshold(document[token].GetDouble());
      }
      token = (char*)"max_penetration_recovery_speed";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        mphysicalSystem.SetMaxPenetrationRecoverySpeed(document[token].GetDouble());
      }
      token = (char*)"min_bounce_speed";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        mphysicalSystem.SetMinBounceSpeed(document[token].GetDouble());
      }
      token = (char*)"flowmeter_xmin";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->flowmeter_xmin = document[token].GetDouble();
      }
      token = (char*)"flowmeter_xmax";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->flowmeter_xmax = document[token].GetDouble();
      }
      token = (char*)"flowmeter_y";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->flowmeter_y = document[token].GetDouble();
      }
      token = (char*)"flowmeter_width";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->flowmeter_width = document[token].GetDouble();
      }
      token = (char*)"flowmeter_bins";
      if (document.HasMember(token)) {
        if (!document[token].IsNumber()) {
          throw(ChException("Invalid number after '" + std::string(token) + "'"));
        }
        this->flowmeter_bins = document[token].GetInt();
      }
      token = (char*)"splitters_collide";
      if (document.HasMember(token)) {
        if (!document[token].IsBool()) {
          throw(ChException("Invalid true/false flag after '" + std::string(token) + "'"));
        }
        this->splitters_collide = document[token].GetBool();
      }
      token = (char*)"CES_forces";
      if (document.HasMember(token)) {
        // Parse the settings of the emitter, emitter positioner etc.
        ParserElectricForcesCES::Parse(this->ces_forces, document[token]);
      }
      token = (char*)"emitter";
      if (document.HasMember(token)) {
        // Parse the settings of the emitter, emitter positioner etc.
        ParserEmitter::Parse(this->emitter,
                             this->mphysicalSystem,
                             this->emitter_positions,
                             this->emitter_rotations,
                             document[token]);
      }
    }
    catch (ChException me) {
      GetLog() << "ERROR loading settings file: \n   " << filename << "\n Reason: " << me.what()
               << "\n\n";
      throw(ChException("Error in parsing."));
    }

    return true;
  }

  //
  // This can be added to store the trajectory on a per-particle basis.
  //
  class ParticleTrajectory : public ChAsset {
   public:
    std::list<ChVector<> > positions;
    std::list<ChVector<> > speeds;
    unsigned int max_points;

    ParticleTrajectory() { max_points = 80; }
  };

  ///
  /// Function that deletes old debris
  /// (to avoid infinite creation that fills memory)
  ///
  void purge_debris(ChSystem& mysystem, double max_age = 5.0) {
    for (unsigned int i = 0; i < mysystem.Get_bodylist()->size(); i++) {
      ChBody* abody = (*mysystem.Get_bodylist())[i];

      bool to_delete = false;

      // Fetch the ElectricParticleProperty asset from the list of
      // assets that have been attached to the object, and retrieve the
      // custom data that have been stored. ***ALEX
      for (unsigned int na = 0; na < abody->GetAssets().size(); na++) {
        ChSharedAssetPtr myasset = abody->GetAssetN(na);
        if (myasset.IsType<ElectricParticleProperty>()) {
          ChSharedPtr<ElectricParticleProperty> electricproperties =
              myasset.DynamicCastTo<ElectricParticleProperty>();
          double particle_birthdate = electricproperties->birthdate;
          double particle_age = mysystem.GetChTime() - particle_birthdate;
          if (particle_age > max_age) {
            to_delete = true;
          }
        }
      }

      if (to_delete) {
        abody->AddRef();                            // dirty trick to convert basic pointer to..
        ChSharedPtr<ChBody> mysharedbody(abody);    // ..shared pointer

        mysystem.Remove(mysharedbody);

        // mysharedbody->RemoveRef(); //***NOT needed - previously needed cause always Add() to POV
        // exporter..
        i--;    // this because if deleted, the rest of the array is shifted back one position..
      }
    }
  }

  ///
  /// Main function of the simulator.
  /// Initialize the simulation, and
  /// Performs the simulation
  /// by running the loop of time integration
  ///
  int simulate() {

    // 3) fetch coordinate values and objects from what was imported from CAD

    // ChCoordsys<> conveyor_csys = CSYSNORM;

    ChSharedPtr<ChMarker> my_marker = mphysicalSystem.SearchMarker("centro_nastro");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find centro_nastro marker from its name in the C::E system! \n";
    else
      conveyor_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    //****Ida

    my_marker = mphysicalSystem.SearchMarker("Splitter1");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find Splitter1 marker from its name in the C::E system! \n";
    else
      Splitter1_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    my_marker = mphysicalSystem.SearchMarker("Splitter2");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find Splitter2 marker from its name in the C::E system! \n";
    else
      Splitter2_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    my_marker = mphysicalSystem.SearchMarker("Spazzola");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find Spazzola marker from its name in the C::E system! \n";
    else
      Spazzola_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    //***Ida

    my_marker = mphysicalSystem.SearchMarker("centro_nozzle");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find centro_nozzle marker from its name in the C::E system! \n";
    else
      nozzle_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    emitter_positions->Outlet() = nozzle_csys;
    emitter_positions->Outlet().rot.Q_from_AngAxis(CH_C_PI_2, VECT_X);    // rotate outlet 90� on x

    my_marker = mphysicalSystem.SearchMarker("centro_cilindro");
    if (my_marker.IsNull())
      GetLog() << "Error: cannot find centro_cilindro marker from its name in the C::E system! \n";
    else
      drum_csys = my_marker->GetAbsCoord();    // fetch both pos and rotation of CAD

    // fetch mrigidBodyDrum pointer! will be used for changing the friction, the collision family,
    // and later to create the motor
    ChSharedPtr<ChBodyAuxRef> mrigidBodyDrum =
        mphysicalSystem.Search("drum-1").DynamicCastTo<ChBodyAuxRef>();
    if (mrigidBodyDrum.IsNull())
      GetLog() << "ERROR: cannot find drum-1 from its name in the C::E system! ! \n";
    else {
      mrigidBodyDrum->GetCollisionModel()->SetFamily(3);
      mrigidBodyDrum->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
      mrigidBodyDrum->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
      mrigidBodyDrum->SetFriction(surface_drum_friction);
      mrigidBodyDrum->SetImpactC(surface_drum_restitution);
      mrigidBodyDrum->SetRollingFriction(surface_drum_rolling_friction);
      mrigidBodyDrum->SetSpinningFriction(surface_drum_spinning_friction);
    }

    //***Ida

    ChSharedPtr<ChBodyAuxRef> mrigidBodySplitter1 =
        mphysicalSystem.Search("Splitter-10").DynamicCastTo<ChBodyAuxRef>();
    if (mrigidBodySplitter1.IsNull())
      GetLog() << "ERROR: cannot find Splitter-10 from its name in the C::E system! ! \n";
    else {
      mrigidBodySplitter1->SetBodyFixed(true);
      mrigidBodySplitter1->GetCollisionModel()->SetFamily(3);    // rivedere
      mrigidBodySplitter1->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(
          1);    // rivedere
      mrigidBodySplitter1->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(
          2);    // rivedere
      mrigidBodySplitter1->SetFriction(0.1f);
      mrigidBodySplitter1->SetCollide(this->splitters_collide);    // deactivate collision?
    }

    ChSharedPtr<ChBodyAuxRef> mrigidBodySplitter2 =
        mphysicalSystem.Search("Splitter2-1").DynamicCastTo<ChBodyAuxRef>();
    if (mrigidBodySplitter2.IsNull())
      GetLog() << "ERROR: cannot find Splitter2-1 from its name in the C::E system! ! \n";
    else {
      mrigidBodySplitter2->SetBodyFixed(true);
      mrigidBodySplitter2->GetCollisionModel()->SetFamily(3);    // rivedere
      mrigidBodySplitter2->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(
          1);    // rivedere
      mrigidBodySplitter2->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(
          2);    // rivedere
      mrigidBodySplitter2->SetFriction(0.1f);
      mrigidBodySplitter2->SetCollide(this->splitters_collide);    // deactivate collision?
    }

    ChSharedPtr<ChBodyAuxRef> mrigidBodySpazzola =
        mphysicalSystem.Search("Spazzola-1").DynamicCastTo<ChBodyAuxRef>();
    if (mrigidBodySpazzola.IsNull())
      GetLog() << "ERROR: cannot find Spazzola-1 from its name in the C::E system! ! \n";
    else {
      mrigidBodySpazzola->GetCollisionModel()->SetFamily(1);                             // rivedere
      mrigidBodySpazzola->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);    // rivedere
      mrigidBodySpazzola->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);    // rivedere
      mrigidBodySpazzola->SetFriction(0.9f);
    }

    ChSharedPtr<ChBodyAuxRef> mrigidBodyConveyor =
        mphysicalSystem.Search("conveyor-1").DynamicCastTo<ChBodyAuxRef>();
    if (mrigidBodyConveyor.IsNull())
      GetLog() << "ERROR: cannot find conveyor from its name in the C::E system! ! \n";
    else {
      mrigidBodyConveyor->GetCollisionModel()->SetFamily(2);
      mrigidBodyConveyor->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
      mrigidBodyConveyor->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
      mrigidBodyConveyor->SetFriction(surface_plate_friction);
      mrigidBodyConveyor->SetImpactC(surface_plate_restitution);
      mrigidBodyConveyor->SetRollingFriction(surface_plate_rolling_friction);
      mrigidBodyConveyor->SetSpinningFriction(surface_plate_spinning_friction);
    }

    //
    // Create a truss (absolute fixed reference body, for connecting the rotating cyl.)
    //

    ChSharedPtr<ChBody> mtruss(new ChBody);
    mtruss->SetBodyFixed(true);

    // Finally, do not forget to add the body to the system:
    mphysicalSystem.Add(mtruss);

    //**Ida

    ChSharedPtr<ChBody> mtruss2(new ChBody);
    mtruss2->SetBodyFixed(true);

    // Finally, do not forget to add the body to the system:
    mphysicalSystem.Add(mtruss2);

    //
    // Create a motor constraint between the cylinder and the truss
    //

    ChSharedPtr<ChLinkEngine> mengine;

    if (!mrigidBodyDrum.IsNull()) {
      mengine = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
      ChSharedPtr<ChBody> mdrum(mrigidBodyDrum);
      mengine->Initialize(mdrum, mtruss, drum_csys);

      mengine->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
      if (ChSharedPtr<ChFunction_Const> mfun =
              (mengine->Get_spe_funct().DynamicCastTo<ChFunction_Const>()))
        mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]

      // Finally, do not forget to add the body to the system:
      mphysicalSystem.Add(mengine);
    }

    //***Ida

    ChSharedPtr<ChLinkEngine> mengine2;

    if (!mrigidBodySpazzola.IsNull()) {
      mengine2 = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
      ChSharedPtr<ChBody> mSpazzola(mrigidBodySpazzola);
      mengine2->Initialize(mSpazzola, mtruss2, Spazzola_csys);

      mengine2->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
      if (ChSharedPtr<ChFunction_Const> mfun =
              (mengine2->Get_spe_funct().DynamicCastTo<ChFunction_Const>()))
        mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]

      // Finally, do not forget to add the body to the system:
      mphysicalSystem.Add(mengine2);
    }

    //
    // What to do by default on ALL newly created particles?
    // A callback executed at each particle creation can be attached to the emitter..
    //

    // a- define a class that implement your custom PostCreation method...
    class MyCreatorForAll : public ChCallbackPostCreation {
     public:
      virtual void PostCreation(ChSharedPtr<ChBody> mbody,
                                ChCoordsys<> mcoords,
                                ChRandomShapeCreator& mcreator) {
        // Set the friction properties (using a shared ChSurfaceMaterial
        mbody->SetMaterialSurface(asurface_material);

        // Disable gyroscopic forces for increased integrator stabilty
        mbody->SetNoGyroTorque(true);
      }

      ChSharedPtr<ChMaterialSurface> asurface_material;
    };
    // b- create the callback object...
    MyCreatorForAll* mcreation_callback = new MyCreatorForAll;
    // c- set callback own data that he might need...
    mcreation_callback->asurface_material = this->surface_particles;
    // d- attach the callback to the emitter!
    emitter.SetCallbackPostCreation(mcreation_callback);

    //
    // PROCESS THE FLOW with these tools:
    //

    // Create also a ChParticleProcessor configured as a
    // counter of particles that flow into a rectangle with a statistical distribution to plot:
    //  -create the trigger:
    double flowmeter_length = this->flowmeter_xmax - this->flowmeter_xmin;
    ChSharedPtr<ChParticleEventFlowInRectangle> distrrectangle(
        new ChParticleEventFlowInRectangle(flowmeter_length, flowmeter_width));
    distrrectangle->rectangle_csys =
        ChCoordsys<>(drum_csys.pos + ChVector<>(this->flowmeter_xmin + 0.5 * flowmeter_length,
                                                this->flowmeter_y,
                                                0),          // position of center rectangle
                     Q_from_AngAxis(-CH_C_PI_2, VECT_X));    // rotate rectangle so that its Z is up
    distrrectangle->margin = 0.05;
    //  -create the counter, with 20x10 resolution of sampling, on x y
    //    This is defined in ProcessFlow.h and distinguishes plastic from metal
    ChSharedPtr<ProcessFlow> countdistribution(new ProcessFlow(this->flowmeter_bins, 1));
    //  -create the processor and plug in the trigger and the counter:
    ChParticleProcessor processor_distribution;
    processor_distribution.SetEventTrigger(distrrectangle);
    processor_distribution.SetParticleEventProcessor(countdistribution);

    // Create a remover, i.e. an object that takes care
    // of removing particles that are inside or outside some volume.
    // The fact that particles are handled with shared pointers means that,
    // after they are removed from the ChSystem, they are also automatically
    // deleted if no one else is referencing them.
    ChSharedPtr<ChParticleEventFlowInRectangle> distrrectangle2(
        new ChParticleEventFlowInRectangle(0.20, 0.30));
    distrrectangle2->rectangle_csys = distrrectangle->rectangle_csys;
    distrrectangle2->margin = 0.05;
    ChSharedPtr<ChParticleProcessEventRemove> removal_event(new ChParticleProcessEventRemove);
    ChParticleProcessor processor_remover;
    processor_remover.SetEventTrigger(distrrectangle2);
    processor_remover.SetParticleEventProcessor(removal_event);

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    int savenum = 0;

    ChFileutils::MakeDirectory("screenshots");
    ChFileutils::MakeDirectory("output");

    while (true) {
      if (mphysicalSystem.GetChTime() > this->Tmax)
        break;

      totframes++;

      // Apply the forces caused by electrodes of the CES machine:

      ces_forces.apply_forces(
          &mphysicalSystem,    // contains all bodies
          drum_csys,           // pos and rotation of axis of drum (not rotating reference!)
          drumspeed_radss,     // speed of drum
          totframes);

      // Continuosly create debris that fall on the conveyor belt
      this->emitter.EmitParticles(mphysicalSystem, timestep);    //***TEST***

      GetLog() << "Total mass=" << this->emitter.GetTotCreatedMass() << "   "
               << "Total n.part=" << this->emitter.GetTotCreatedParticles() << "   "
               << "Average kg/s=" << this->emitter.GetTotCreatedMass() / mphysicalSystem.GetChTime()
               << "\n";

      // Use the processor to count particle flow in the rectangle section:
      processor_distribution.ProcessParticles(mphysicalSystem);

      // Continuosly check if some particle must be removed:
      processor_remover.ProcessParticles(mphysicalSystem);

      // Maybe the user played with the slider and changed the speed of drum...
      if (!mengine.IsNull())
        if (ChSharedPtr<ChFunction_Const> mfun =
                (mengine->Get_spe_funct().DynamicCastTo<ChFunction_Const>()))
          mfun->Set_yconst(-drumspeed_radss);    // angular speed in [rad/s]
    }

    // At the end ot the T max simulation time,
    // save output distributions to disk (non normalized for unit area/volume),
    // they can be a nxm matrix of 2d bins or a n-vector of 1d bins
    GetLog() << "\n saving output distributions... \n ";

    ChStreamOutAsciiFile file_for_metal("out_distribution_metal.txt");
    countdistribution->mmass_metal.StreamOUTdenseMatlabFormat(file_for_metal);
    ChStreamOutAsciiFile file_for_plastic("out_distribution_plastic.txt");
    countdistribution->mmass_plastic.StreamOUTdenseMatlabFormat(file_for_plastic);

    GetLog() << "\n Simulation Terminated. \n ";

    return 0;
  }

};    // end of class

#endif
