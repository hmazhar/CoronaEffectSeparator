#ifndef CHRONOMODELS_INOUT_H
#define CHRONOMODELS_INOUT_H

// Include bullet collision model
#include "collision/ChCModelBullet.h"
#include "collision/ChCCollisionSystemBullet.h"

// Include Parallel system and descriptor
#include "chrono_parallel/physics/ChSystemParallel.h"

using namespace chrono;
using namespace chrono::collision;
using namespace geometry;

#include <iostream>
#include <sstream>
#include <cstring>
#include <zlib.h>

std::string get_binary_file_contents(const std::string filename) {
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in) {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return (contents);
  }
  throw(errno);
}
std::string getfile_contents(const std::string filename) {
  std::ifstream in(filename, std::ios::in);
  if (in) {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return (contents);
  }
  throw(errno);
}

class CSVGen {
 public:
  CSVGen() {
    delim = ",";
    binary = false;
    gz_file = 0;
  }
  ~CSVGen() {}

  void OpenFile(std::string filename, bool bin = false) {
    if (bin) {
      binary = true;
      gz_file = gzopen(filename.c_str(), "wb");
    } else {
      ofile.open(filename.c_str(), std::ios::out);
    }
  }

  void CloseFile() {
    if (binary) {
      unsigned long int file_size = sizeof(char) * ss.str().size();
      gzwrite(gz_file,
              (void*)&file_size,
              sizeof(file_size));  // writing size of file
      gzwrite(gz_file, (void*)(ss.str().data()), file_size);
      gzclose(gz_file);
    } else {
      ofile << ss.str();
      ofile.close();
    }
  }
  template <class T>
  void operator<<(const T& token) {
    WriteToken(token);
  }
  void WriteToken(real token) { ss << token << delim; }
  void WriteToken(real2 token) { ss << token.x << delim << token.y << delim; }
  void WriteToken(real3 token) {
    ss << token.x << delim << token.y << delim << token.z << delim;
  }
  void WriteToken(real4 token) {
    ss << token.w << delim << token.x << delim << token.y << delim << token.z
       << delim;
  }
  void WriteToken(std::string token) { ss << token << delim; }
  void endline() { ss << std::endl; }

  std::string delim;
  std::ofstream ofile;
  std::stringstream ss;
  gzFile gz_file;
  bool binary;
};


void DumpAllObjectsWithGeometryPovray(ChSystem* mSys,
                                      std::string filename,
                                      bool binary = false) {
  CSVGen csv_output;
  csv_output.OpenFile(filename.c_str(), binary);
  for (int i = 0; i < mSys->Get_bodylist()->size(); i++) {
    ChBody* abody = mSys->Get_bodylist()->at(i);
    const Vector pos = abody->GetPos();
    const Vector vel = abody->GetPos_dt();
    Quaternion b_rot = abody->GetRot();
    Vector pos_final, rad_final;
    ShapeType type = SPHERE;

    for (int j = 0; j < abody->GetAssets().size(); j++) {
      Quaternion rot(1, 0, 0, 0);
      ChSharedPtr<ChAsset> asset = abody->GetAssets().at(j);
      ChVisualization* visual_asset = ((ChVisualization*)(asset.get_ptr()));
      Vector center = visual_asset->Pos;
      center = b_rot.Rotate(center);
      pos_final = pos + center;
      Quaternion lrot = visual_asset->Rot.Get_A_quaternion();
      rot = b_rot % lrot;
      rot.Normalize();
      if (asset.IsType<ChSphereShape>()) {
        ChSphereShape* sphere_shape = ((ChSphereShape*)(asset.get_ptr()));
        float radius = sphere_shape->GetSphereGeometry().rad;
        rad_final.x = radius;
        rad_final.y = radius;
        rad_final.z = radius;
        type = SPHERE;
      }

      else if (asset.IsType<ChEllipsoidShape>()) {
        ChEllipsoidShape* ellipsoid_shape =
            ((ChEllipsoidShape*)(asset.get_ptr()));
        rad_final = ellipsoid_shape->GetEllipsoidGeometry().rad;
        type = ELLIPSOID;
      } else if (asset.IsType<ChBoxShape>()) {
        ChBoxShape* box_shape = ((ChBoxShape*)(asset.get_ptr()));
        rad_final = box_shape->GetBoxGeometry().Size;
        type = BOX;
      } else if (asset.IsType<ChCylinderShape>()) {
        ChCylinderShape* cylinder_shape = ((ChCylinderShape*)(asset.get_ptr()));
        double rad = cylinder_shape->GetCylinderGeometry().rad;
        double height = cylinder_shape->GetCylinderGeometry().p1.y -
                        cylinder_shape->GetCylinderGeometry().p2.y;
        rad_final.x = rad;
        rad_final.y = height;
        rad_final.z = rad;
        type = CYLINDER;
      } else if (asset.IsType<ChConeShape>()) {
        ChConeShape* cone_shape = ((ChConeShape*)(asset.get_ptr()));
        rad_final.x = cone_shape->GetConeGeometry().rad.x;
        rad_final.y = cone_shape->GetConeGeometry().rad.y;
        rad_final.z = cone_shape->GetConeGeometry().rad.z;
        type = CONE;
      }

      csv_output << R3(pos_final.x, pos_final.y, pos_final.z);
      csv_output << R4(rot.e0, rot.e1, rot.e2, rot.e3);
      csv_output << R3(vel.x, vel.y, vel.z);

      if (asset.IsType<ChSphereShape>()) {
        csv_output << type;
        csv_output << rad_final.x;
        csv_output.endline();
      } else if (asset.IsType<ChEllipsoidShape>()) {
        csv_output << type;
        csv_output << R3(rad_final.x, rad_final.y, rad_final.z);
        csv_output.endline();
      } else if (asset.IsType<ChBoxShape>()) {
        csv_output << type;
        csv_output << R3(rad_final.x, rad_final.y, rad_final.z);
        csv_output.endline();
      } else if (asset.IsType<ChCylinderShape>()) {
        csv_output << type;
        csv_output << R2(rad_final.x, rad_final.y);
        csv_output.endline();
      } else if (asset.IsType<ChConeShape>()) {
        csv_output << type;
        csv_output << R2(rad_final.x, rad_final.y);
        csv_output.endline();
      } else {
        csv_output << -1;
        csv_output.endline();
      }
    }
  }
  csv_output.CloseFile();
}


void TimingOutput(ChSystem* mSys) {
  double TIME = mSys->GetChTime();
  double STEP = mSys->GetTimerStep();
  double BROD = mSys->GetTimerCollisionBroad();
  double NARR = mSys->GetTimerCollisionNarrow();
  double LCP = mSys->GetTimerLcp();
  double UPDT = mSys->GetTimerUpdate();
  double RESID = 0;
  int REQ_ITS = 0;
  int BODS = mSys->GetNbodies();
  int CNTC = mSys->GetNcontacts();
  if (ChSystemParallelDVI* parallel_sys =
          dynamic_cast<ChSystemParallelDVI*>(mSys)) {
    RESID =
        ((ChLcpSolverParallelDVI*)(mSys->GetLcpSolverSpeed()))->GetResidual();
    REQ_ITS = ((ChLcpSolverParallelDVI*)(mSys->GetLcpSolverSpeed()))
                  ->GetTotalIterations();
    BODS = parallel_sys->GetNbodies();
    CNTC = parallel_sys->GetNcontacts();
  }

  printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n",
         TIME,
         STEP,
         BROD,
         NARR,
         LCP,
         UPDT,
         BODS,
         CNTC,
         REQ_ITS,
         RESID);
}

void TimingOutputVerbose(ChSystemParallelDVI* parallel_sys) {
  double TIME = parallel_sys->GetChTime();
  double STEP = parallel_sys->GetTimerStep();
  double BROD = parallel_sys->GetTimerCollisionBroad();
  double NARR = parallel_sys->GetTimerCollisionNarrow();
  double LCP = parallel_sys->GetTimerLcp();
  double UPDT = parallel_sys->GetTimerUpdate();
  double RESID = 0;
  int REQ_ITS = 0;
  int BODS = parallel_sys->GetNbodies();
  int CNTC = parallel_sys->GetNcontacts();

  RESID = ((ChLcpSolverParallelDVI*)(parallel_sys->GetLcpSolverSpeed()))
              ->GetResidual();
  REQ_ITS = ((ChLcpSolverParallelDVI*)(parallel_sys->GetLcpSolverSpeed()))
                ->GetTotalIterations();
  BODS = parallel_sys->GetNbodies();
  CNTC = parallel_sys->GetNcontacts();
  printf("Time: %7.4f\n", parallel_sys->GetChTime());
  printf("Step: %7.4f\n", parallel_sys->GetTimerStep());
  printf("brod: %7.4f\n", parallel_sys->GetTimerCollisionBroad());
  printf("narr: %7.4f\n", parallel_sys->GetTimerCollisionNarrow());
  printf("lcp : %7.4f\n", parallel_sys->GetTimerLcp());
  printf("updt: %7.4f\n", parallel_sys->GetTimerUpdate());
  printf("resi: %7.4f\n", RESID);
  printf("iter: %7d\n", REQ_ITS);
  printf("bods: %7d\n", BODS);
  printf("cntc: %7d\n", CNTC);
  printf("bldD: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildD"));
  printf("blDA: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildDAllocate"));
  printf("blDC: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildDCompute"));
  printf("bldE: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildE"));
  printf("bldN: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildN"));
  printf("bldM: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("BuildM"));
  printf("bldb: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("Buildb"));
  printf("shur: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime("ShurProduct"));
  printf("solA: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime(
             "ChSolverParallel_solverA"));
  printf("solB: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime(
             "ChSolverParallel_solverB"));
  printf("solC: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime(
             "ChSolverParallel_solverC"));
  printf("solD: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime(
             "ChSolverParallel_solverD"));
  printf("solE: %7.4f\n",
         parallel_sys->data_manager->system_timer.GetTime(
             "ChSolverParallel_solverE"));
}
#endif
