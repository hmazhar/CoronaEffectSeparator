#ifndef PARSERELECTRICFORCESCES_H
#define PARSERELECTRICFORCESCES_H

#include "ElectricForcesCES.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"
#include <fstream>

using namespace chrono;

/// Class for parsing data of the CES emitter
class ParserElectricForcesCES {
 public:
  /// Just a simple static function here, call as ParserElectricForcesCES::Parse(...)
  static void Parse(ElectricForcesCES& ces_forces, rapidjson::Value& mval) {
    char* token;
    if (!mval.IsObject()) {
      throw(ChException("Invalid object after 'CES_forces'"));
    }

    token = (char*)"U";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.U = mval[token].GetDouble();
    }
    token = (char*)"L";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.L = mval[token].GetDouble();
    }
    token = (char*)"alpha_deg";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.alpha = (CH_C_PI / 180) * mval[token].GetDouble();
    }
    token = (char*)"drum_diameter";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.drum_diameter = mval[token].GetDouble();
    }
    token = (char*)"drum_width";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.drum_width = mval[token].GetDouble();
    }
    token = (char*)"electrode_diameter";
    if (mval.HasMember(token)) {
      if (!mval[token].IsNumber()) {
        throw(ChException("Invalid number after '" + std::string(token) + "'"));
      }
      ces_forces.electrode_diameter = mval[token].GetDouble();
    }
  }
};

#endif
