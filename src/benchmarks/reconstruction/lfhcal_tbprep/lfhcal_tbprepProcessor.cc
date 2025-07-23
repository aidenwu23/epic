// Copyright 2023, Friederike Bock
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//  Sections Copyright (C) 2023 Friederike Bock
//  under SPDX-License-Identifier: LGPL-3.0-or-later

#include "lfhcal_tbprepProcessor.h"

#include <DD4hep/Detector.h>
#include <DD4hep/IDDescriptor.h>
#include <DD4hep/Readout.h>
#include <JANA/JApplication.h>
#include <JANA/JEvent.h>
#include <JANA/Services/JGlobalRootLock.h>
#include <RtypesCore.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4hep/CaloHitContributionCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/Vector3d.h>
#include <edm4hep/Vector3f.h>
#include <fmt/core.h>
#include <gsl/pointers>
#include <iostream>
#include <limits>
#include <map>
#include <podio/RelationRange.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>

#include "clusterizer_MA.h"
#include "services/geometry/dd4hep/DD4hep_service.h"
#include "services/log/Log_service.h"
#include "services/rootfile/RootFile_service.h"

// =================================
// Notes & Implementation Details
// =================================
// 1. Event.h compatability:
//    - Fields like RunNumber, Vov, Vop, TimeStamp are ommitted in the simulation, and therefore not in the processor
//
// 2. BeginRun timestamp:
//    - fNanoSec field is excluded because it's empty in the calibrated files
//    - Only the integer seconds field (fSec -> BeginRun) is retained.
//
// 3. Placeholder fields:
//    - ROtype, Trigger Bit (ltpr), and Trigger Primitive (ltrbit) are placeholders for now as per Fredi's instructions
//
// 4. Tree Structure:
//    - The output file structure lacks an intermediate singular branch (typically same name as the tree) between the event
//      tree and leaves, unsure if this will cause significant issues.
// 
// 5. BeamName:
//    - BeamName is not included in the event structure in the clusterizer header file, strings seemed to cause some issues,
//      unsure how to fix it.
//
// 6. For further reference, 
//    - Structures from Event.h are defined as eventsStruct in clusterizer_MA.h
//    - Structures from Tile.h are defined as towersStrct in clusterizer_MA.h

//******************************************************************************************//
// InitWithGlobalRootLock
//******************************************************************************************//
void lfhcal_tbprepProcessor::Init() {
  std::string plugin_name = ("lfhcal_tbprep");

  // ===============================================================================================
  // Get JANA application and seup general variables
  // ===============================================================================================
  auto* app = GetApplication();

  m_log = app->GetService<Log_service>()->logger(plugin_name);

  // Ask service locator for the DD4hep geometry
  auto dd4hep_service = app->GetService<DD4hep_service>();

  // Ask service locator a file to write histograms to
  auto root_file_service = app->GetService<RootFile_service>();

  // Get TDirectory for histograms root file
  auto globalRootLock = app->GetService<JGlobalRootLock>();
  globalRootLock->acquire_write_lock();
  auto* file = root_file_service->GetHistFile();
  globalRootLock->release_lock();

  // ===============================================================================================
  // Create a directory for this plugin. And subdirectories for series of histograms
  // ===============================================================================================
  m_dir_main = file->mkdir(plugin_name.c_str());

  // ===============================================================================================
  // Trees
  // ===============================================================================================
  if (enableTree) {
    event_tree = new TTree("event_tree", "event_tree");
    event_tree->SetDirectory(m_dir_main);

    // ===============================================================================================
    // Event.h structures
    // ===============================================================================================
    beamEnergy                 = 0.0;                   // BeamEnergy (GeV)
    beamPDG                    = 0;                     // BeamID, PDG code of the beam particle 
    eventID                    = 0;                     // EventID
    t_tower_ROtype             = new int[maxNTowers];   // ROtype, currently a place holder (0=Undef, 1=Hgcroc, 2=Caen)
    eventTime                  = 0;                     // BeginRun? TimeStamp? (some sort of UNIX time in seconds)
    beamPosX                   = 0.0;                   // BeamPosX (mm)
    beamPosY                   = 0.0;                   // BeamPosY (mm)

    event_tree->Branch("BeamEnergy", &beamEnergy, "BeamEnergy/D");
    event_tree->Branch("BeamID", &beamPDG, "BeamID/I");
    event_tree->Branch("BeamName", &beamName);          // BeamName, converted from PDG code
    event_tree->Branch("EventID", &eventID, "EventID/I");
    event_tree->Branch("ROtype", t_tower_ROtype, "ROtype/I");
    event_tree->Branch("BeginRun", &eventTime, "BeginRun/L");
    event_tree->Branch("BeamPosX", &beamPosX, "BeamPosX/D");
    event_tree->Branch("BeamPosY", &beamPosY, "BeamPosY/D");

    // ===============================================================================================
    // Tile.h structures
    // ===============================================================================================
    t_ltpr                     = new float[maxNTowers];           // Trigger primitive
    t_ltrbit                   = new unsigned char[maxNTowers];   // Trigger bit
    t_cellID                   = new uint64_t[maxNTowers];        // Unconverted cell IDs
    t_cellID_TB                = new uint64_t[maxNTowers];        // Converted cell IDs
    t_lFHCal_towers_cellE      = new float[maxNTowers];           // Energy (GeV)
    t_lFHCal_towers_cellT      = new float[maxNTowers];           // Time (nanoseconds)

    event_tree->Branch("cell_size", &t_cell_size, "cell_size/I"); // @size for Tiles?
    event_tree->Branch("ltpr", t_ltpr, "ltpr[cell_size]/F");
    event_tree->Branch("ltrbit", t_ltrbit, "ltrbit[cell_size]/b");
    event_tree->Branch("cellID", t_cellID, "cellID[cell_size]/l");
    event_tree->Branch("cellID_TB", t_cellID_TB, "cellID_TB[cell_size]/l");  
    event_tree->Branch("tower_LFHCAL_E", t_lFHCal_towers_cellE, "tower_LFHCAL_E[cell_size]/F"); 
    event_tree->Branch("tower_LFHCAL_T", t_lFHCal_towers_cellT, "tower_LFHCAL_T[cell_size]/F");  
  }

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  auto detector = dd4hep_service->detector();
  std::cout << "--------------------------\nID specification:\n";
  try {
    m_decoder = detector->readout("LFHCALHits").idSpec().decoder();
    std::cout << "1st: " << m_decoder << std::endl;
    iLx      = m_decoder->index("towerx");
    iLy      = m_decoder->index("towery");
    iLz      = m_decoder->index("layerz");
    iPassive = m_decoder->index("passive");
    
    std::cout << "full list: "
              << " " << m_decoder->fieldDescription() << std::endl;
  } catch (...) {
    std::cout << "2nd: " << m_decoder << std::endl;
    m_log->error("readoutClass not in the output");
    throw std::runtime_error("readoutClass not in the output.");
  }
}

//******************************************************************************************//
// ProcessSequential
//******************************************************************************************//
void lfhcal_tbprepProcessor::Process(const std::shared_ptr<const JEvent>& event) {
  // Set eventTime and eventID. Readout Type, trigger bit, and trigger primitive are currently placeholders.
  static eventsStruct tempstructE;
  tempstructE.eventTime = static_cast<long>(std::time(nullptr));
  tempstructE.eventID = static_cast<int>(event->GetEventNumber());
  int tower_ROtype = 0; 
  float ltpr = 0.0f;
  unsigned char ltrbit = 0; 

  // ===============================================================================================
  // process MC particles
  // ===============================================================================================
  const auto& mcParticles = *(event->GetCollection<edm4hep::MCParticle>("MCParticles"));
  tempstructE.beamEnergy = 0.0;
  tempstructE.beamPDG = 0;
  beamName = "unknown";
  tempstructE.beamPosX = 0.0;
  tempstructE.beamPosY = 0.0;
  auto pdgToName = [](int pdg) -> std::string {
    switch (pdg) {
      // Common particles
      case 11: return "e-";
      case -11: return "e+";
      case 12: return "nu_e";
      case -12: return "anti-nu_e";
      case 13: return "mu-";
      case -13: return "mu+";
      case 14: return "nu_mu";
      case -14: return "anti-nu_mu";
      case 15: return "tau-";
      case -15: return "tau+";
      case 16: return "nu_tau";
      case -16: return "anti-nu_tau";
      case 22: return "gamma";
      case 111: return "pi0";
      case 211: return "pi+";
      case -211: return "pi-";
      case 130: return "K0_L";
      case 310: return "K0_S";
      case 321: return "K+";
      case -321: return "K-";
      case 2112: return "neutron";
      case -2112: return "anti-neutron";
      case 2212: return "proton";
      case -2212: return "anti-proton";
      // Rarer baryons and hyperons
      case 3122: return "Lambda";
      case -3122: return "anti-Lambda";
      case 3222: return "Sigma+";
      case -3222: return "anti-Sigma+";
      case 3212: return "Sigma0";
      case -3212: return "anti-Sigma0";
      case 3112: return "Sigma-";
      case -3112: return "anti-Sigma-";
      case 3322: return "Xi0";
      case -3322: return "anti-Xi0";
      case 3312: return "Xi-";
      case -3312: return "anti-Xi-";
      case 3334: return "Omega-";
      case -3334: return "anti-Omega-";
      // Nuclei and special cases
      case 1000010020: return "deuteron";
      case 1000010030: return "triton";
      case 1000020030: return "He3";
      case 1000020040: return "alpha";
      default: return "unknown/not added";
    }
  };
  for (auto mcparticle : mcParticles) {
    if (mcparticle.getGeneratorStatus() != 1) continue;
    tempstructE.beamEnergy = mcparticle.getEnergy();
    tempstructE.beamPDG = mcparticle.getPDG();
    beamName = pdgToName(tempstructE.beamPDG);
    auto vertex = mcparticle.getVertex();
    tempstructE.beamPosX = vertex.x;
    tempstructE.beamPosY = vertex.y;
    break; 
  }

  // ===============================================================================================
  // read rec hits & fill structs
  // ===============================================================================================
  const auto& recHits = *(event->GetCollection<edm4eic::CalorimeterHit>(nameRecHits));
  int nCaloHitsRec    = 0;
  std::vector<towersStrct> input_tower_recSav;
  // process rec hits
  for (const auto caloHit : recHits) {
    uint64_t cellID = caloHit.getCellID();
    float energy    = caloHit.getEnergy();
    float time      = caloHit.getTime();

    auto detector_module_x = m_decoder->get(cellID, "moduleIDx");
    auto detector_module_y = m_decoder->get(cellID, "moduleIDy");
    auto detector_passive  = m_decoder->get(cellID, iPassive);
    auto detector_layer_x  = m_decoder->get(cellID, iLx);
    auto detector_layer_y  = m_decoder->get(cellID, iLy);
    int detector_layer_rz  = -1;
    int detector_layer    = m_decoder->get(cellID, iLz);
    if (isLFHCal) {
      detector_layer_rz = m_decoder->get(cellID, "rlayerz");
    }
    if (detector_passive > 0) {
      continue;
    }

    // Process test beam cell IDs
    uint64_t cellID_TB = 0;
    if (isLFHCal) {
      cellID_TB = ReconstructCellID(detector_module_x, detector_module_y,
                                     detector_layer_x, detector_layer_y, detector_layer);
    }

    nCaloHitsRec++;

    //loop over input_tower_recSav and find if there is already a tower with the same cellID
    bool found = false;
    for (auto& tower : input_tower_recSav) {
      if (tower.cellID == static_cast<decltype(tower.cellID)>(cellID)) {
        tower.energy += energy;
        found = true;
        break;
      }
    }
    if (!found) {
      towersStrct tempstructT;
      tempstructT.cellID  = cellID;
      tempstructT.energy  = energy;
      tempstructT.time    = time;
      if (isLFHCal) {
        tempstructT.cellIDz = detector_layer_rz;
        tempstructT.cellID_TB = cellID_TB;
      }
      // Set the placeholder variables in the struct
      tempstructT.tower_ROtype = tower_ROtype;
      tempstructT.ltpr = ltpr;
      tempstructT.ltrbit = ltrbit;
      input_tower_recSav.push_back(tempstructT);
    }
  }

  m_log->trace("LFHCal mod: nCaloHits rec {}", nCaloHitsRec);
  
  // ===============================================================================================
  // Write event tree & clean-up variables
  // ===============================================================================================
  if (enableTree) {
    t_cell_size = (int)input_tower_recSav.size();

    for (int iCell = 0; iCell < (int)input_tower_recSav.size(); iCell++) {
      t_lFHCal_towers_cellE[iCell]      = (float)input_tower_recSav.at(iCell).energy;
      t_lFHCal_towers_cellT[iCell]      = (float)input_tower_recSav.at(iCell).time;
      t_cellID[iCell] = input_tower_recSav.at(iCell).cellID;
      t_cellID_TB[iCell] = input_tower_recSav.at(iCell).cellID_TB;
      t_tower_ROtype[iCell] = input_tower_recSav.at(iCell).tower_ROtype;
      t_ltpr[iCell] = input_tower_recSav.at(iCell).ltpr;
      t_ltrbit[iCell] = input_tower_recSav.at(iCell).ltrbit;
    }

    event_tree->Fill();

    t_cell_size = 0;

    for (Int_t itow = 0; itow < maxNTowers; itow++) {
      t_lFHCal_towers_cellE[itow]      = 0;
      t_lFHCal_towers_cellT[itow]      = 0;
      t_cellID[itow] = 0;
      t_cellID_TB[itow] = 0;
      t_tower_ROtype[itow] = 0;
      t_ltpr[itow] = 0.0f;
      t_ltrbit[itow] = 0;
    }
  }
}

//******************************************************************************************//
// Finish
//******************************************************************************************//
void lfhcal_tbprepProcessor::Finish() {
  std::cout << "------> LFHCal processor finished" << std::endl;
  // Do any final calculations here.

  if (enableTree) {
    delete[] t_lFHCal_towers_cellE;
    delete[] t_lFHCal_towers_cellT;
    delete[] t_cellID;
    delete[] t_cellID_TB;
    delete[] t_tower_ROtype;
    delete[] t_ltpr;
    delete[] t_ltrbit;
  }
}

//******************************************************************************************//
// Reconstruct cell ID based on bit info in NewStructure
//******************************************************************************************//
uint64_t lfhcal_tbprepProcessor::ReconstructCellID(auto detector_module_x, auto detector_module_y,
                                                   auto detector_layer_x, auto detector_layer_y, int detector_layer) {
  if (!isLFHCal) return 0; // Only works for LFHCal
  
  auto moduleID = detector_module_x * 4 + detector_module_y;
  uint64_t cellID_TB = (moduleID << 9) + (detector_layer_y << 8) + (detector_layer_x << 6) + detector_layer;

  return cellID_TB;
}