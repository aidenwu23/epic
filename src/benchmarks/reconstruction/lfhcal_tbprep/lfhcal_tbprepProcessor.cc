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

#include "clusterizer_MA.h"
#include "services/geometry/dd4hep/DD4hep_service.h"
#include "services/log/Log_service.h"
#include "services/rootfile/RootFile_service.h"

/* 
Notes/TODO
  - Current output file excludes RunNumber, Vov, Vop, TimeStamp because these are not available in simulation
  - fNanoSec in BeginRun is excluded because its equivalent branch in calibrated files is empty (leaving only fSec, which is simplified to BeginRun)
  - Meaning and usage of @size for Tiles remain unclear and are thus not included
    - For now its likely equivalent to tower_LFHCAL_N, not confirmed
  - ROtype, triggerBit, and triggerPrimitive are currently placeholders, implemented as per Fredi's instructions
  - Most lFHCal_tower code (from lfhcal_studies) are commented out but retained for potential future use
  - The output file structure lacks an intermediate singular branch between the event tree and leaves, unsure if this is will cause issues
*/

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

    // Beam energy and PDG
    beamEnergy = 0.0;
    beamPDG = 0;
    event_tree->Branch("BeamEnergy", &beamEnergy, "BeamEnergy/D");
    event_tree->Branch("BeamID", &beamPDG, "BeamID/I");
    event_tree->Branch("BeamName", &beamName);
    // Event ID branch
    eventID = 0;
    event_tree->Branch("EventID", &eventID, "EventID/I");
    // Placeholder for readout type (0=undef, 1=hgcroc, 2=caen)
    readoutType = 0;
    event_tree->Branch("ROtype", &readoutType, "ROtype/I");
    // Event timestamp branch (UNIX time in seconds)
    eventTime = 0;
    event_tree->Branch("BeginRun", &eventTime, "BeginRun/L");
    // Particle gun positions
    beamPosX = 0.0;
    beamPosY = 0.0;
    event_tree->Branch("BeamPosX", &beamPosX, "BeamPosX/D");
    event_tree->Branch("BeamPosY", &beamPosY, "BeamPosY/D");
    // Placeholders for trigger bit and trigger primitive 
    triggerBit = 0;
    triggerPrimitive = 0;
    event_tree->Branch("triggerBit", &triggerBit, "triggerBit/I");
    event_tree->Branch("triggerPrimitive", &triggerPrimitive, "triggerPrimitive/I");
    // Cell IDs
    event_tree->Branch("tower_LFHCAL_N", &t_lFHCal_towers_N, "tower_LFHCAL_N/I"); // from lfhcal_studies
    t_cellID_TB = new uint64_t[maxNTowers];
    t_cellID = new uint64_t[maxNTowers];
    event_tree->Branch("cellID", t_cellID, "cellID[tower_LFHCAL_N]/l");
    event_tree->Branch("cellID_TB", t_cellID_TB, "cellID_TB[tower_LFHCAL_N]/l");

    // lFHCal_towers taken from lfhcal_studies, kept for potential future use
    /*
    t_lFHCal_towers_cellE      = new float[maxNTowers];
    t_lFHCal_towers_cellT      = new float[maxNTowers];
    t_lFHCal_towers_cellIDx    = new short[maxNTowers];
    t_lFHCal_towers_cellIDy    = new short[maxNTowers];
    t_lFHCal_towers_cellIDz    = new short[maxNTowers];
    t_lFHCal_towers_cellTrueID = new int[maxNTowers];

    event_tree->Branch("tower_LFHCAL_E", t_lFHCal_towers_cellE, "tower_LFHCAL_E[tower_LFHCAL_N]/F");
    event_tree->Branch("tower_LFHCAL_T", t_lFHCal_towers_cellT, "tower_LFHCAL_T[tower_LFHCAL_N]/F");
    event_tree->Branch("tower_LFHCAL_ix", t_lFHCal_towers_cellIDx,
                       "tower_LFHCAL_ix[tower_LFHCAL_N]/S");
    event_tree->Branch("tower_LFHCAL_iy", t_lFHCal_towers_cellIDy,
                       "tower_LFHCAL_iy[tower_LFHCAL_N]/S");
    event_tree->Branch("tower_LFHCAL_iz", t_lFHCal_towers_cellIDz,
                       "tower_LFHCAL_iz[tower_LFHCAL_N]/S");
    event_tree->Branch("tower_LFHCAL_trueID", t_lFHCal_towers_cellTrueID,
                       "tower_LFHCAL_trueID[tower_LFHCAL_N]/I");
    */
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
  // Set eventTime and eventID. ReadoutType, triggerBit, and triggerPrimitive are currently placeholders.
  eventTime = static_cast<long>(std::time(nullptr));
  eventID = static_cast<int>(event->GetEventNumber());
  readoutType = 0;
  triggerBit = 0;
  triggerPrimitive = 0;

  // ===============================================================================================
  // process MC particles
  // ===============================================================================================
  const auto& mcParticles = *(event->GetCollection<edm4hep::MCParticle>("MCParticles"));
  beamEnergy = 0.0;
  beamPDG = 0;
  beamName = "unknown";
  beamPosX = 0.0;
  beamPosY = 0.0;
  auto pdgToName = [](int pdg) -> std::string {
    switch (pdg) {
      case 11: return "e-";
      case -11: return "e+";
      case 13: return "mu-";
      case -13: return "mu+";
      case 22: return "gamma";
      case 211: return "pi+";
      case -211: return "pi-";
      case 2212: return "proton";
      case -2212: return "anti-proton";
      case 2112: return "neutron";
      case -2112: return "anti-neutron";
      default: return "unknown/not added";
    }
  };
  for (auto mcparticle : mcParticles) {
    if (mcparticle.getGeneratorStatus() != 1) continue;
    beamEnergy = mcparticle.getEnergy();
    beamPDG = mcparticle.getPDG();
    beamName = pdgToName(beamPDG);
    auto vertex = mcparticle.getVertex();
    beamPosX = vertex.x;
    beamPosY = vertex.y;
    break; 
  }

  // ===============================================================================================
  // read rec hits & fill structs
  // ===============================================================================================
  const auto& recHits = *(event->GetCollection<edm4eic::CalorimeterHit>(nameRecHits));
  int nCaloHitsRec    = 0;
  std::vector<towersStrct> input_tower_recSav;
  std::vector<uint64_t> input_cellID_TB;
  std::vector<uint64_t> input_cellID;
  // process rec hits
  for (const auto caloHit : recHits) {
    float x         = caloHit.getPosition().x / 10.;
    float y         = caloHit.getPosition().y / 10.;
    float z         = caloHit.getPosition().z / 10.;
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

    // calc cell IDs
    long cellIDx = -1;
    long cellIDy = -1;
    if (isLFHCal) {
      cellIDx = 54LL * 2 - detector_module_x * 2 + detector_layer_x;
      cellIDy = 54LL * 2 - detector_module_y * 2 + detector_layer_y;
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
      // Kept for potential future use
      tempstructT.energy  = energy;
      tempstructT.time    = time;
      tempstructT.posx    = x;
      tempstructT.posy    = y;
      tempstructT.posz    = z;
      tempstructT.cellID  = cellID;
      tempstructT.cellIDx = cellIDx;
      tempstructT.cellIDy = cellIDy;
      if (isLFHCal) {
        tempstructT.cellIDz = detector_layer_rz;
      }
      tempstructT.tower_trueID = 0; //TODO how to get trueID?
      input_tower_recSav.push_back(tempstructT);

      // Cell IDs
      input_cellID_TB.push_back(cellID_TB);
      input_cellID.push_back(cellID);
    }
  }
  m_log->trace("LFHCal mod: nCaloHits rec {}", nCaloHitsRec);
  
  // ===============================================================================================
  // Write event tree & clean-up variables
  // ===============================================================================================
  if (enableTree) {
    t_lFHCal_towers_N = (int)input_tower_recSav.size();
    // See comment in Init
    /*
    for (int iCell = 0; iCell < (int)input_tower_recSav.size(); iCell++) {
      t_lFHCal_towers_cellE[iCell]      = (float)input_tower_recSav.at(iCell).energy;
      t_lFHCal_towers_cellT[iCell]      = (float)input_tower_recSav.at(iCell).time;
      t_lFHCal_towers_cellIDx[iCell]    = (short)input_tower_recSav.at(iCell).cellIDx;
      t_lFHCal_towers_cellIDy[iCell]    = (short)input_tower_recSav.at(iCell).cellIDy;
      t_lFHCal_towers_cellIDz[iCell]    = (short)input_tower_recSav.at(iCell).cellIDz;
      t_lFHCal_towers_cellTrueID[iCell] = (int)input_tower_recSav.at(iCell).tower_trueID;
      t_cellID[iCell] = input_cellID.at(iCell);
      t_cellID_TB[iCell] = input_cellID_TB.at(iCell);
    }
    */
    for (int iCell = 0; iCell < (int)input_tower_recSav.size(); iCell++) {
      t_cellID[iCell] = input_cellID.at(iCell);
      t_cellID_TB[iCell] = input_cellID_TB.at(iCell);
    }

    event_tree->Fill();

    t_lFHCal_towers_N = 0;
    // See comment in Init
    /*
    for (Int_t itow = 0; itow < maxNTowers; itow++) {
      t_lFHCal_towers_cellE[itow]      = 0;
      t_lFHCal_towers_cellT[itow]      = 0;
      t_lFHCal_towers_cellIDx[itow]    = 0;
      t_lFHCal_towers_cellIDy[itow]    = 0;
      t_lFHCal_towers_cellIDz[itow]    = 0;
      t_lFHCal_towers_cellTrueID[itow] = 0;
      t_cellID[itow] = 0;
      t_cellID_TB[itow] = 0;
    }
    */
    for (Int_t itow = 0; itow < maxNTowers; itow++) {
      t_cellID[itow] = 0;
      t_cellID_TB[itow] = 0;
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
    // See comment in Init
    /*
    delete[] t_lFHCal_towers_cellE;
    delete[] t_lFHCal_towers_cellT;
    delete[] t_lFHCal_towers_cellIDx;
    delete[] t_lFHCal_towers_cellIDy;
    delete[] t_lFHCal_towers_cellIDz;
    delete[] t_lFHCal_towers_cellTrueID;
    */
    delete[] t_cellID;
    delete[] t_cellID_TB;
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