// Copyright 2023, Friederike Bock
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//  Sections Copyright (C) 2023 Friederike Bock
//  under SPDX-License-Identifier: LGPL-3.0-or-later

#include <DDSegmentation/BitFieldCoder.h>
#include <JANA/JEvent.h>
#include <JANA/JEventProcessor.h>
#include <JANA/Utils/JTypeInfo.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <spdlog/logger.h>
#include <memory>
#include <string>

class lfhcal_tbprepProcessor : public JEventProcessor {
public:
  lfhcal_tbprepProcessor() { SetTypeName(NAME_OF_THIS); }

  void Init() override;
  //     void InitWithGlobalRootLock() override;
  //     void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
  void Process(const std::shared_ptr<const JEvent>& event) override;
  void Finish() override;
  //     void FinishWithGlobalRootLock() override;
  TDirectory* m_dir_main;
  
  // Simple event tree for per-event analysis, mostly unused
  bool enableTree = true;
  TTree* event_tree;
  const int maxNTowers = 65000;
  static const int maxCells= 65000;
  uint64_t t_lFHCal_cells_cellID[maxCells];
  float t_lFHCal_cells_energy[maxCells];
  float* t_lFHCal_towers_cellE;
  float* t_lFHCal_towers_cellT;

  int t_cell_size;

  uint64_t* t_cellID_TB;
  uint64_t* t_cellID;
  int* t_tower_ROtype;         // readout type
  float* t_ltpr;               // trigger primitive
  unsigned char* t_ltrbit;     // trigger bit

  double beamEnergy;
  int beamPDG;
  std::string beamName;
  double beamPosX;
  double beamPosY;
  int eventID;
  long eventTime; // UNIX time in seconds
  int tileSize; // TODO: what is this?

  // ReconstructCellID function declaration
  uint64_t ReconstructCellID(auto detector_module_x, auto detector_module_y,
                             auto detector_layer_x, auto detector_layer_y, int detector_layer);

  bool isLFHCal           = true;
  std::shared_ptr<spdlog::logger> m_log;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  std::string nameSimHits       = "LFHCALHits";
  std::string nameRecHits       = "LFHCALRecHits";
  std::string nameClusters      = "LFHCALClusters";
  std::string nameProtoClusters = "LFHCALIslandProtoClusters";
  short iPassive;
  short iLx;
  short iLy;
  short iLz;
};
