// Copyright 2023, Friederike Bock
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//  Sections Copyright (C) 2023 Friederike Bock
//  under SPDX-License-Identifier: LGPL-3.0-or-later

#include <vector>
#include <TVector3.h>
#include <string>

struct towersStrct {
  towersStrct()
      : energy(0)
      , time(0)
      , cellID(0)
      , cellID_TB(0)
      , cellIDz(-1)
      , tower_ROtype(-1)
      , ltpr(0.0)
      , ltrbit(0) {}
  float energy;
  float time;
  int cellID;
  uint64_t cellID_TB;
  int cellIDz;
  int tower_ROtype;
  float ltpr;
  unsigned char ltrbit;
};

struct eventsStruct {
  double beamEnergy;
  int beamPDG;
  int eventID;
  long eventTime;
  double beamPosX;
  double beamPosY;
};
