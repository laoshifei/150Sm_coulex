#ifndef SET_H
#define SET_H

// Historical constants only. The current build executable uses CLI options
// and the complete pars/time_alignment profile, not this header.

#define OPENS3ATT
#define OPENS3COR

constexpr int nMaxS3Ring = 24;
constexpr int nMaxS3Sector = 32;

constexpr int TimeWindow = 200; // ns

constexpr double CutSiEnergy = 2500.; // keV 
constexpr double CutS3EnergyDiff = 2500.; // keV 

#endif
