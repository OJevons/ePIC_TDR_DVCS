// ---------------------------------------------------------------------------
//  DVCSBinning.hh
//
//  Flexible (ragged) binning reader for the ePIC DVCS analysis.
//
//  Parses a plain-text binning file (see bins.txt) that defines, per Q2 bin,
//  a B0-region binning and an RP-region binning. Within each region the |t|
//  edges may depend on the xB slice ("ragged" binning), and the two regions
//  may use different xB/|t| schemes. The NUMBER of bins in Q2, xB and |t| is
//  therefore taken entirely from the file and can vary between configurations.
//
//  Expected file format (blank lines and any non-matching header/separator
//  lines are ignored):
//
//    Q2 bin 1 - [0.99 - 1.4] GeV2
//    (1) B0 region -
//    |t| = {0.45, 2.}
//    x = {0., 1.}
//    (2) RP region -
//    xB = {0, 0.00032, 0.00046, 1}
//    xB: {0, 0.00032}      |t| = {0, 0.08, 0.14, 0.28, 2}
//    xB: {0.00032, 0.00046}|t| = {0, 0.08, 0.14, 0.26, 2}
//    xB: {0.00046, 1}      |t| = {0, 0.1, 0.16, 0.26, 2}
//    ... (repeat for each Q2 bin) ...
//
//  Notes:
//   - Q2 bins are assumed contiguous: the edge list is built from the low edge
//     of the first bin and the high edge of every bin.
//   - The B0 region is expected to have a single xB slice ("x = {..}") with a
//     single "|t| = {..}" definition.
//   - In the RP region, one "xB: {a, b} |t| = {..}" line is expected per xB
//     slice, in the same order as the "xB = {..}" edge list.
// ---------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <iostream>

class DVCSBinning {
 public:
  enum Region { kB0 = 0, kRP = 1, kNRegions = 2 };

  DVCSBinning() {}

  // Load and parse a binning file. Returns true on success.
  bool load(const std::string& path);

  bool ok() const { return fLoaded; }

  // --- Q2 (shared between regions) ---
  int    nQ2()          const { return (int)fQ2Edges.size() - 1; }
  double q2Low (int iq2) const { return fQ2Edges[iq2];   }
  double q2High(int iq2) const { return fQ2Edges[iq2+1]; }

  // --- xB (per region, per Q2 bin) ---
  int    nXB   (int reg, int iq2)          const { return (int)fXBEdges[reg][iq2].size() - 1; }
  double xBLow (int reg, int iq2, int ixB) const { return fXBEdges[reg][iq2][ixB];   }
  double xBHigh(int reg, int iq2, int ixB) const { return fXBEdges[reg][iq2][ixB+1]; }

  // --- |t| (per region, per Q2 bin, per xB slice) ---
  int    nT   (int reg, int iq2, int ixB)         const { return (int)fTEdges[reg][iq2][ixB].size() - 1; }
  double tLow (int reg, int iq2, int ixB, int it) const { return fTEdges[reg][iq2][ixB][it];   }
  double tHigh(int reg, int iq2, int ixB, int it) const { return fTEdges[reg][iq2][ixB][it+1]; }

  // --- Bin finders (return -1 if the value falls outside the scheme) ---
  int findQ2(double q2) const;
  int findXB(int reg, int iq2, double xB) const;
  int findT (int reg, int iq2, int ixB, double t) const;

  // Pretty-print the loaded scheme (for logging / sanity checks).
  void print() const;

 private:
  bool fLoaded{false};
  std::vector<double> fQ2Edges;                                                 // [nQ2+1]
  std::array<std::vector<std::vector<double>>, kNRegions>              fXBEdges; // [reg][iq2]       -> xB edges
  std::array<std::vector<std::vector<std::vector<double>>>, kNRegions> fTEdges;  // [reg][iq2][ixB]  -> |t| edges

  bool validate() const;
  static std::vector<double> parseBraces(const std::string& s); // extract numbers inside the first {...}
  static std::string trim(const std::string& s);
};

// ---------------------------------------------------------------------------
//  Implementation (header-only: inline)
// ---------------------------------------------------------------------------
inline std::string DVCSBinning::trim(const std::string& s) {
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

inline std::vector<double> DVCSBinning::parseBraces(const std::string& s) {
  std::vector<double> v;
  const size_t lb = s.find('{');
  const size_t rb = s.find('}', lb);
  if (lb == std::string::npos || rb == std::string::npos) return v;
  std::string inside = s.substr(lb + 1, rb - lb - 1);
  std::stringstream ss(inside);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    const std::string t = trim(tok);
    if (t.empty()) continue;
    try { v.push_back(std::stod(t)); }
    catch (...) { std::cerr << "[DVCSBinning] WARNING: could not parse number '" << t << "'\n"; }
  }
  return v;
}

inline bool DVCSBinning::load(const std::string& path) {
  fLoaded = false;
  fQ2Edges.clear();
  for (int r = 0; r < kNRegions; ++r) { fXBEdges[r].clear(); fTEdges[r].clear(); }

  std::ifstream in(path);
  if (!in.is_open()) {
    std::cerr << "[DVCSBinning] ERROR: cannot open binning file '" << path << "'\n";
    return false;
  }

  int iq2    = -1;  // current Q2 bin index
  int region = -1;  // current region (kB0 / kRP)

  std::string raw;
  while (std::getline(in, raw)) {
    const std::string line = trim(raw);
    if (line.empty()) continue;

    // Q2 bin header: "Q2 bin N - [lo - hi] GeV2"
    if (line.rfind("Q2 bin", 0) == 0 && line.find('[') != std::string::npos) {
      const size_t lb = line.find('[');
      const size_t rb = line.find(']', lb);
      if (rb == std::string::npos) continue;
      const std::string inside = line.substr(lb + 1, rb - lb - 1); // "lo - hi"
      const size_t dash = inside.find('-');
      if (dash == std::string::npos) continue;
      double lo = 0., hi = 0.;
      try {
        lo = std::stod(trim(inside.substr(0, dash)));
        hi = std::stod(trim(inside.substr(dash + 1)));
      } catch (...) { std::cerr << "[DVCSBinning] WARNING: bad Q2 header '" << line << "'\n"; continue; }

      if (fQ2Edges.empty()) fQ2Edges.push_back(lo);
      fQ2Edges.push_back(hi);
      iq2    = (int)fQ2Edges.size() - 2; // index of the bin just opened
      region = -1;
      for (int r = 0; r < kNRegions; ++r) { fXBEdges[r].resize(iq2 + 1); fTEdges[r].resize(iq2 + 1); }
      continue;
    }

    // Region markers
    if (line.find("High-t") != std::string::npos) { region = kB0; continue; }
    if (line.find("Low-t") != std::string::npos) { region = kRP; continue; }

    // Ignore anything before we know both the Q2 bin and region
    if (iq2 < 0 || region < 0) continue;

    // RP per-slice line: "xB: {a, b}   |t| = {..}"  -> take the |t| group
    if (line.rfind("xB:", 0) == 0) {
      const size_t tp = line.find("|t|");
      if (tp == std::string::npos) continue;
      fTEdges[region][iq2].push_back(parseBraces(line.substr(tp)));
      continue;
    }
    // RP xB edge list: "xB = {..}"
    if (line.rfind("xB", 0) == 0) {
      fXBEdges[region][iq2] = parseBraces(line);
      continue;
    }
    // B0 |t| edges: "|t| = {..}" (single xB slice)
    if (line.rfind("|t|", 0) == 0) {
      fTEdges[region][iq2].clear();
      fTEdges[region][iq2].push_back(parseBraces(line));
      continue;
    }
    // B0 xB edges: "x = {..}"
    if (line.rfind("x", 0) == 0) {
      fXBEdges[region][iq2] = parseBraces(line);
      continue;
    }
  }

  fLoaded = validate();
  return fLoaded;
}

inline bool DVCSBinning::validate() const {
  if (nQ2() < 1) { std::cerr << "[DVCSBinning] ERROR: no Q2 bins parsed\n"; return false; }
  const char* rn[2] = {"B0", "RP"};
  for (int r = 0; r < kNRegions; ++r) {
    if ((int)fXBEdges[r].size() != nQ2() || (int)fTEdges[r].size() != nQ2()) {
      std::cerr << "[DVCSBinning] ERROR: region " << rn[r] << " Q2-size mismatch\n";
      return false;
    }
    for (int q = 0; q < nQ2(); ++q) {
      const int nx = (int)fXBEdges[r][q].size() - 1;
      if (nx < 1) {
        std::cerr << "[DVCSBinning] ERROR: region " << rn[r] << " Q2 bin " << q << " has no xB bins\n";
        return false;
      }
      if ((int)fTEdges[r][q].size() != nx) {
        std::cerr << "[DVCSBinning] ERROR: region " << rn[r] << " Q2 bin " << q
                  << ": number of |t| slices (" << fTEdges[r][q].size()
                  << ") != number of xB bins (" << nx << ")\n";
        return false;
      }
      for (int x = 0; x < nx; ++x) {
        if ((int)fTEdges[r][q][x].size() < 2) {
          std::cerr << "[DVCSBinning] ERROR: region " << rn[r] << " Q2 bin " << q
                    << " xB slice " << x << " has no |t| bins\n";
          return false;
        }
      }
    }
  }
  return true;
}

inline int DVCSBinning::findQ2(double q2) const {
  for (int i = 0; i < nQ2(); ++i)
    if (q2 >= fQ2Edges[i] && q2 < fQ2Edges[i + 1]) return i;
  return -1;
}

inline int DVCSBinning::findXB(int reg, int iq2, double xB) const {
  if (iq2 < 0 || iq2 >= nQ2()) return -1;
  const std::vector<double>& e = fXBEdges[reg][iq2];
  for (int i = 0; i + 1 < (int)e.size(); ++i)
    if (xB >= e[i] && xB < e[i + 1]) return i;
  return -1;
}

inline int DVCSBinning::findT(int reg, int iq2, int ixB, double t) const {
  if (iq2 < 0 || iq2 >= nQ2())        return -1;
  if (ixB < 0 || ixB >= nXB(reg, iq2)) return -1;
  const std::vector<double>& e = fTEdges[reg][iq2][ixB];
  for (int i = 0; i + 1 < (int)e.size(); ++i)
    if (t >= e[i] && t < e[i + 1]) return i;
  return -1;
}

inline void DVCSBinning::print() const {
  const char* rn[2] = {"B0", "RP"};
  std::cout << "[DVCSBinning] Loaded scheme: " << nQ2() << " Q2 bin(s)\n";
  for (int q = 0; q < nQ2(); ++q) {
    std::cout << "  Q2 bin " << q << "  [" << q2Low(q) << ", " << q2High(q) << "] GeV^2\n";
    for (int r = 0; r < kNRegions; ++r) {
      std::cout << "    " << rn[r] << " region: " << nXB(r, q) << " xB bin(s)\n";
      for (int x = 0; x < nXB(r, q); ++x) {
        std::cout << "      xB [" << xBLow(r, q, x) << ", " << xBHigh(r, q, x) << "] : "
                  << nT(r, q, x) << " |t| bin(s)\n";
      }
    }
  }
}
