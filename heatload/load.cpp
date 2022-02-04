
#include <assert.h>
#include <string>

#include "adios2.h"

#define GET(X, i, j) X[i * 9 + j]

void load_data(adios2::IO reader_io, adios2::Engine reader,
               std::vector<long> &igid, std::vector<float> &iphase,
               std::vector<long> &egid, std::vector<float> &ephase) {
  // Clear vector
  igid.clear();
  egid.clear();
  iphase.clear();
  ephase.clear();

  // Inquire variables
  auto var_igid = reader_io.InquireVariable<long>("igid");
  auto var_egid = reader_io.InquireVariable<long>("egid");
  auto var_iphase = reader_io.InquireVariable<float>("iphase");
  auto var_ephase = reader_io.InquireVariable<float>("ephase");

  var_igid.SetSelection({{0}, {var_igid.Shape()[0]}});
  var_egid.SetSelection({{0}, {var_egid.Shape()[0]}});
  var_iphase.SetSelection(
      {{0, 0}, {var_iphase.Shape()[0], var_iphase.Shape()[1]}});
  var_ephase.SetSelection(
      {{0, 0}, {var_ephase.Shape()[0], var_ephase.Shape()[1]}});

  reader.Get<long>(var_igid, igid);
  reader.Get<long>(var_egid, egid);
  reader.Get<float>(var_iphase, iphase);
  reader.Get<float>(var_ephase, ephase);
}