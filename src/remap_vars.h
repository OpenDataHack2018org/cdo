#ifndef REMAP_VARS_H
#define REMAP_VARS_H

class gradientsType
{
  public:
  std::vector<double> grad_lat;
  std::vector<double> grad_lon;
  std::vector<double> grad_latlon;

  void init(size_t size)
    {
      grad_lat.resize(size);
      grad_lon.resize(size);
      grad_latlon.resize(size);
    }

  gradientsType(size_t size) { init(size); }
  gradientsType() {  }
};

enum struct RemapType
{
  UNDEF,
  BILINEAR,
  BICUBIC,
  DISTWGT,
  CONSERV,
  CONSERV_YAC
};

enum struct NormOpt
{
  NONE,
  DESTAREA,
  FRACAREA
};

typedef struct
{
  bool option;
  size_t max_links;
  size_t num_blks;
  size_t *num_links;
  size_t **src_add;
  size_t **dst_add;
  size_t **w_index;
} remaplink_t;

struct remapVarsType
{
  long links_per_value;
  bool sort_add;
  bool pinit;                       // true: if the pointers are initialized
  size_t max_links;                 // current size of link arrays
  size_t num_links;                 // actual number of links for remapping
  size_t num_wts;                   // num of weights used in remapping
  RemapType mapType;                // identifier for remapping method
  NormOpt normOpt;                  // option for normalization (conserv only)
  size_t resize_increment;          // default amount to increase array size

  std::vector<size_t> src_cell_add; // source grid address for each link
  std::vector<size_t> tgt_cell_add; // target grid address for each link
  std::vector<double> wts;          // map weights for each link [max_links*num_wts]

  remaplink_t links;
};


void remap(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array,
           gradientsType &gradients);
void remap_laf(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array);
void remap_sum(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array);
void remapVarsInit(RemapType mapType, size_t src_grid_size, size_t tgt_grid_size, remapVarsType &rv);
void remapVarsResize(remapVarsType &rv, int64_t increment);
void remapVarsReorder(remapVarsType &rv);
void remapVarsFree(remapVarsType &rv);
void remapVarsCheckWeights(const remapVarsType &rv);

#endif /* REMAP_VARS_H */
