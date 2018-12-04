
#define schnaps_real float   // !!!!!!!!!!!!!!  to be improved !!!!!!!!!

// user functions (defined in another source file and compiled at runtime)

//! \brief numflux function
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
extern void (*NumFlux)(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);

//! \brief boundary flux function
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
extern void (*BoundaryFlux)(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);

//! \brief source function
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] w :  state
//! \param[out] source : the source
extern void (*Source)(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);
  
//! \brief init data function
// !\param[in] x : space position
//! \param[out] w : init state at point x
extern void (*InitData)(schnaps_real *x, schnaps_real *w);

//! \brief imposed data function
//!\param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
extern void (*ImposedData)(const schnaps_real *x, const schnaps_real t, schnaps_real *w);

