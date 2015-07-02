#ifdef PARALUTION

  void paralution_fortran_solve_coo(int, int, int, char*, char*, char*, 
    char*, int*, int*, double*, double*, double, 
    double, double, int, int, int, int, double*);

//void paralution_fortran_solve_coo(int, int, int, char*, char*, char*, 
//  char*, const int*, const int*, const double*, const double*, double, 
//  double, double, int, int, int, int, double*);

  void paralution_fortran_solve_csr(int, int, int, char*, char*, char*, 
    char*, int*, int*, double*, double*, double, 
    double, double, int, int, int, int, double*, int*, double*, int*);

#endif /* PARALUTION */
