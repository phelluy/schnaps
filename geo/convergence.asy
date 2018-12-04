import graph;
import markers;

// Reads the output of convergence_test.py and produces a pdf or eps file.

size(200, 150, IgnoreAspect);

scale(Log,Log);

string filename;
filename = getstring("external data");
file fin = input(filename).line();
real[][] a = fin.dimension(0,0);
a = transpose(a);
real[] deg = a[0];
real[] nraf = a[1];
real[] dof = a[2];
real[] error = a[3];
string xvals =  getstring("DOF, nraf, or dx");
real dx[];
if(a.length > 4) {
  dx = a[4];
} else {
  // Values for the CEMRACS14 test case (to be deprecated with new
  // python script where the dx values are saved).
  for(int i = 0; i < error.length; ++i) {
    real val;
    if(deg[i] == 1) 
      val = 0.083333 / nraf[i];
    if(deg[i] == 2) 
      val = 0.055556 / nraf[i];
    if(deg[i] == 3) 
      val = 0.041667 / nraf[i];
    if(deg[i] == 4) 
      val = 0.033333 / nraf[i];
    dx.push(val);
  }
}


// Find the dimensions in the file
real[] dimlist;
for(int i = 0; i < deg.length; ++i) {
  int d = round(deg[i]);
  bool found = false;
  for(int j = 0; j < dimlist.length; ++j) {
    if(d == dimlist[j]) {
      found = true;
      break;
    }
  }
  if(!found)
    dimlist.push(d);
}

// Loop over the dimensions
for(int d = 0; d < dimlist.length; ++d) { 
  //write(d);
  real dim = dimlist[d];
  real[] derror = {};
  real[] dnraf = {};
  real[] ddof = {};
  real[] ddx = {};
  for(int i = 0; i < deg.length; ++i) {
    if(deg[i] == dim) {
      derror.push(error[i]);
      dnraf.push(nraf[i]);
      ddof.push(dof[i]);
      ddx.push(dx[i]);
    }
  }
  //write(derror);

  real[] x;
  if(xvals == "DOF") x = ddof;
  if(xvals == "nraf") x = dnraf;
  if(xvals == "dx") x = ddx;
  
  int last = derror.length-1;
  if(last >= 1) {
    real order = log(derror[last-1] / derror[last]) / log(2.0);

    pen p = Pen(d);
    if(d == 1) p += dashed;
    if(d == 2) p += Dotted;
    if(d == 3) p += longdashdotted;
   
    draw(graph(x, derror), p, 
	 "Degree " + string(dim) + " convergence rate: " + string(order, 3),
	 MarkFill[0]);
  }
}

if(xvals == "DOF")
  xaxis("Degrees of Freedom", BottomTop, LeftTicks);
if(xvals == "nraf")
    xaxis("nraf", BottomTop, LeftTicks);
if(xvals == "dx")
    xaxis("$dx$", BottomTop, LeftTicks);
yaxis("L2 error", LeftRight, RightTicks);

attach(legend(),point(plain.E),20plain.E);
