import graph;
import markers;

// Reads the output of convergence_test.py and produces a pdf or eps file.

size(200, 150, IgnoreAspect);

scale(Log,Log);

string filename;
filename = getstring("external data");
file fin = input(filename).line();
real[][] a=fin.dimension(0,0);
a = transpose(a);
real[] deg = a[0];
real[] nraf = a[1];
real[] error = a[2];

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
  for(int i = 0; i < deg.length; ++i) {
    if(deg[i] == dim) {
      derror.push(error[i]);
      dnraf.push(nraf[i]);
    }
  }
  //write(derror);
  
  int last = derror.length-1;
  if(last >= 1) {
    real order = log(derror[last-1] / derror[last]) / log(2.0);
    draw(graph(dnraf, derror), Pen(d), 
	 "deg " + string(dim) + " convergence rate: " + string(order, 3),
	 MarkFill[0]);
  }
}

xaxis("nraf", BottomTop, LeftTicks);
yaxis("error", LeftRight, RightTicks);

attach(legend(),point(plain.E),20plain.E);
