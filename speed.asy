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
real[] dof = a[2];
real[] time = a[3];
real[] timeperrk = a[5];
real[] err = a[6];
string xvals =  getstring("DOF or nraf");
string yvals =  getstring("err, timeperrk, or time");

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
  real[] dtimeperrk = {};
  real[] dtime = {};
  real[] dnraf = {};
  real[] ddof = {};
  real[] derr = {};
  for(int i = 0; i < deg.length; ++i) {
    if(deg[i] == dim) {
      dtimeperrk.push(timeperrk[i]);
      dtime.push(time[i]);
      dnraf.push(nraf[i]);
      ddof.push(dof[i]);
      derr.push(err[i]);
    }
  }

  real[] x = xvals == "DOF" ? ddof : dnraf;
  real[] y;
  if(yvals == "timeperrk")
    y = dtimeperrk;
  if(yvals == "time")
    y = dtime;
  if(yvals == "err")
    y = derr;

  pen p = Pen(d);
  if(d == 1) p += dashed;
  if(d == 2) p += Dotted;
  if(d == 3) p += longdashdotted;
  
  int last = dtimeperrk.length-1;
  if(last >= 1) {
    draw(graph(x, y), p, "deg " + string(dim), MarkFill[0]);
  }
}

if(xvals == "DOF")
  xaxis("Degrees of Freedom", BottomTop, LeftTicks);
else
    xaxis("nraf", BottomTop, LeftTicks);
if(yvals == "timeperrk")
  yaxis("Computation time", LeftRight, RightTicks);
if(yvals == "time")
  yaxis("Computation time", LeftRight, RightTicks);
if(yvals == "err")
  yaxis("L2 error", LeftRight, RightTicks);

//attach(legend(),point(plain.E),20plain.E);
