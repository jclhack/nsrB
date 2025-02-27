(*nsr_Bdd_20250226.mm: version for UIUC
  nsr_Bdd_20210301.mm: better match to proposed sample geometry
  nsr_Bdd_20210223.mm: classical field from polarized object with distance; compare field from magnetized 
  cylinder
  nsr_Bdd_20210222.mm: classical field from polarized object with distance
  fourier_ptb_Beff_20200701.mm: ptb version, more accurate geometry
  fourier_ptb_Beff_20200630.mm: ptb version
  fourier_ariadne_Beff_20200527.mm: better plots and scaling
  fourier_ariadne_Beff_20200525.mm: plot ariadne effective field and calculate amplitudes.
  fourier_ariadne_signal_20200523.mm: plot ariadne torque and calculate amplitudes.
  fourier_force_probe_20200422.mm: calculate fourier amplitude of force at probe, clean up code.
  fourier_stats_v1_20200422_2.mm: replace obsolete ErrorBarPlots with About function
  fourier_stats_v1_20200422.mm: first working version on Carbonate (path specification in OpenRead cmd)
  fourier_stats.mm: reads virtual work and phase data from MC, plots, computes Fourier amplitudes different ways
*)

Clear[mu0,g1,s1,n1,mtm,q,xs,ys,zs,vs,k,fstring,ntext,nmark,pntext,TSSP,VRawDat,Tclose,dimf,pRawDat,xmRawDat,ymRawDat,zmRawDat,xeRawDat,yeRawDat,zeRawDat,pmin,pmx,bxmin,bxmax,xRaw,xRawPlt,bymin,bymax,yRaw,yRawPlt,bzmin,bzmax,zRaw,zRawPlt,zRawPltl,prad,xmRawDatm,ymRawDatm,zmRawDatm,ampxmm,sampxmc,camxmc,ampxmc,xd,yd,zd,vd,N1,mmtm,zmfa,zmfaPlt,zmfaPltl,M,zmna,zmnaPlt,zmnaPltl];

(*Field parameters*)
mu0  = 4. Pi 10^-7;  (*perm. free space, N/(A^2)*)
g1   = 2.8025 10^10; (*e gyro ratio, 1/(sT)*)
s1   = (1.05 10^-34)/2.0;     (*e spin, Js*)
n1   = 4.0 10^26;    (*sample (e) spin density, 1/m^3*)
mtm  = 1.0 10^-6;    (*microns to m*)
q    = mtm^3 n1 mu0 g1 s1/(4.0 Pi); (*MC to field*volume conversion factor, gives field*volume in Tm^3*)

(*source geometry*)
xs = 0.1 1.68826 10^-2; (*source x-size, m*) 
ys = 0.1 1.68826 10^-2; (*source y-size, m*) 
zs = 0.1 1.0 10^-2;   (*source z-size, m*) 
vs = xs ys zs;    (*source volume, m^3*)
k = q/vs;         (*MC to field conversion factor, gives field in T*)

(*Enter virtual work data file as strings*)
fstring = "W_phase.txt"

Print["Virtual work data file:"]
Print[fstring]

(*Open data files for reading*)
ntext = OpenRead["C:/Users/jcl/mrg/mm/"<>fstring]

(*Find last line in header*)
nmark = Find[ntext, "Bzdev"]

(*Find position of last byte in this line*)
pntext = StreamPosition[ntext]

(*Set pointer to this position*)
TSSP = SetStreamPosition[ntext, pntext]

(*Read rest of file as list of ordered number pairs*)
VRawDat = ReadList[ntext, Number, RecordLists -> True];

(*Close files*)
Tclose = Close[ntext];

(*Get number of phase points*)
dimf = Dimensions[VRawDat][[1]]

(*Separate X and Y data*)
pRawDat  = 1.0 10^-6 VRawDat[[All,1]]
xmRawDat = k VRawDat[[All,2]]
ymRawDat = k VRawDat[[All,3]]
zmRawDat = k VRawDat[[All,4]]
xeRawDat = Abs[k] VRawDat[[All,5]]
yeRawDat = Abs[k] VRawDat[[All,6]]
zeRawDat = Abs[k] VRawDat[[All,7]]

(*Plot*)
pmin=0.0 Min[pRawDat];
pmax=1.05 Max[pRawDat];
bxmin=Min[xmRawDat]-Max[xeRawDat];
bxmax=Max[xmRawDat]+Max[xeRawDat];
xRaw    = Table[{pRawDat[[i]],Around[xmRawDat[[i]],xeRawDat[[i]]]},{i,1,dimf}];
xRawPlt = ListPlot[xRaw, PlotStyle -> {PointSize[.020]},
         PlotRange -> {{pmin,pmax},{bxmin,bxmax}}, PlotLabel -> "Bx (transverse) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_x (T)"},Frame->True];

bymin=Min[ymRawDat]-Max[yeRawDat];
bymax=Max[ymRawDat]+Max[yeRawDat];
yRaw    = Table[{pRawDat[[i]],Around[ymRawDat[[i]],yeRawDat[[i]]]},{i,1,dimf}];
yRawPlt = ListPlot[yRaw, PlotStyle -> {PointSize[.020]},
         PlotRange -> {{pmin,pmax},{bymin,bymax}}, PlotLabel -> "By (transverse) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_y (T)"}, Frame->True];

bzmin=Min[zmRawDat]-Max[zeRawDat];
bzmax=1.5 (Max[zmRawDat]+Max[zeRawDat]);
zRaw    = Table[{pRawDat[[i]],Around[zmRawDat[[i]],zeRawDat[[i]]]},{i,1,dimf}];
zRawPlt = ListPlot[zRaw, PlotStyle -> {PointSize[.020]},
         PlotRange -> {{pmin,pmax},{bzmin,bzmax}}, PlotLabel -> "Bz (longitudinal) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_z (T)"}, Frame->True];

zRawPltl = ListLogPlot[zRaw, PlotStyle -> {PointSize[.020]},
         PlotRange -> {{pmin,pmax},{bzmin,bzmax}}, PlotLabel -> "Bz (longitudinal) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_z (T)"}, Frame->True];


(*Show*)
Show[xRawPlt, DisplayFunction -> (Print[#];&)]
Print[" "]
Print[" "]
Show[yRawPlt, DisplayFunction -> (Print[#];&)]
Print[" "]
Print[" "]
Show[zRawPlt, DisplayFunction -> (Print[#];&)]

(*Analytical*)
(*sample geometry*)
xd = 1.68826 10^-2; (*sample x-size, m*) 
yd = 1.68826 10^-2; (*sample y-size, m*) 
zd = 1.0 10^-2;   (*sample z-size, m*) 
vd = xd yd zd;    (*sample volume, m^3*)

(*field from cylinder on axis [Jackson 5-19]*)
M   = n1 g1 s1; (*magnetization*)
ae  = xd/Sqrt[Pi]; (*effective sample radius, m*)
Bza = -(mu0 M/2)((z-zd/2)/((ae^2+(z-zd/2)^2)^(1/2))-(z+zd/2)/((ae^2+(z+zd/2)^2)^(1/2))); (*field, T*)

BzaPlt = Plot[Bza, {z,pmin,pmax}, PlotStyle->{Green}, PlotRange -> {{pmin,pmax},{bzmin,bzmax}}, PlotLabel -> "Bz analytical (along source travel) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_z (T)"}, Frame->True];

BzaPltl = LogPlot[Bza, {z,pmin,pmax}, PlotStyle->{Green}, PlotRange -> {{pmin,pmax},{bzmin,bzmax}}, PlotLabel -> "Bz analytical (along source travel) vs. position", GridLines -> Automatic, FrameLabel->{"z (m)", "B_z (T)"}, Frame->True];


Show[BzaPlt, DisplayFunction -> (Print[#];&)]

Show[BzaPlt, zRawPlt, DisplayFunction -> (Print[#];&)]

Show[BzaPltl, zRawPltl, DisplayFunction -> (Print[#];&)]



