// new generation real-time PS code based on kp_rtACQ_39.c

/*
_01	fisrt implementation
_02	npoints, droppts1, droppts2 has the same meaning as np (they all must be even numbers, forced by setlimit)
_04	new version using kp_cycles and kp_npoints for setting up the chunks and loops
	tested, works fine; droppts1 and droppts2 needs processing macro kp_rt2D_proc_v3
_05	adding flags BIRD, BIRDmode, and making compatibility with kp_makePS7 (bw_a; pw180_a; etc..)
_08	fine; half completed with options
_09	working on ifzero PFGs; fine
_10	kp_scyc new implementation
_11	tau_p and tau_r
_12	altJ='y'/'n'/'k' to alternate the order of hard180 and J-refocusing element in the loop 
	problems; use altJ='n'
_13	works well
_15	add options to change pulses only in real-time loop for MRC
	altJ='k' has been enforced
_16	optional ACQ_tpwr, etc.. for MRC paper
*/


/*----------------------
Developed By NMR Group
School of Chemistry
University of Manchester
United Kingdom
Apr 2017
----------------------*/


#include <standard.h>

static int	ph1[32] = {1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3},		//v1
		ph2[32] = {0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2},		//v2
		ph3[32] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},		//v3
		ph4[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},		//v4
		ph5[32] = {1,3,1,3,3,1,3,1,3,1,3,1,1,3,1,3,3,1,3,1,1,3,1,3,1,3,1,3,3,1,3,1};		//oph			  			  
static int	ph7[32] = {0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},		//v7 - 1st 90 of bird and the hard 180
		ph8[32] = {1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2},		//v8 - simpulse 180 of bird
		ph9[32] = {2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3};		//v9 - 2nd 90 of bird 

static int	nph1[32] = {1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,3,1,1,3,3},		//v1
		nph2[32] = {0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2},		//v2
		nph3[32] = {0,0,0,0,2,2,2,2,0,0,0,0,2,2,2,2,0,0,0,0,2,2,2,2,0,0,0,0,2,2,2,2},		//v3
		nph4[32] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},		//v4
		nph5[32] = {1,3,3,1,3,1,1,3,3,1,1,3,1,3,3,1,1,3,3,1,3,1,1,3,3,1,1,3,1,3,3,1};		//oph 			  
static int	nph7[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},		//v7 - 1st 90 of bird and the hard 180
		nph8[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},		//v8 - simpulse 180 of bird
		nph9[32] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};		//v9 - 2nd 90 of bird 



//v10;v11 reserved for the real-time loop
//v12;v13;v14 reserved for ifzero PFGs


//v6 for auto-incremented supercycle additional to the scan-to-scan 
//v16,v18 used for small-angle phases
static int	 ph6m4[ 4] = {0,0,2,2},
		 ph6t5[ 5] = {0,5,2,5,0},
		 ph6t4[20] = {0,5,2,5,0,0,5,2,5,0,6,11,8,11,6,6,11,8,11,6},
		ph6m16[16] = {0,0,2,2, 0,2,2,0, 2,2,0,0, 2,0,0,2},
		 ph6t7[ 7] = {0,7,20,17,20,7,0},
		 ph6t9[ 9] = {0,1,12,11,18,11,12,1,0},
		   ph6[ 4] = {0,0,0,0};	


pulsesequence()
{

//HSQC part
double  evolcorr=2.0*pw+4.0e-6,	
	j1xh = getval("j1xh"),
	tau = (double)(1.0e-6*(floor)(1.0e6/(4.0*j1xh))), 	// delay as integer multiple of 1.0us
	taug=2.0*tau,

	ACQ_j1xh = getval("ACQ_j1xh"),
	ACQ_tau = (double)(1.0e-6*(floor)(1.0e6/(4.0*ACQ_j1xh))), 	// delay as integer multiple of 1.0us

	multh = getval("multh");

int	phase1 = (int)(getval("phase")+0.5),
	ZZgsign=1.0, 
	icosel;
//BIRD
double	//d1_corr=0.0,
	rof3=getval("rof3"), 		//delay for receiver off - can be zero if ddrpm='r'

	pw_HBIP = getval("pw_HBIP"),
	pw_XBIP = getval("pw_XBIP"), 
	pwr_HBIP = getval("pwr_HBIP"),	
	pwr_XBIP = getval("pwr_XBIP"), 	

	//tau_a = getval("tau_a"),
	tau_p = getval("tau_p"),
	tau_r=0.0, 	// calculated by the sequence from tau_p
	
	tpwr = getval("tpwr"),
	tpwrf = getval("tpwrf"),
	pw = getval("pw"),

	ACQ_tpwr = getval("ACQ_tpwr"),
	ACQ_pw = getval("ACQ_pw"),
	ACQ_pwxlvl = getval("ACQ_pwxlvl"),
	ACQ_pwx = getval("ACQ_pwx"),

//hard180 refocusing
	pp = getval("pp"),
	pplvl = getval("pplvl"),
	pplvlf = getval("pplvlf"),
//sel180 refocusing of the active spin
	pw180_a = getval("pw180_a"),
	pwr180_a = getval("pwr180_a"),


	dmf = getval("dmf"),
	dres = getval("dres"),

	ACQ_gt1 = getval("ACQ_gt1"),  
	ACQ_gzlvl1 = getval("ACQ_gzlvl1"),
	ACQ_gzlvl2 = getval("ACQ_gzlvl2"),
	ACQ_gzlvl3 = getval("ACQ_gzlvl3"),
	ACQ_gzlvl4 = getval("ACQ_gzlvl4"),
	ACQ_gstab = getval("ACQ_gstab"),
	

	offRes_DEC = getval("offRes_DEC"), 	//not used if zero; set real dof for HPWR DEC
	offRes_X180 = getval("offRes_X180"), 	//not used if zero; set real dof for interchunk pulses
	
	calH = getval("calH"), 		//calibrate pw90
	calX = getval("calX"), 
	cal_pwlvl = getval("cal_pwlvl"), 
	cal_pwxlvl = getval("cal_pwxlvl"), 
	cal_pw = getval("cal_pw"), 	//if cal_pwlvl!=tpwr then use this with calH
	cal_pwx = getval("cal_pwx"), 	//if cal_pwxlvl!=pwxlvl then use this with calX


	//gt1 = getval("gt1"),		//CTP during ACQ for BIRD or NHsel
	//gzlvl1 = getval("gzlvl1"), 	
	//gt2 = getval("gt2"),		//CTP during ACQ for hard 180
	//gzlvl2 = getval("gzlvl2"), 	
	gstab = getval("gstab"),


	npoints = getval("kp_npoints"),	//number of lines (np/2) of each pure shift data points
	droppts1 = getval("droppts1"),	//droppts at the begining of each chunk except the 1st chunk
	droppts2 = getval("droppts2"),	//droppts at the end of each chunk
	kp_cycles =(double) (floor)( getval("kp_cycles") );

	F_initval(kp_cycles,v10);
int	kpph = (int)(getval("kpph"));

char   	shp_HBIP[MAXSTR], shp_XBIP[MAXSTR], 
	lkgate_flg[MAXSTR], 
	homo[MAXSTR],
	kp_scyc[MAXSTR],	// 'n'/'y' 
	dseq[MAXSTR], dm[MAXSTR],
	shp_a[MAXSTR], 
	hsgradaxis[MAXSTR],	//gradient axis for homospoil
	gradaxis[MAXSTR];	//gradient axis during ACQ ("t" implement simultanous XY)

	getstr("shp_HBIP",shp_HBIP);
	getstr("shp_XBIP",shp_XBIP);
	getstr("shp_a",shp_a);		//created by macro kp_makePS8 using bw_a/kp_phincr_a
        getstr("lkgate_flg",lkgate_flg);
	getstr("homo",homo);
	getstr("kp_scyc",kp_scyc);
	getstr("dseq",dseq);
	getstr("dm",dm);
        getstr("hsgradaxis",hsgradaxis);		
        getstr("gradaxis",gradaxis);	

//extensions for AD
//double	pwx180 = getval("pwx180"),
//	pwxlvl180 = getval("pwxlvl180"),
//	pwx180r = getval("pwx180r"),
//	pwxlvl180r = getval("pwxlvl180r");
//char	pwx180ad[MAXSTR],
//	pwx180adR[MAXSTR];
//	getstr("pwx180ad", pwx180ad);
//	getstr("pwx180adR", pwx180adR);

//gradients
double 	gtE = getval("gtE"),		//HSQC encoding
	gzlvlE = getval("gzlvlE"),
	gtD = getval("gtD"),		//HSQC decoding
	gzlvlD = getval("gzlvlD"),

	hsglvl = getval("hsglvl"),
        hsgt = getval("hsgt"),
        hsgstab = getval("hsgstab");	

char   	BIRD[MAXSTR],          	 // Flag to choose gHSQC/rtgHSQC-BIRD ('n'/'y')
	kp_hsqc[MAXSTR],
	BIRDmode[MAXSTR]; 	//Flag to choose hard/bip/wurst2i ('h'/'b'/'w') 13C inversion pulse within BIRD and NHselective refocusing via 'n'

	getstr("BIRD",BIRD);
getstr("kp_hsqc",kp_hsqc);
	getstr("BIRDmode",BIRDmode);

char    sspul[MAXSTR], nullflg[MAXSTR], PFGflg[MAXSTR];

  	getstr("sspul",sspul);
	getstr("nullflg",nullflg);
  	getstr("PFGflg",PFGflg);


//d1 correction was normalized to BIRDmode='b' and only supports HSQC/bip-BIRD/NHsel
	
//	if (BIRD[0]=='n') d1_corr=(np/npoints)*(pw*4.0+1.0/j1xh+pw_XBIP+tauA+tauB+tauC+2.0*(gt1+gt2+2.0*gstab)+8.0*rof1);
	//if ((BIRD[0]=='n') && (BIRDmode[0]=='b')) d1_corr=(np/npoints)*(pw*4.0+1.0/j1xh+pw_XBIP+tauA+tauB+tauC+2.0*(gt1+gt2+2.0*gstab)+8.0*rof1);
	//if ((BIRD[0]=='n') && (BIRDmode[0]=='n')) d1_corr=(np/npoints)*(pw180_a+2.0*tauA+(pw*2.0+tauA+tauB+tauC)+2.0*(gt1+gt2+2.0*gstab)+4.0*rof1);	
	
//	if ((BIRD[0]=='y') && (BIRDmode[0]=='n')) d1_corr=(np/npoints)*(4.0*rof1+(pw*4.0+1.0/j1xh+pw_XBIP+tauA+tauB+tauC+2.0*(gt1+gt2+2.0*gstab))-(pw180_a+2.0*tauA+(pw*2.0+tauA+tauB+tauC)+2.0*(gt1+gt2+2.0*gstab)));	
//	if ((BIRD[0]=='y') && (BIRDmode[0]=='b')) d1_corr=0.0;
	



if (BIRDmode[0]=='n')
{
	tau_r = (tau_p -4.0*ACQ_gt1-4.0*ACQ_gstab- -pp-4.0*rof1 -2.0*rof1-pw180_a ) /2.0;
}
if (BIRDmode[0]=='h')
{
	if (2.0*ACQ_pwx>2.0*ACQ_pw) 
	{
	tau_r = (tau_p -4.0*ACQ_gt1-4.0*ACQ_gstab- -pp-4.0*rof1 -2.0*rof1-6.0*ACQ_pwx-2.0*ACQ_pw -4.0*tau) /2.0;
	}
	else
	{
	tau_r = (tau_p -4.0*ACQ_gt1-4.0*ACQ_gstab- -pp-4.0*rof1 -2.0*rof1-4.0*ACQ_pwx-4.0*ACQ_pw  -4.0*tau) /2.0;
	}
}
if (BIRDmode[0]=='b')
{
	if (pw_XBIP>pw_HBIP) 
	{
	tau_r = (tau_p -4.0*ACQ_gt1-4.0*ACQ_gstab- -pp-4.0*rof1 -2.0*rof1-3.0*pw_XBIP-2.0*ACQ_pw -4.0*tau) /2.0;
	}
	else
	{
	tau_r = (tau_p -4.0*ACQ_gt1-4.0*ACQ_gstab- -pp-4.0*rof1 -2.0*rof1-2.0*pw_XBIP-pw_HBIP-2.0*ACQ_pw -4.0*tau) /2.0;
	}
}

	
/*	if (d1_corr > 0.2*d1)
	{
	abort_message("Correction for relaxation delay is more than 20 percent! Check ACQ parameters...\n");
	}
*/

	
//evolcorr and mult declarations
	evolcorr=2.0*pw+4.0e-6;
	  
	if (multh > 0.5)
   	taug = 2.0*tau;
   	else
    	taug = gtE + gstab + 2*GRADIENT_DELAY;
   	ZZgsign=-1;
   	if (multh == 2) ZZgsign=1;
   	icosel = 1;

//rtgHSQC-BIRD phases
if (kpph>0)
{
	settable(t1,kpph,ph1);
	settable(t2,kpph,ph2);
	settable(t3,kpph,ph3);
	settable(t4,kpph,ph4);
	settable(t5,kpph,ph5);
	settable(t7,kpph,ph7);
	settable(t8,kpph,ph8);
	settable(t9,kpph,ph9);
}
if (kpph==0)
{
	settable(t1,32,ph1);
	settable(t2,32,ph2);
	settable(t3,32,ph3);
	settable(t4,32,ph4);
	settable(t5,32,ph5);
	settable(t7,32,ph7);
	settable(t8,32,ph8);
	settable(t9,32,ph9);
}
if (kpph<0)
{
	settable(t1,-kpph,nph1);
	settable(t2,-kpph,nph2);
	settable(t3,-kpph,nph3);
	settable(t4,-kpph,nph4);
	settable(t5,-kpph,nph5);
	settable(t7,-kpph,nph7);
	settable(t8,-kpph,nph8);
	settable(t9,-kpph,nph9);
}

	
	getelem(t7, ct, v7);
  	getelem(t8, ct, v8);
	getelem(t9, ct, v9);

  	getelem(t1, ct, v1);		
  	getelem(t2, ct, v2);
  	getelem(t3, ct, v3);
  	getelem(t4, ct, v4);
	getelem(t5, ct, oph);
  	

   	initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v5);
   	if ((phase1 == 2) || (phase1 == 5))
     	icosel = -1;

   	add(v2,v5,v2);	
   	add(oph,v5,oph);


//no phase sequencing
if (kp_scyc[0]=='n')
{
	settable(t6,4,ph6);
}
//m4 phase sequencing
if ((kp_scyc[0]=='m') && (kp_scyc[1]=='4'))
{
settable(t6,4,ph6m4);
}
//m16 phase sequencing
if ((kp_scyc[0]=='m') && (kp_scyc[1]=='1') && (kp_scyc[2]=='6'))
{
settable(t6,16,ph6m16);
}
//t5 phase sequencing
if ((kp_scyc[0]=='t') && (kp_scyc[1]=='5'))
{
settable(t6,5,ph6t5);
obsstepsize(30.0);
decstepsize(30.0);
}
//t5m4 phase sequencing
if ((kp_scyc[0]=='t') && (kp_scyc[1]=='4'))
{
settable(t6,20,ph6t4);
obsstepsize(30.0);
decstepsize(30.0);
}
//t7 phase sequencing
if ((kp_scyc[0]=='t') && (kp_scyc[1]=='7'))
{
settable(t6,7,ph6t7);
obsstepsize(15.0);
decstepsize(15.0);
}
//t9 phase sequencing
if ((kp_scyc[0]=='t') && (kp_scyc[1]=='9'))
{
settable(t6,9,ph6t9);
obsstepsize(15.0);
decstepsize(15.0);
}
assign(zero,v16);
getelem(t6,v16,v16);

	if (kp_scyc[0]=='t')
	{
		if ((kp_scyc[1]=='4') || (kp_scyc[1]=='5') )
		{
		mult(three,v8,v18);
		add(v18,v16,v16);		
		}
		if ((kp_scyc[1]=='7') || (kp_scyc[1]=='9') )
		{
		mult(three,v8,v18);
		mult(two,v18,v18);
		add(v18,v16,v16);
		}
	}


//ifzero PFGs
sub(v10,one,v14);

/* BEGIN PULSE SEQUENCE */
status(A);

//correct the total experiment time for HSQC-ref experiment  
//	if ((d1_corr>0.0) && (BIRD[0]=='n')) delay(d1_corr);
//	if ((d1_corr>0.0) && (BIRD[0]=='y') && (BIRDmode[0]=='n')) delay(d1_corr);

	if (sspul[A] == 'y')
	{
//lock gating
if (lkgate_flg[0] == 'y')  lk_hold(); /* turn lock sampling off */		
  
      if (PFGflg[A] == 'y')
      {
	obspower(tpwr);
	obspwrf(tpwrf);
	
	delay(5.0e-5);
	if (hsgt>0.0)
	{
	  if (hsgradaxis[0]=='t')  
	  { 
	    rgradient('x',hsglvl);
	    rgradient('y',hsglvl); 
	    delay(hsgt); 
	    rgradient('x',0.0); 
	    rgradient('y',0.0); 
	  }   
	  else 
	  {
	    if ((hsgradaxis[0]=='x') || (hsgradaxis[0]=='y')) 
	    { 
	    rgradient(hsgradaxis[0],hsglvl); 
	    delay(hsgt); 
	    rgradient(hsgradaxis[0],0.0); 
	    } 
	    else zgradpulse(hsglvl,hsgt); 
	  }
	}
        rgpulse(pw,zero,rof1,rof1);
	if (hsgt>0.0)
	{
	  if (hsgradaxis[0]=='t')  
	  { 
	    rgradient('x',hsglvl);
	    rgradient('y',hsglvl); 
	    delay(hsgt); 
	    rgradient('x',0.0); 
	    rgradient('y',0.0); 
	  }   
	  else 
	  {
	    if ((hsgradaxis[0]=='x') || (hsgradaxis[0]=='y')) 
	    { 
	    rgradient(hsgradaxis[0],hsglvl); 
	    delay(hsgt); 
	    rgradient(hsgradaxis[0],0.0); 
	    } 
	    else zgradpulse(hsglvl,hsgt); 
	  }
	}

      }
       	else
       	{
        obspower(tpwr-12);
	delay(5.0e-5);
        rgpulse(500*pw,zero,rof1,rof1);
        rgpulse(500*pw,one,rof1,rof1);
	}
	}

	if (nullflg[0]=='y') obspower(tpwr);
	else obspower(cal_pwlvl);
	decpower(pwxlvl);
	txphase(zero);
        decphase(zero);
	obsoffset(tof);
	decoffset(dof);
	
if (lkgate_flg[0] == 'y')  lk_sample(); /* turn lock sampling on */	
	delay(d1);
if (lkgate_flg[0] == 'y')  lk_hold(); /* turn lock sampling off */		
	delay(5.0e-5);

status(B);  

if (kp_hsqc[0]=='y')
{

/****** null flag starts here *****/

      if (getflag("nullflg"))
      {
        rgpulse(0.5*pw,zero,rof1,rof1);
	txphase(zero);
	delay(2.0*tau);
	simpulse(2.0*pw,2.0*pwx,zero,zero,rof1,rof1); 
	txphase(two);
        delay(2.0*tau);
        rgpulse(1.5*pw,two,rof1,rof1);
	txphase(zero);
	
/* purgeing gradient */	
	if (hsgt>0.0)
	{
	  if (hsgradaxis[0]=='t')  
	  { 
	    rgradient('x',hsglvl);
	    rgradient('y',hsglvl); 
	    delay(hsgt); 
	    rgradient('x',0.0); 
	    rgradient('y',0.0); 
	  }   
	  else 
	  {
	    if ((hsgradaxis[0]=='x') || (hsgradaxis[0]=='y')) 
	    { 
	    rgradient(hsgradaxis[0],hsglvl); 
	    delay(hsgt); 
	    rgradient(hsgradaxis[0],0.0); 
	    } 
	    else zgradpulse(hsglvl,hsgt); 
	  }
	}
	if (cal_pwlvl!=tpwr) obspower(cal_pwlvl);
	delay(hsgstab);
      }

/****************************gHSQC or gHSQC part of pure shift starts here *************************/

//option to proton calibrations	  

	rgpulse(cal_pw*calH,zero,rof1,rof1); 	
	if (cal_pwlvl!=tpwr) obspower(tpwr);
	delay(tau);
	simpulse(2.0*pw,2.0*pwx,zero,zero,rof1,rof1); 
	txphase(v1);
	delay(tau);
     	rgpulse(pw,v1,rof1,rof1);
	
/* purgeing gradient */	
	if (hsgt>0.0)
	{
	  if (hsgradaxis[0]=='t')  
	  { 
	    rgradient('x',hsglvl);
	    rgradient('y',hsglvl); 
	    delay(2.0*hsgt); 
	    rgradient('x',0.0); 
	    rgradient('y',0.0); 
	  }   
	  else 
	  {
	    if ((hsgradaxis[0]=='x') || (hsgradaxis[0]=='y')) 
	    { 
	    rgradient(hsgradaxis[0],hsglvl); 
	    delay(2.0*hsgt); 
	    rgradient(hsgradaxis[0],0.0); 
	    } 
	    else zgradpulse(hsglvl,2.0*hsgt); 
	  }
	}
   	decphase(v2);
	if (cal_pwxlvl!=pwxlvl) decpower(cal_pwxlvl);
	delay(hsgstab);
//option to pwx calibrations	
	decrgpulse(cal_pwx*calX, v2, rof1, 2.0e-6); 
	txphase(zero);
	decphase(zero);
	if (cal_pwxlvl!=pwxlvl) decpower(pwxlvl);

	delay(d2/2.0);				// First half of t1 evolution
	rgpulse(2.0*pw,zero,2.0e-6,2.0e-6);
	delay(d2/2.0);				// Second half of t1 evolution

     	zgradpulse(gzlvlE,gtE);
	delay(taug - gtE);
	simpulse(multh*pw,2.0*pwx,zero,zero,rof1,rof1);
	delay(taug + evolcorr); 

	decrgpulse(pwx,v4,2.0e-6,rof1);
	
/* purgeing gradient */	
	if (hsgt>0.0)
	{
	  if (hsgradaxis[0]=='t')  
	  { 
	    rgradient('x',ZZgsign*0.6*hsglvl);
	    rgradient('y',ZZgsign*0.6*hsglvl); 
	    delay(1.2*hsgt); 
	    rgradient('x',0.0); 
	    rgradient('y',0.0); 
	  }   
	  else 
	  {
	    if ((hsgradaxis[0]=='x') || (hsgradaxis[0]=='y')) 
	    { 
	    rgradient(hsgradaxis[0],ZZgsign*0.6*hsglvl); 
	    delay(1.2*hsgt); 
	    rgradient(hsgradaxis[0],0.0); 
	    } 
	    else zgradpulse(ZZgsign*0.6*hsglvl,1.2*hsgt); 
	  }
	}
	txphase(v3);
	delay(hsgstab);
/**/	
     	rgpulse(pw,v3,rof1,rof1);
	delay(tau - (2.0*pw/PI) - 2.0*rof1);

	simpulse(2.0*pw,2.0*pwx,zero,zero,rof1,rof1); 

	zgradpulse(icosel*gzlvlD,gtD);
	decpower(dpwr);	

	if (offRes_DEC!=0.0) decoffset(dof+offRes_DEC);

	if (BIRD[0]=='y')
	{
	delay(tau - gtD - 50.0e-9 - droppts1/sw/2.0);
	}
	else
	{
	delay(tau - gtD - 50.0e-9);
	}

}
else
{
decpower(dpwr);
rgpulse(pw,oph,rof1,rof1);

	if (pwr180_a!=tpwr)  obspower(pwr180_a);
zgradpulse(gzlvlD,gtD);
delay(0.001);
delay(droppts1/sw/2.0);
	shaped_pulse(shp_a,pw180_a,zero,rof1,rof1);
	if (pwr180_a!=tpwr) obspower(tpwr);
zgradpulse(gzlvlD,gtD);
delay(0.001);
}	



/********************************gHSQC part stops and BIRD Acquisition starts here**************************/

	//delay(1.0/(getval("fb")*1.3))
							
if (BIRD[0]=='y')
{
setacqmode(WACQ|NZ);	

//setup phase sequencing
	if (kp_scyc[0]=='t')
	{
	dcplrphase(v16);
	}
	else
	{
	decphase(v16);
	}

	obsblank();	
	delay(rof2);
	startacq(alfa);
	acquire(droppts1,1.0/sw);
	if (dm[2]=='y')
	{
	decprgon(dseq, 1.0/dmf, dres);
	decon(); decunblank();
	}
//	acquire(npoints/2.0+(droppts1+droppts2),1.0/sw);
	acquire(npoints/2.0,1.0/sw);
//	recoff();					 	
	if (dm[2]=='y')
	{
	decoff();
	decprgoff();
	}
	acquire(droppts2,1.0/sw);
	recoff();					 	

loop(v10,v11);


	obspower(pplvl);	
	if (pplvlf!=tpwrf) obspwrf(pplvlf);
	obsunblank();
	//change offset
	if (( (offRes_X180!=0.0) || (offRes_DEC!=0.0)) && (offRes_X180!=offRes_DEC)) decoffset(dof+offRes_X180);

//setup phase sequencing
	getelem(t6,v11,v16);
	if (kp_scyc[0]=='t')
	{
	xmtrphase(v16);
	dcplrphase(v16);
	}
	else
	{
	add(v16,v7,v7);
	add(v16,v8,v8);
	add(v16,v9,v9);
	}
	mod4(v11,v12); //PFG levels

	delay(rof1-50.0e-9);

//apply the next PFG level
if (ACQ_gt1>0.0)
{

	if (gradaxis[0]=='t')  
	{ 
		ifzero(v12);
		 rgradient('x',ACQ_gzlvl1);
		 rgradient('y',ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient('x',ACQ_gzlvl2);
		 rgradient('y',ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient('x',ACQ_gzlvl3);
		 rgradient('y',ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient('x',ACQ_gzlvl4);
		 rgradient('y',ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
	}   
	else 
	{    
    		if ((gradaxis[0]=='x') || (gradaxis[0]=='y')) 
		{ 
		ifzero(v12);
		 rgradient(gradaxis[0],ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		} 
		else 
		{
		ifzero(v12);
		 zgradpulse(ACQ_gzlvl1,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 zgradpulse(ACQ_gzlvl2,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 zgradpulse(ACQ_gzlvl3,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 zgradpulse(ACQ_gzlvl4,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		}
	}   
}	//endif ACQ_gt1>0.0

	rgpulse(pp,v7,rof1,rof1);

	delay(rof1-50.0e-9);

//apply the next PFG level
if (ACQ_gt1>0.0)
{

	if (gradaxis[0]=='t')  
	{ 
		ifzero(v12);
		 rgradient('x',ACQ_gzlvl1);
		 rgradient('y',ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient('x',ACQ_gzlvl2);
		 rgradient('y',ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient('x',ACQ_gzlvl3);
		 rgradient('y',ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient('x',ACQ_gzlvl4);
		 rgradient('y',ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
	}   
	else 
	{    
    		if ((gradaxis[0]=='x') || (gradaxis[0]=='y')) 
		{ 
		ifzero(v12);
		 rgradient(gradaxis[0],ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		} 
		else 
		{
		ifzero(v12);
		 zgradpulse(ACQ_gzlvl1,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 zgradpulse(ACQ_gzlvl2,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 zgradpulse(ACQ_gzlvl3,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 zgradpulse(ACQ_gzlvl4,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		}
	}   
}	//endif ACQ_gt1>0.0


//droppts1 and droppts2
delay((droppts1+droppts2)/sw/2.0);

	if (BIRDmode[0]== 'n')
	{
	obspower(pwr180_a);	
	if (tpwrf!=4095.0) obspwrf(4095.0);
	}
	else
	{
	obspower(ACQ_tpwr);	
	if (tpwrf!=4095.0) obspwrf(tpwrf);
	}

	
	delay(tau_r);


//apply the next PFG level
if (ACQ_gt1>0.0)
{
	if (gradaxis[0]=='t')  
	{ 
		ifzero(v12);
		 rgradient('x',ACQ_gzlvl3);
		 rgradient('y',ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient('x',ACQ_gzlvl4);
		 rgradient('y',ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient('x',ACQ_gzlvl1);
		 rgradient('y',ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient('x',ACQ_gzlvl2);
		 rgradient('y',ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
	}   
	else 
	{    
    		if ((gradaxis[0]=='x') || (gradaxis[0]=='y')) 
		{ 
		ifzero(v12);
		 rgradient(gradaxis[0],ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		} 
		else 
		{
		ifzero(v12);
		 zgradpulse(ACQ_gzlvl3,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 zgradpulse(ACQ_gzlvl4,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 zgradpulse(ACQ_gzlvl1,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 zgradpulse(ACQ_gzlvl2,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		}
	}   
}	//endif ACQ_gt1>0.0




if (BIRDmode[0]== 'h')
{
  	delay(2.0*ACQ_pwx);
	rgpulse(ACQ_pw,v7,rof1,rof1);
	decpower(ACQ_pwxlvl);		
	delay(2.0*ACQ_tau-2.0*rof1-50.0e-9);
	simpulse(2.0*ACQ_pw,2.0*ACQ_pwx,v8,v8,rof1,rof1);
	delay(2.0*ACQ_tau-2.0*rof1);
	rgpulse(ACQ_pw,v9,rof1,0.0);
	decrgpulse(2.0*ACQ_pwx,v8,rof1,0.0);
}
if (BIRDmode[0]== 'b')
{
	delay(pw_XBIP);
	rgpulse(ACQ_pw,v7,rof1,rof1);
	decpower(pwr_XBIP);
	if (pwr_HBIP!=ACQ_tpwr)
	{
	obspower(pwr_HBIP);
	delay(2.0*ACQ_tau-2.0*rof1-100.0e-9);
	}
	else
	{
	delay(2.0*ACQ_tau-2.0*rof1-50.0e-9);
	}				
	simshaped_pulse(shp_HBIP,shp_XBIP,pw_HBIP,pw_XBIP,v8,v8,rof1,rof1);
	if (pwr_HBIP!=ACQ_tpwr)
	{
	obspower(ACQ_tpwr);
	delay(2.0*ACQ_tau-2.0*rof1-50.0e-9);
	}
	else
	{
	delay(2.0*ACQ_tau-2.0*rof1);
	}	
	rgpulse(ACQ_pw,v9,rof1,0.0);
	decshaped_pulse(shp_XBIP,pw_XBIP,v8,rof1,0.0);
}
if (BIRDmode[0]== 'n')
{
	shaped_pulse(shp_a,pw180_a,v7,rof1,rof1);
}


//apply the next PFG level
if (ACQ_gt1>0.0)
{
	if (gradaxis[0]=='t')  
	{ 
		ifzero(v12);
		 rgradient('x',ACQ_gzlvl3);
		 rgradient('y',ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient('x',ACQ_gzlvl4);
		 rgradient('y',ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient('x',ACQ_gzlvl1);
		 rgradient('y',ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient('x',ACQ_gzlvl2);
		 rgradient('y',ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient('x',0.0); 
		 rgradient('y',0.0); 
		 delay(ACQ_gstab);
		endif(v13);
	}   
	else 
	{    
    		if ((gradaxis[0]=='x') || (gradaxis[0]=='y')) 
		{ 
		ifzero(v12);
		 rgradient(gradaxis[0],ACQ_gzlvl3); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl4); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl1); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 rgradient(gradaxis[0],ACQ_gzlvl2); 
		 delay(ACQ_gt1); 
		 rgradient(gradaxis[0],0.0); 
		 delay(ACQ_gstab);
		endif(v13);
		} 
		else 
		{
		ifzero(v12);
		 zgradpulse(ACQ_gzlvl3,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v12);
		ifrtEQ(v12,one,v13);
		 zgradpulse(ACQ_gzlvl4,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,two,v13);
		 zgradpulse(ACQ_gzlvl1,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		ifrtEQ(v12,three,v13);
		 zgradpulse(ACQ_gzlvl2,ACQ_gt1);
		 delay(ACQ_gstab);
		endif(v13);
		}
	}   
}	//endif ACQ_gt1>0.0


/*change offset*/
	if (( (offRes_X180!=0.0) || (offRes_DEC!=0.0)) && (offRes_X180!=offRes_DEC)) decoffset(dof+offRes_DEC);

	delay(tau_r -rof1-rof1-rof3-50.0e-9);


// change DEC power
	decpower(dpwr);	
	obsblank();

ifrtEQ(v11,v14,v13);
	delay(rof1);
	rcvron(); 	//this includes rof3
	delay(rof1);

acquire(droppts1,1.0/sw);
	if (dm[2]=='y')
	{
	decprgon(dseq, 1.0/dmf, dres);
	decon(); decunblank();
	}
//	acquire(npoints/2.0+(droppts1+droppts2),1.0/sw);
	acquire(npoints/2.0,1.0/sw);
//	recoff();					 	
	if (dm[2]=='y')
	{
	decoff();
	decprgoff();
	}
acquire(droppts2,1.0/sw);
recoff();					 	

elsenz(v13);
	delay(rof1);
	rcvron(); 	//this includes rof3
	delay(rof1);
acquire(droppts2,1.0/sw);
	if (dm[2]=='y')
	{
	decprgon(dseq, 1.0/dmf, dres);
	decon(); decunblank();
	}
//	acquire(npoints+(droppts1+droppts2),1.0/sw);
	acquire(npoints,1.0/sw);		
//	recoff();					 	
	if (dm[2]=='y')
	{
	decoff();
	decprgoff();
	}
acquire(droppts2,1.0/sw);
recoff();					 	

endif(v13);

endloop(v11);

delay(0.0001);
	endacq();
	if (lkgate_flg[0] == 'y')  lk_sample(); // turn lock sampling on

}


/************************ ACQ for conventional gHSQC ***********************************/
						
else 
{
	if (homo[A]=='y')
	{   
	status(C);	
	}
	else  
	{  

	setacqmode(WACQ|NZ);	

	obsblank();	
	delay(rof2);
	startacq(alfa);
	status(C);
	acquire(np,1.0/sw);
	recoff();					 	
	endacq();

	if (lkgate_flg[0] == 'y')  lk_sample(); /* turn lock sampling on */	
	}	  
}

}




