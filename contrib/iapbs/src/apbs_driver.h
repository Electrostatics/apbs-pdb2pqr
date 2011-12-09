/**
 *
 *  @file    apbs_driver.h
 *  @author  Robert Konecny
 *  @brief   Header file for the main iAPBS driver
 *
 *  @version $Id: apbs_driver.h 389 2010-03-29 20:18:15Z rok $
 *
 */


/*! \def NATOMS 
    \brief Maximum number of atoms.
*/
#define NATOMS 150000

/**
 * @brief  Wrapper iAPBS function
 * @author Robert Konecny
 *
 * @param nat Number of atoms
 * @param x Atomic coordinate (x)
 * @param y Atomic coordinate (y)
 * @param z Atomic coordinate (z)
 * @param radius Atomic radii
 * @param charge Atomic charges
 * @param r_param Input APBS parameters (real values)
 * @param i_param Input APBS parameters (integer values)
 * @param grid Grid spacing
 * @param dime Grid dimensions
 * @param pdime Grid of processors to be used in calculation
 * @param glen Grid side lengths
 * @param center Grid center
 * @param cglen Coarse grid side lengths
 * @param fglen Fine grid side lengths
 * @param ccenter Coarse grid center
 * @param fcenter Fine grid center
 * @param ofrac Overlap fraction between procs
 * @param dbg Debug verbosity flag
 * @param ionq Mobile ion charge
 * @param ionc Mobile ion concentration
 * @param ionr Mobile ion radius
 * @param esEnergy Electrostatic energy
 * @param npEnergy Non-polar energy
 * @param apbsdx Total electrostatic force per atom (x coordinate)
 * @param apbsdy Total electrostatic force per atom (y coordinate)
 * @param apbsdz Total electrostatic force per atom (z coordinate)
 * @param apbsqfx Fixed charge force (x)
 * @param apbsqfy Fixed charge force (y)
 * @param apbsqfz Fixed charge force (z)
 * @param apbsibx Ionic boundary force (x)
 * @param apbsiby Ionic boundary force (y)
 * @param apbsibz Ionic boundary force (z)
 * @param apbsnpx Non-polar force (x)
 * @param apbsnpy Non-polar force (y)
 * @param apbsnpz Non-polar force (z)
 * @param apbsdbx Dielectric boundary force (x)
 * @param apbsdby Dielectric boundary force (y)
 * @param apbsdbz Dielectric boundary force (z)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int apbsdrv_(int *nat, 
		      double x[NATOMS], 
		      double y[NATOMS], 
		      double z[NATOMS], 
		      double radius[NATOMS], 
		      double charge[NATOMS], 
		      double r_param[9], 
		      int i_param[25],
		      double grid[3],
		      int dime[3], 
		      int pdime[3], 
		      double glen[3], 
		      double center[3], 
		      double cglen[3],
		      double fglen[3], 
		      double ccenter[3], 
		      double fcenter[3], 
		      double *ofrac, 
		      int *dbg, 
		      double ionq[MAXION], 
		      double ionc[MAXION], 
		      double ionr[MAXION], 
		      double esEnergy[NOSH_MAXCALC],
		      double npEnergy[NOSH_MAXCALC],
		      double apbsdx[NATOMS], 
		      double apbsdy[NATOMS], 
		      double apbsdz[NATOMS],
		      double apbsqfx[NATOMS], 
		      double apbsqfy[NATOMS],
		      double apbsqfz[NATOMS],
		      double apbsibx[NATOMS], 
		      double apbsiby[NATOMS],
		      double apbsibz[NATOMS],
		      double apbsnpx[NATOMS], 
		      double apbsnpy[NATOMS],
		      double apbsnpz[NATOMS],
		      double apbsdbx[NATOMS], 
		      double apbsdby[NATOMS],
		      double apbsdbz[NATOMS]);

/**
 * @brief  Calculate forces from MG solution
 * @author Robert Konecny (based on forceMG)
 *
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param  mgparm  MG-specific parmaeters
 * @param  pmg  MG object
 * @param  nforce 0 => no forces, 1 => net forces, >1 => number of
 *                       forces (1 per atom)
 * @param  atomForce Pointer to array of force objects
 * @param  alist  List of atom lists
 * @param  debug verbosity flag
 * @return  1 if successful, 0 otherwise */
VEXTERNC int iforceMG(Vmem *mem, 
		      NOsh *nosh, 
		      PBEparm *pbeparm,  
		      MGparm *mgparm,
		      Vpmg *pmg, 
		      int *nforce, 
		      AtomForce **atomForce, 
		      Valist *alist[NOSH_MAXMOL], 
		      int debug);

/**
* @brief  Combine and pretty-print energy data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
//VEXTERNC double getElecEnergy(

double getElecEnergy(
	Vcom *com, /** Communications object */
	NOsh *nosh, /** Parameters from input file */
	double totEnergy[NOSH_MAXCALC], /** Array of energies
					  from different calculations */
	int iprint /** Index of energy statement to print */
	);


/**
* @brief Creates APBS input string
* @author Robert Konecny
* @param r_param Input APBS parameters (real values)
* @param i_param Input APBS parameters (integer values)
* @param grid Grid spacing
* @param dime Grid dimensions
* @param ionq Mobile ion charge
* @param ionc Mobile ion concentration
* @param ionr Mobile ion radius
* @param glen Grid side lengths
* @param center Grid center
* @param cglen Coarse grid side lengths
* @param fglen Fine grid side lengths
* @param ccenter Coarse grid center
* @param fcenter Fine grid center
* @param ofrac Overlap fraction between procs
* @param pdime Grid of processors to be used in calculation
* @param debug Debug verbosity level.
* @return the input string */
char *setupString(double r_param[9], int i_param[25], double grid[3], 
	int dime[3], double ionq[MAXION], double ionc[MAXION],
	double ionr[MAXION], double glen[3], double center[3],
	double cglen[3], double fglen[3],
	double ccenter[3], double fcenter[3], double *ofrac, 
	int pdime[3], int debug);
