#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "modules.h"
#include "cpbconcz2.h"

void print_array_values( const double *values );

/*
 * Helper functions.  Most of these are artifacts of refactoring.
 */


void
initXYZR( double xyzr[MAXATOMS][XYZRWIDTH] )
{
    for(int i = 0; i < XYZRWIDTH; i++)
    {
        for(int j = 0; j < MAXATOMS; j++)
        {
            xyzr[j][i] = 0.0;
        }
    }
}


template<typename T>
void 
initValues( T* array, int numelems, T initvalue )
{
    for(int i = 0; i < numelems; i++)
    {
        array[i] = initvalue;
    }
}

/*
 * Initialize all elements of the given array to hold the initvalue.
 * Can easily be used to initialize multidimensional arrays as well.
void
initReals( double *array, int numelems, double initvalue )
{
    for(int i = 0; i < numelems; i++)
    {
        array[i] = initvalue;
    }
}
 */

void 
writePressGamaRms()
{
    std::ofstream f;
    std::string spaces("    ");
    f.open("press_gama_rms.txt");
    f << spaces << "pressure" << spaces << "gamma" << spaces << "rmse" << spaces << "rmse_error";
    f << std::endl;
    f.close();
}

void
writeOutputMolFiles(int nmol)
{
    std::string spaces("    ");
    for(int i = 1; i <= nmol; i++)
    {
        std::ofstream f;
        std::stringstream s;
        s << "output_mol" << 9+i << ".txt"; // e.g. output_mol10.txt
        f.open( s.str().c_str() );
        f << spaces << "pressure" << spaces << "gamma" << spaces << "nonpolar" << spaces << "electro " << spaces << " totl_energy" << "  area" << spaces << "volume" << spaces << "uattint" << spaces << "expt" << spaces << "error";
        f << std::endl;
        f.close();
    }
    
}

void 
writeSupportingFiles(int nmol)
{
    writePressGamaRms();
    writeOutputMolFiles(nmol);
}


int
calculateMDArrayOffset( std::vector<int> dims, std::vector<int> elementIdx )
{
    int offset = 0, dimlen = dims.size();
    for(int i = 0; i < dimlen; i++)
    {
        int tmp = 1;
        for(int j = i + 1; j < dimlen; j++)
        {
            tmp *= dims[j];
        }
        offset += tmp * elementIdx[i];
    }
    return offset;
}

/*
 * Given an n-dimensional array arrayBase, with dimensions specified by dims,
 * sets the element at the indices given by elementIdx to the value val
 */
template <typename T> 
void
setMDArrayElement( T* arrayBase, std::vector<int> dims, std::vector<int> elementIdx, T val )
{
    int offset =  calculateMDArrayOffset( dims, elementIdx );
    arrayBase[offset] = val;
}

template <typename T>
T
getMDArrayElement( T* arrayBase, std::vector<int> dims, std::vector<int> elementIdx)
{
    int offset =  calculateMDArrayOffset( dims, elementIdx );
    return arrayBase[offset];

}

template<typename T>
std::vector<T>
arrayToVector( T array[], int numelems )
{
    return std::vector<T>( (T*)array, array + numelems );
}

int
numberOfLines( std::string fileName )
{
    std::ifstream file(fileName.c_str());
    return std::count(std::istreambuf_iterator<char>(file),
              std::istreambuf_iterator<char>(), '\n');
}


void
normalizeSurfuAndEps( double* surfu, double* eps, int nx, int ny, int nz, double epsilons, double epsilonp )
{
    int aDim[] = {nz, ny, nx};
    std::vector<int> dims( aDim, aDim + sizeof(aDim) / sizeof(int));
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                int aElem[] = {k,j,i};
                std::vector<int> elems( aElem , aElem + 3 );
                if(getMDArrayElement( surfu, dims, elems ) > 1000.0)
                {
                    setMDArrayElement( surfu, dims, elems, 1000.0 );
                }
                if(getMDArrayElement( surfu, dims, elems ) < 0.0)
                {
                    setMDArrayElement( surfu, dims, elems, 0.0 );
                }
                double s = getMDArrayElement( surfu, dims, elems );
                double newEps = epsilonp + (epsilons-epsilonp) * ((1000.0 - s) / 1000.0);
                setMDArrayElement (eps, dims, elems, newEps );
            }
        }
    }

}

/*
 * Compute soleng1 and soleng2 (solvation).  This is an artifact of refactoring the original F90 code.
 * Parameters:
 *     soleng:        The soleng variable to set (compute).  For soleng1, the phi array should be phi; for 
 *             soleng2, the phiarray should be phivoc. 
 *     phi:        Either phi or phivoc array, depending on which soleng variable we are computing
 *     phiDims:    Vector holding the dimensions of the phi array
 *     nchr:        The number of atoms
 *     loc_qt:        The loc_qt array.  This is a [3][8][nchr] int array, but must be passed as a pointer,
 *             since nchr is a variable.
 *    charget:    charget array.  This is an [8][nchr] int array.
 */
void
computeSoleng( double& soleng, double* phi, std::vector<int> phiDims, double* charget, int nchr, int *loc_qt )
{
    int locqtDims[] = { 3,8,nchr };
    std::vector<int> locqtDimv = arrayToVector(locqtDims, 3);
    int chargetDims[] = {8,nchr};
    std::vector<int> chargetDimv = arrayToVector(chargetDims, 2);
    soleng = 0.0;

    for(int iind = 0; iind < nchr; iind++)
    {
        for(int jind = 0; jind < 8; jind++)
        {
            int idx[] = {0,jind,iind}; 
            std::vector<int>  idxv = arrayToVector(idx, 3);

            int i1 = getMDArrayElement( loc_qt, locqtDimv, idxv );
            idxv[0] = 1;
            int j1 = getMDArrayElement( loc_qt, locqtDimv, idxv );
            idxv[0] = 2;
            int k1 = getMDArrayElement( loc_qt, locqtDimv, idxv );


            int chargetIdx[] = {jind, iind}; 
            std::vector<int> chargetIdxv = arrayToVector(chargetIdx, 2);

            double c = getMDArrayElement( charget, chargetDimv, chargetIdxv );
            if( c != 0 )
            {
                int idx[] = {k1 - 1,j1 - 1,i1 - 1};
                std::vector<int> idxv = arrayToVector(idx, 3);
                
                soleng += 0.5 * getMDArrayElement( phi, phiDims, idxv ) * c;
            }


        }
    }
}

double maxrms( double* sumpot, double* expv, double* elec, double gama, int nmol )
{
    double rms = 0.0;    
    double snum[50]; // constants are f90 legacy
    double err[50];
    for(int i = 0; i < nmol; i++)
    {
        snum[i] = gama * sumpot[i] + elec[i];
        err[i] = snum[i] - expv[i];
        rms += err[i] * err[i];
    }
    return sqrt(rms/nmol);
}

int loadData(std::ifstream& f,    //<i
         int imord, //<i
         int ffmodel, //<i
         double radexp, //<i
         double* expv,  //<o expv[imord]
         double xyzr[MAXATOMS][XYZRWIDTH], //<o
         double* pqr,//<o
         double* ljepsilon  //<o
         ){
    std::string molecule;
    double val;
    f >> molecule >> val;
    std::cout << "expSolv:\t" << val << "\t" << molecule << std::endl;
    expv[imord] = val;

    std::string atomFileName = molecule + ".xyzr";
    std::ifstream atomFile;
    atomFile.open( atomFileName.c_str() );

    int natm = 0;
    initXYZR(xyzr);
    switch(ffmodel)
    {
        case 1:
            while( !atomFile.eof() )
            {
                atomFile >> xyzr[natm][0] >> xyzr[natm][1] >> 
                    xyzr[natm][2] >> xyzr[natm][3] >> pqr[natm];
                atomFile.ignore(); // skip remainder of line
                xyzr[natm][3] *= radexp;
                natm++;
            }
            break;
        case 2:
        default:
            initValues(ljepsilon, MAXATOMS, 0.0);
            while( !atomFile.eof() )
            {
                atomFile >> xyzr[natm][0] >> xyzr[natm][1] >> 
                    xyzr[natm][2] >> xyzr[natm][3] >> pqr[natm] >> ljepsilon[natm];
                atomFile.ignore(10000000, '\n'); // skip remainder of line
                xyzr[natm][3] *= radexp;
                natm++;
            /*
                atomFile >> a->pos[0] >> a->pos[1] >> a->pos[2] >> a->radius >> pqr[natm] >> epsilon[natm];
                a->radius *= radexp;
                natm++;
                */
            }

            break;
    }
    natm--;
    return natm;
}

GeoflowOutput geoflowSolvation(double xyzr[MAXATOMS][XYZRWIDTH], int natm, double dcel, int ffmodel, double extvalue, double* pqr, int maxstep, double crevalue, int iadi, double tottf, double* ljepsilon, double alpha, int igfin, double epsilons, double epsilonp, int idacsl, double tol, int iterf, double tpb, int itert, double pres, double gama, double tauval, double prob, int vdwdispersion, double sigmas, double density, double epsilonw){
    double elec;
    double area = 0.0, volume = 0.0, attint = 0.0;
    double *phi, *phix, *phivoc, *surfu;

    comdata.dcel = dcel;
    comdata.deltax = comdata.deltay = comdata.deltaz = dcel;
    comdata.pi = acos(-1.0);
    
    lj.ffmodel = ffmodel;
    lj.tauval = tauval; 
    lj.prob = prob;
    lj.vdwdispersion = vdwdispersion;
    lj.sigmas = sigmas;
    lj.density = density;
    lj.epsilonw = epsilonw;
    lj.roro = density / gama;
    double potcoe = 1 / gama;
    lj.conms = pres / gama;

    domainini(xyzr, natm, extvalue);
 
    int nchr = natm;
    int width = 3*8*nchr;
    double corlocqt[width],
        charget[8*nchr];
    int loc_qt[width];
    initValues((int*)loc_qt, width, 0);
    initValues((double*)corlocqt, width, 0.0);
    initValues((double*)charget, width, 0.0);
    for(int iatm = 1; iatm <= nchr; iatm++)
    {
        chargedist((double*)xyzr, pqr, nchr, (double*)charget, (double*)corlocqt, (int*)loc_qt, iatm);
    }
 
    comdata.xc = new double[comdata.nx];
    comdata.yc = new double[comdata.ny];
    comdata.zc = new double[comdata.nz];
 
    for(int i = 0; i < comdata.nx; i++)
    {
        comdata.xc[i] = comdata.xleft + (i-1)*comdata.dcel;
    }
 
    for(int i = 0; i < comdata.ny; i++)
    {
        comdata.yc[i] = comdata.yleft + (i-1)*comdata.dcel;
    }
 
    for(int i = 0; i < comdata.nz; i++)
    {
        comdata.zc[i] = comdata.zleft + (i-1)*comdata.dcel;
    }
 
    /* These all correspond to arrays in original f90 code.
     * These are actually 3d arrays, but we define them as a 1d array 
     * since the dimensions are not known at compile time: 
     
    double surfu[comdata.nz][comdata.ny][comdata.nx],
        phi[comdata.nz][comdata.ny][comdata.nx],
        phix[comdata.nz][comdata.ny][comdata.nx],
             eps[comdata.nz][comdata.ny][comdata.nx],
             phivoc[comdata.nz][comdata.ny][comdata.nx],
             bg[comdata.nz][comdata.ny][comdata.nx],
             bguni[comdata.nz][comdata.ny][comdata.nx];
         */
         int arrayLengths = comdata.nz * comdata.ny * comdata.nx;
         int aDim[] = { comdata.nz, comdata.ny, comdata.nx };
         std::vector<int> dims( aDim, aDim + sizeof(aDim) / sizeof(int) );
 
         int solvWidth = maxstep + 1;
         double solv[solvWidth];
         initValues( (double*)solv, solvWidth, 0.0 );
         double diffEnergy = 100;
 
         int iloop = 0, ipath = 0; double tott = 0.0;
         phi= new double[arrayLengths];
         phix = new double[arrayLengths];
         phivoc = new double[arrayLengths];
         surfu = new double[arrayLengths] ;
         initValues( phi, arrayLengths, 0.0 );
         initValues( phix, arrayLengths, 0.0 );
         initValues( phivoc, arrayLengths, 0.0 );
         initValues( surfu, arrayLengths, 0.0 );
 
         double *eps = new double[arrayLengths],
             *bg = new double[arrayLengths],
             *bguni = new double[arrayLengths];
         while( (iloop < maxstep) && (diffEnergy > crevalue) ){
             iloop++;
            
             double deltat=0;//this is wrong for adi...
             initValues( eps, arrayLengths, 0.0 );
             if(!iadi){
                 deltat = comdata.dcel * comdata.dcel / 4.5;
             }
             if(ipath == 0) /* will always be true, but this is how
                  * it was in the original F90 code.  
                  * TODO: rewrite this to remove the ipath variable (?) */
             {
                 tott = fmax(tottf - iloop + 1, 1.0);    
             }else{
                 if (iloop == 1){
                     tott = tottf;
                 }else{
                     tott = iloop * comdata.dcel * comdata.dcel / 4.5;
                 }
             }
             area = volume = attint = 0.0;
             yhsurface(xyzr, ljepsilon, natm, tott, deltat, (double*)phix, (double*)surfu, iloop, area, volume, attint, alpha, iadi, igfin);
             normalizeSurfuAndEps( surfu, eps, comdata.nz, comdata.ny, comdata.nx, epsilons, epsilonp );
             if(iloop == 1)
             {
                 seteqb( (double*)bg, xyzr, (double*)pqr, natm, (double*)charget, (double*)corlocqt, &epsilons);
             }
          
             if(idacsl == 1)
             {
                  for(int i = 1; i <= comdata.nx; i++){
                  for(int j = 1; j <= comdata.ny ; j++){
                  for(int k = 1; k <= comdata.nz;k++){
                      double x = xvalue(i),    
                             y = yvalue(j),
                             z = zvalue(k);
                      int idx[] = { k - 1, j - 1, i - 1};
                      std::vector<int> idxv( idx, idx + sizeof(idx) / sizeof(int) );
                      setMDArrayElement( eps, dims, idxv, x+y+z );
      
                      int ijk = (i-1) * comdata.nz * comdata.ny + (j-1) * comdata.nz + k;
                      if( i == 1 || i == comdata.nx ||
                          j == 1 || j == comdata.ny ||
                          k == 1 || k == comdata.nz )
                          {
                              bg[ijk - 1] = cos(x) * cos(y) * cos(z);    
                          }else{
                              bg[ijk - 1] = 
(-sin(x)*cos(y)*cos(z)-cos(x)*sin(y)*cos(z)-cos(x)*cos(y)*sin(z))-3.0*(x+y+z)*cos(x)*cos(y)*cos(z);
                     }
             }}}
        }

        int iter = 1000; double fpb, titer = 0.0;
        pbsolver( eps, phi, bg, comdata.nx, comdata.ny, comdata.nz, comdata.dcel, tol, iter);
        if(iloop == 1 )
        {
            fpb = titer;
            iterf=iter;
        }
        tpb = tpb + iter;
        itert += iter;


        if(idacsl == 1)
        {
            //maxerr = 0.0;    
            for(int i = 1; i <= comdata.nx; i++){
            for(int j = 1; j <= comdata.ny ; j++){
            for(int k = 1; k <= comdata.nz;k++){
                double     x = xvalue(i),
                     y = yvalue(j),
                    z = zvalue(k);    
                int idx[] = { k - 1, j - 1, i - 1 };
                std::vector<int> idxv( idx, idx + sizeof(idx) / sizeof(int) );
                double phival = getMDArrayElement( phi, dims, idxv );
                double err = abs(phival) - cos(x)*cos(y)*cos(z);
              //  maxerr = max(maxerr, err);

            }}}
            exit(0);
        }

        if(iloop==1)
        {
            for(int ii = 1; ii <= comdata.nx; ii++){
            for(int jj = 1; jj <= comdata.ny; jj++){
            for(int kk = 1;  kk <= comdata.nz; kk++){
                int ijk = (ii-1)*comdata.ny*comdata.nz + (jj-1)*comdata.nz + kk;
                if( ii<2 || ii>comdata.nx-1 || jj<2 || jj>comdata.ny-1 || kk<2 || kk>comdata.nz-1 )
                  {
                      bguni[ijk - 1] = bg[ijk - 1]*epsilons/epsilonp;
                  }else{
                      bguni[ijk - 1] = bg[ijk - 1];
                  }
                int idx[] = { kk-1,jj-1,ii-1 };
                std::vector<int> idxv(idx, idx + sizeof(idx) / sizeof(int) );
                setMDArrayElement( eps, dims, idxv, epsilonp );
            }}}
        }

        initValues( eps, arrayLengths, epsilonp );

        pbsolver( eps, phivoc, bg, comdata.nx, comdata.ny, comdata.nz, comdata.dcel, tol,iter);

        double weit = 1.0/(2.0 * comdata.dcel);
        for(int ix = 2; ix <= comdata.nx - 1; ix++){
        for(int iy = 2; iy <= comdata.ny - 1; iy++){
        for(int iz = 2; iz <= comdata.nz - 1; iz++){
            double phixx,phixy,phixz;
            {
                int idx1[] = {iz-1,iy-1,ix},
                    idx2[] = {iz-1,iy-1,ix-2};

                std::vector<int> idx1v = arrayToVector( idx1, 3 ),
                         idx2v = arrayToVector( idx2, 3 );
                phixx = getMDArrayElement( phi, dims, idx1v ) -
                           getMDArrayElement( phi, dims, idx2v ) * weit;
            }
            {
                int idx1[] = {iz-1,iy,ix-1},
                    idx2[] = {iz-1,iy-2,ix-1};

                std::vector<int> idx1v = arrayToVector( idx1, 3 ),
                         idx2v = arrayToVector( idx2, 3 );
                phixy = getMDArrayElement( phi, dims, idx1v ) -
                           getMDArrayElement( phi, dims, idx2v ) * weit;
        
            }
            {
                int idx1[] = {iz,iy-1,ix-1},
                    idx2[] = {iz-2,iy-1,ix-1};

                std::vector<int> idx1v = arrayToVector( idx1, 3 ),
                         idx2v = arrayToVector( idx2, 3 );
                phixz = getMDArrayElement( phi, dims, idx1v ) -
                           getMDArrayElement( phi, dims, idx2v ) * weit;
        
            }
            int idx[] = {iz-1,iy-1,ix-1};

            std::vector<int> idxv = arrayToVector( idx, 3 );
            double p = 0.5*(epsilons - epsilonp) * (phixx*phixx + phixy*phixy + phixz*phixz) * potcoe;
            setMDArrayElement( phix, dims, idxv, p );
        

        }}}

        // solvation
        double soleng1, soleng2;
        soleng1 = soleng2 = 0.0;
        computeSoleng( soleng1, phi,    dims, (double*)charget, nchr, (int*)loc_qt );
        computeSoleng( soleng2, phivoc, dims, (double*)charget, nchr, (int*)loc_qt );
        solv[iloop - 1] = (soleng1 - soleng2) * 332.0716; 
        elec = solv[iloop - 1];
        solv[iloop - 1] = elec + gama * (area + volume * lj.conms + attint * lj.roro);
        if(iloop > 1)
        {
            diffEnergy  = abs( (solv[iloop - 1] - solv[iloop - 2]) );
        }
    }



    double sumpot = area + volume*lj.conms + attint*lj.roro;
    double nonpolarSolvation = sumpot*gama;
    double totalSolvation = nonpolarSolvation + elec; 

    return (GeoflowOutput){area, volume, attint, sumpot, totalSolvation, nonpolarSolvation, elec};
}

void
processAtomsFile( std::string fileName, int ffmodel, int radexp, double extvalue, int maxstep, double crevalue, int iadi, double tottf, double alpha, int igfin, double epsilons, double epsilonp, int idacsl, double tol, double &tpb, int &iterf, int &itert, double pres, double gama, double dcel, double tauval, double prob, int vdwdispersion, double sigmas, double density, double epsilonw)
{
    std::ifstream f;
    f.open(fileName.c_str());

    int nmol = numberOfLines( fileName );
    for(int imord = 0; imord < nmol ; imord++)
    {
        double xyzr[MAXATOMS][XYZRWIDTH];
        double ljepsilon[MAXATOMS];
        double pqr[MAXATOMS];
        double expv[100]; // hardcoded 100 constant, as per interface to F90 code :(

        int natm = loadData(f, imord, ffmodel, radexp, expv, xyzr, pqr, ljepsilon);
        
        GeoflowOutput gf = geoflowSolvation(xyzr, natm, dcel, ffmodel, extvalue, pqr, maxstep, crevalue, iadi, tottf, ljepsilon, alpha, igfin, epsilons, epsilonp, idacsl, tol, iterf, tpb, itert, pres, gama, tauval, prob, vdwdispersion, sigmas, density, epsilonw);

        std::cout << "totalSolv:\t" << gf.totalSolvation << "\t";
        std::cout << "nonpolar: " << gf.nonpolarSolvation << "\t";
        std::cout << "electro: " << gf.elecSolvation << "\n" << std::endl;
    }
//    std::cout << "RMS Error: " << maxrms( (double*)sumpot, (double*)expv, (double*)elec, gama, nmol) << std::endl;
}

void 
setPresStep( double &pres_step, double pres )
{
    if( pres < 0.001 )
    {
        pres_step = 0.0001;    
    }else if( (pres > 0.001) && (pres < 0.01) )
    {
        pres_step = 0.001;
    }else // if( pres >= 0.001)
    {
        pres_step = 0.005;
    }
}

void 
setGamaStep( double &gama_step, double gama )
{
    if( (gama >= 0.00001) && (gama < 0.0001) )
    {
        gama_step = 0.00001;
    }else if(( gama >= 0.0001) && (gama < 0.001) )
    {
        gama_step = 0.0001;
    }else if(( gama >= 0.001) && ( gama < 0.01 ) )
    {
        gama_step = 0.001;
    }else if((gama > 0.01) && (gama <= 0.055) )
    {
        gama_step = 0.005;
    }else // gama > 0.055
    {
        gama_step = 0.005;
    }
}

void
pbconcz2(
    // These parameters correspond directly to those read in via the datafiles (fort.12 and 17set.txt)
    // in the original Fortran code    
    int     nmol,
    double  pres_i,
    double  gama_i,
    int    npiter,
    int    ngiter,
    double    tauval,
    double    prob,
    int    ffmodel, // 1 for ZAP-9/AM1-BCCv1; 2 for OPLS/AA
    double    sigmas,     // Angstrom (radius of water molecule based on LJ parameter sigma)
    double    epsilonw,// epsilon parameter of 0 (kcal/mol) of water molecule
    int    vdwdispersion,// 1(on) or 0(off)- previously called REPULSIVE
    double    extvalue, // (distance atom surface and box boundary)
//    int    iprec,    // flag to indicate the usage of preconditioner iprec =1 (yes); 0 (no)
//    int    istep,
    int    iadi,    // 0 for explicit scheme; 1 for ADI scheme
    double        alpha, //  weight of previous solution to change the next solution in geometry flow
//    int        ipbin, //  start guess for PB 1; inherit '0'
    double        tol, 
    double    tottf, //  total time
    double        dcel,
    int        maxstep,
    double        epsilons,
    double        epsilonp,
    int        radexp,
    double        crevalue,
    int        idacsl, //  0 for solvation force calculation; 1 or accuracy test
    double         density     //  (use 0.03346) 
)
{
    double pres = pres_i;
//    writeSupportingFiles(nmol);

    double pres_step;
    for(int indpres = 1; indpres <= npiter; indpres++)
    {
        setPresStep( pres_step, pres );    
        if(indpres > 1 )
        {
            pres += pres_step;
        }

        double gama = gama_i, 
            gama_step;
        for(int indgama = 1; indgama <= ngiter; indgama++)
        {    
            setGamaStep( gama_step, gama );
            std::cout << "GAMA_STEP= " << gama_step << std::endl;
            if( indgama > 1 )
            {    
                gama += gama_step;
            }
        }

        std::cout << "GAMA = " << gama << std::endl;
//        std::cout << "RORO = " << lj.roro << std::endl;
//        std::cout << "POTCOE = " << potcoe << std::endl;
        std::cout << "PRES = " << pres << std::endl;
//        std::cout << "CONMS = " << lj.conms << std::endl;
        
        nmol = 0;

        int igfin = 1;
        double tpb = 0.0;
        int iterf = 0, itert = 0;

        processAtomsFile( "17set.txt",  ffmodel, radexp, extvalue, maxstep, crevalue, iadi, tottf, alpha, igfin, epsilons, epsilonp, idacsl, tol, tpb, iterf, itert, pres, gama, dcel, tauval, prob, vdwdispersion, sigmas, density, epsilonw);
    }
}

//int
//main()
//{
//    pbconcz2(17,    // nmol
//         0.03,    // pres_i
//         0.08,    // gama_i
//         1,    // npiter
//         1,    // ngiter
//         1.40,    // tauval
//         0.0,    // prob
//         2,    // ffmodel
//         1.5828,// sigmas
//         0.1554,// epsilonw
//         1,    // vdwdispersion
//         1.90,    // extvalue
////         0,    // iprec
////         10,    // istep
//         0,    // iadi
//        0.50,     // ALPHA
////        1,    // IPBIN
//        1e-5,    // TOL
//        3.5,    // TOTTF
//        0.25,    // DCEL
//        20,    // MAXSTEP
//        80.00,     // EPSILONS
//        3.00,     // EPSILONP
//        1,     // RADEXP
//        0.01,     // CREVALUE
//        0,     // idacsl
//        0.03346 //density (use 0.03346) 
//        );
//}

