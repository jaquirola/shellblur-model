#include <xsTypes.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <XSUtil/Utils/XSutility.h>
#include <cfortran.h>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB6(DOSHBLUR,doshblur,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
#define DOSHBLUR(ear,ne,param,ifl,photar,tmp) \
    CCALLSFSUB6(DOSHBLUR,doshblur,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV, \
      ear,ne,param,ifl,photar,tmp)

extern "C" void shblur(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Wrapper to handle dynamic memory allocation.  Actual work
   // is done in doshblur

   using namespace XSutility;

   int ne = static_cast<int>(energyArray.size()) - 1;
   auto_array_ptr<float> apTmp(new float[ne]);

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   auto_array_ptr<float> apEar(ear);
   auto_array_ptr<float> apPars(pars);
   auto_array_ptr<float> apPhotar(photar);
   auto_array_ptr<float> apPhoter(photer);

   DOSHBLUR(ear, ne, pars, spectrumNumber, photar, apTmp.get());

   XSFunctions::floatFluxToStl<float>(photar, photer, ne, false, flux, fluxErr);       
}
