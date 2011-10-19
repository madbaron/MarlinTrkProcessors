#ifndef SimplePlanarDigiProcessor_h
#define SimplePlanarDigiProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>

#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;


/** ======= SimplePlanarDigiProcessor ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The plannar geometry should be either VXD, SIT or SET described using ZPlannarLayout
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param PointResolutionRPhi_Inner Point resolution perpendicular to the ladder (in mm) <br>
 * (default value 0.004) <br>
 * @param PointResolutionZ_Inner Point resolution along the ladder (in mm) <br>
 * (default value 0.004) <br>
 * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
 * (default value false) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
 * (default value ILDDetID::VXD) <br>
 * <br>
 * 
 */
class SimplePlanarDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new SimplePlanarDigiProcessor ; }
  
  
  SimplePlanarDigiProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  // find phi in correct range, taken from gear::VXDParameters
  double correctPhiRange( double Phi ) const ;  
  
  //  void fillMapTVolumeChildren(TGeoVolume* volume) ;
  
protected:
  
  std::string _inColName ;
  
  std::string _outColName ;
  
  int _sub_det_id ;
  
  int _nRun ;
  int _nEvt ;
  
  float _pointResoRPhi ;
  float _pointResoZ, ;
  
  bool _ladder_Number_encoded_in_cellID;
  
  gsl_rng* _rng ;
  
  
} ;

#endif



