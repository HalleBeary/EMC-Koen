#include "treeStructure.hpp"

#include "resetParticleLocation.hpp"
#include "EMCGenericWriter.hpp"
#include "free_flight_scatter.hpp"
#include "scat_table.hpp"
#include "acoustic.hpp"
#include "polaroptical.hpp"
#include "output.hpp"
#include "init_kspace.hpp"

/*

#ifdef _OPENMP
#include "ScalFMM/Core/FFmmAlgorithmThread.hpp"
#include "ScalFMM/Core/FFmmAlgorithmPeriodic.hpp"
#include "ScalFMM/Core/FFmmAlgorithmSectionTask.hpp"
#else
#include "ScalFMM/Core/FFmmAlgorithm.hpp"
#endif


#ifdef _OPENMP
//using FmmClass = FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
using FmmClass    = FFmmAlgorithmSectionTask<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> ;
using FmmClassPer = FFmmAlgorithmPeriodic<FReal,OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#else
using FmmClass  = FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#endif

*/


int main(int argc, char* argv[])
{ 
  // Define the mover and arranger classes dealing with moving and rearranging particles
  typedef FBasicParticleContainerIndexedMover<FReal,OctreeClass, ContainerClass> MoverClass;
  typedef FOctreeArranger<FReal, OctreeClass, ContainerClass, MoverClass> ArrangerClass;

  const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZQ100.bfma" );
  const std::string filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
  const int TreeHeight       = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
  const int SubTreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  bool periodicCondition     = false ; //! Periodic FMM. Repeat simulation volume in all directions. 
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options))
      periodicCondition = true;
  const unsigned int aboveTree = FParameters::getValue(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options, -1);

#ifdef _OPENMP
  const unsigned int NbThreads = periodicCondition ? 1: 
                                 FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << NbThreads << " threads.\n" << std::endl;
#else
  const int NbThreads =  1;
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif

  // Set up timers
  FTic time;
  FTic counter;

  // Print out some information on the input file and octree setup
  outputFileOctreeInfo(TreeHeight, SubTreeHeight, aboveTree, NbThreads, filename, periodicCondition);

  //! Currently calling loader twice, first time to populate octree, second time to add momenta. Reason: Right now you need populated octree to call container particle class which you need to add momenta.
  FFmaGenericLoader<FReal> loader(filename); // Open particle file
  geometryClass DevGeometry(&loader);
  
  OctreeClass tree(TreeHeight, SubTreeHeight, DevGeometry.BoxWidth, DevGeometry.BoxCenter);
  // Read particles and insert them in octree
  time.tic();
  //loadParticlesInOctree(&loader, TreeHeight, SubTreeHeight, &tree);
  time.tac();
  std::cout << "Done  " << "(@Creating and Inserting Particles = " << time.elapsed() << " s) ." << std::endl;


  EMCGenericLoader<FReal> testloader(filename);
  FSize Nb = testloader.getNumberOfParticles();
  FPoint<FReal> momenta;
  FReal *TotalMomenta;
  TotalMomenta = new FReal[3*Nb] ;
  LoadParticles(&testloader, TreeHeight, SubTreeHeight, &tree, &momenta, TotalMomenta);





//------------------------------------------------------------------------------------------------------
  //initialize meterial /scattering parameters

  mat_paramClass mat_par;
  scat_paramClass scat_par(0,1); // Acoustic, polar scattering

  ScatteringTable(&scat_par, &mat_par);


  // -------------------------------------initialization  k-space ---------------------------------------
  
  OctreeClass::Iterator octreeIterator(&tree);
  octreeIterator.gotoBottomLeft();

  int i, N, H;
  N = loader.getNumberOfParticles(); // total number of particles.
  i = 0; // counter
  H = 0*N; 
  
  ContainerClass *particles2;
	

    do{  
    particles2 = octreeIterator.getCurrentListTargets();

    for(FSize idxPart = 0; idxPart < particles2->getNbParticles() ; ++idxPart){ 

//TODO K_inreader and kspace init need to be in the same function

       k_inreader(particles2, TotalMomenta,  &scat_par, &mat_par, idxPart, i);
    //  kspace_init(particles2, &scat_par, &mat_par, idxPart); // Add initial momenta
	
      if (i > 0 && i <= H)
        ChargeChanger(particles2, idxPart);

     // realspace()

      i++;
      }
  
  } while(octreeIterator.moveRight());

  
  delete[] TotalMomenta;


  // -- Rearrange particles into octree after grid initialization .

//  ArrangerClass arrange2(&tree);
//  arrange2.rearrange();

 
   
//------------------------------------------------------------------------------------------------------
  int totaltime = 10001;
  int tpulse = totaltime / 10;
  int twrite = totaltime / 100;
  double pulsesize = 0.9;

  for(int tsim = 0 ; tsim < totaltime; tsim ++) // total simulation time in fs
  {
    std::cout << "Start of " << tsim << "'th FMM run" << std::endl;
 
    std::cout << "eps_infty: " << mat_par.get_eps_infty() << "\n \n";


   ////////////////////////////////////////////////////////////////////
    //
    //    Execute FMM Algorithm
    //
    ////////////////////////////////////////////////////////////////////

 // { // -----------------------------------------------------
    std::cout << "\n" << interpolationType << "  FMM (ORDER= "<< ORDER << ") ... " << std::endl;
    
    // reset all particle forces
    resetParticleLocation(&tree);

    const MatrixKernelClass    MatrixKernel;
    time.tic();

    std::unique_ptr<KernelClass> kernelsNoPer(new KernelClass(TreeHeight, DevGeometry.BoxWidth, DevGeometry.BoxCenter,&MatrixKernel));
    //
    FmmClass algoNoPer(&tree, kernelsNoPer.get());
    // periodic FMM algorithm
    FmmClassPer algoPer(&tree, aboveTree);
    KernelClass kernelsPer(algoPer.extendedTreeHeight(), algoPer.extendedBoxWidth(),
                           algoPer.extendedBoxCenter(),&MatrixKernel);
    algoPer.setKernel(&kernelsPer);
    //
    FAbstractAlgorithm * algorithm  = nullptr;
    FAlgorithmTimers   * timer      = nullptr;
    if(! periodicCondition) { // Non periodic case
        algorithm  = &algoNoPer ;
        timer      = &algoNoPer ;
      }
    else {                    // Periodic case
        algorithm  = &algoPer ;
        timer      = &algoPer ;
	std::cout << "\n periodic is on \n" ;
      }

    //


 // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //  algorithm->execute(FFmmFarField);
    
 // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //  algorithm->execute(FFmmNearField);

    algorithm->execute();

    time.tac();
    std::cout << "Timers Far Field \n"
              << "P2M " << timer->getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
              << "M2M " << timer->getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
              << "M2L " << timer->getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
              << "L2L " << timer->getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
              << "P2P and L2P " << timer->getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
	      << std::endl;

    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;

  //}
  // -----------------------------------------------------
std::cout << " Particle positions updated according to new field\n" << std::endl;
  // Arrange Section--------------------------------------------------------------------------------*//
  
  ArrangerClass arrange(&tree);
  
        std::cout << "Working on particles ..." << std::endl;
        counter.tic();

        {
         //   OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            
            do{  
                ContainerClass* particles = octreeIterator.getCurrentListTargets();
                // Advance each particle
               //----------------------------- MONTE-CARLO PROCEDURE ----------------------------------------------------------



                for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){ 

                  free_flight_scatter(DevGeometry, particles, &scat_par, &mat_par, idxPart, tsim, tpulse); 

		 


		
                } // end for Advance each particle
            } while(octreeIterator.moveRight());
        }
        counter.tac();
        std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;

      
        if (tsim == tpulse){
          octreeIterator.gotoBottomLeft();
          i = 0;

          do{
          
          ContainerClass* particles3 = octreeIterator.getCurrentListTargets();

          for(FSize idxPart = 0; idxPart <  particles3->getNbParticles() ; ++idxPart){
	    
	    if (particles3->getPhysicalValues()[idxPart] < 0 ) // Only electrons get excited in pulse
              deltaE(particles3, DevGeometry, &scat_par, &mat_par, idxPart);


	    i++;

            }
            
          }while(octreeIterator.moveRight()  &&  i < pulsesize*N);
        
      //  std::cout << "Pulse initiated. " << std::endl;
        }
	
	// ------------------------- writer -------------------


        if(tsim % twrite  == 0){

  
	  std::string nombre;
		  

          nombre = "../../../../data/p274598/Output/pulseL250K2-10001-3e17_holes-0.2-10/" + std::to_string(tsim) + ".fma";
          EMCGenericWriter<FReal> writeur(nombre);
          writeur.writeDataFromOctree(&tree, loader.getNumberOfParticles());
          
        }

	
        // --------------------------------------------
	
      
        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Arrange ..." << std::endl;
        counter.tic();
        arrange.rearrange();
        counter.tac();
        std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
        //    }
    std::cout << "Test ..." << std::endl;
    counter.tic();

  //-------------------------------------------------------------------------------*//
  // Perhaps output statement of number of fmm runs.
   std::cout <<std::endl<<"End of FMM run "<< tsim << std::endl<<std::endl;
  }   // -end of for-loop FMM
  //---------------------------------------------------------------------------

  //----------------------------------------------------------------------------

  // ------------------------------------------
  // Write output using EMCGenericWriter class
  // ------------------------------------------

  if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
    std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
    EMCGenericWriter<FReal> writer(name) ;
    writer.writeDataFromOctree(&tree,loader.getNumberOfParticles());
  }

  //------------------------------------------------------------------------------------------------------

  
  return 0;
}
