// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh> 
#include <core/pose/Pose.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.hh>           
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh> 
#include <core/pack/task/PackerTask.hh> 
#include <core/pack/task/TaskFactory.hh> 
#include <core/pack/pack_rotamers.hh>

#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>

// core/import/import_pose
int main( int argc, char ** argv ) {
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
        std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
    core::Real score = (sfxn)->score(*mypose);

    std::cout << score << " ";

    // i) MoveMap: Define degrees of freedom to optimize
    core::kinematics::MoveMap mm;
    mm.set_bb( true ); // Open backbone dihedral angles (phi/psi)
    mm.set_chi( true ); // Open side-chain dihedral angles (chi)

    // ii) MinimizerOptions: Define minimization parameters
    core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true /*use_nblist*/ );

    // iii) AtomTreeMinimizer: The minimization engine
    core::optimization::AtomTreeMinimizer atm;
    
    // iv) Copy Pose: Declare copy_pose outside for massive speedup
    core::pose::Pose copy_pose;

    core::Real pert1 = numeric::random::gaussian() * 10.0;
    core::Real pert2 = numeric::random::gaussian() * 10.0; 
    core::Size randres = static_cast< core::Size > ( numeric::random::uniform() * mypose->size() + 1 );
    core::Real orig_phi = mypose->phi( randres );
    core::Real orig_psi = mypose->psi( randres );
    mypose->set_phi( randres, orig_phi + pert1 );
    mypose->set_psi( randres, orig_psi + pert2 );


    // protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
    protocols::moves::MonteCarlo mc(*mypose, *sfxn, 100);
    mc.boltzmann(*mypose);
    for (int i = 0; i < 1000; i++) {
        if (i % 100 == 99) {
            mypose->dump_pdb("out" + std::to_string(i) + ".pdb");
        }
        // the_observer->pymol().apply( *mypose);       
        core::Size randres = static_cast< core::Size > ( numeric::random::uniform() * mypose->size() + 1 ); 
        core::Real pert_phi = numeric::random::gaussian() * 5.0;
        core::Real pert_psi = numeric::random::gaussian() * 5.0; 

        core::Real new_phi = mypose->phi( randres ) + pert_phi;
        core::Real new_psi = mypose->psi( randres ) + pert_psi;

        // 2. PACKING (Optimize side-chains on the fixed, perturbed backbone)
        // a. Create a new PackerTask object (one-time use)
        core::pack::task::PackerTaskOP repack_task =
            core::pack::task::TaskFactory::create_packer_task( *mypose );
        
        // b. Restrict task to repacking (no redesign)
        repack_task->restrict_to_repacking();
        
        // c. Run the packer
        core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

        
        // 3. MINIMIZATION (Optimize backbone and side-chains)
        // a. Copy current pose (optimization to avoid PyMOL slowdown)
        copy_pose = *mypose;
        
        // b. Run AtomTreeMinimizer on the copy
        atm.run( copy_pose, mm, *sfxn, min_opts );
        
        // c. Copy minimized pose back to mypose
        *mypose = copy_pose;

        mypose->set_phi( randres, new_phi );
        mypose->set_psi( randres, new_psi );
        
        mc.boltzmann(*mypose);
        
        core::Real current_score = mc.last_accepted_score();
        std::cout << current_score << "\n";

    }
    return 0;
} 
