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
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>

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

    core::Real pert1 = numeric::random::gaussian();
    core::Real pert2 = numeric::random::gaussian();
    core::Size randres = static_cast< core::Size > ( numeric::random::uniform() * mypose->size() + 1 );
    core::Real orig_phi = mypose->phi( randres );
    core::Real orig_psi = mypose->psi( randres );
    mypose->set_phi( randres, orig_phi + pert1 );
    mypose->set_psi( randres, orig_psi + pert2 );

    protocols::moves::MonteCarlo mc(*mypose, *sfxn, 100);
    mc.boltzmann(*mypose);
    float mc_score = INFINITY;
    float next_mc_score = mc.last_accepted_score();
    for (int i = 0; i < 10; i++) {
        std::cout << next_mc_score << "\n";
        if (mc_score - next_mc_score < 0) {
            pert1 = numeric::random::gaussian();
            pert2 = numeric::random::gaussian();
            randres = static_cast< core::Size > ( numeric::random::uniform() * mypose->size() + 1 );

            mypose->set_phi( randres, orig_phi + pert1 );
            mypose->set_psi( randres, orig_psi + pert2 );
        }
        mc.boltzmann(*mypose);
        mc_score = next_mc_score;
        next_mc_score = mc.last_accepted_score();
    }
    return 0;
} 
