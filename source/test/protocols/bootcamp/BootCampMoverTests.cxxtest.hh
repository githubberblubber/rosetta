// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>
#include <protocols/moves/MoverFactory.hh>

class BootCampMoverTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}
	
	void test_hello_world() {
		TS_ASSERT( true );
	}

    void test_bootcamp_mover_init() {
        // get_instance();
        protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover(protocols::bootcamp::BootCampMover::mover_name());
        utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( base_mover_op );
    }
};