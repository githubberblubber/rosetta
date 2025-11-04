// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  p/r/o/t/o/c/o/l/s///b/o/o/t/c/a/m/p/QueueTests.cxxtest.hh
/// @brief  queue testing
/// @author jonwu (jonwu@terpmail.umd.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/bootcamp/Queue.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("QueueTests");


class QueueTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
    }    

	void tearDown() {

	}

	void test_first() {
       TS_TRACE( "Running my first unit test!" );
       TS_ASSERT( true );

	}

    void test_enqueue() {
        queue = protocols::bootcamp::Queue();
        TS_ASSERT( queue.size() == 0 );
        TS_ASSERT( queue.is_empty() == true);
        queue.enqueue("abc");
        TS_ASSERT( queue.size() == 1);
        TS_ASSERT( queue.is_empty() == false);
    }

    void test_fifo_order() {
        queue = protocols::bootcamp::Queue();
        TS_TRACE( "Testing FIFO order of the queue." );
        std::string s1 = "First";
        std::string s2 = "Second";
        std::string s3 = "Third";

        queue.enqueue(s1);
        queue.enqueue(s2);
        queue.enqueue(s3);
        TS_ASSERT_EQUALS( queue.size(), 3 );

        TS_ASSERT_EQUALS( queue.dequeue(), s1 ); 
        TS_ASSERT_EQUALS( queue.size(), 2 );

        TS_ASSERT_EQUALS( queue.dequeue(), s2 ); 
        TS_ASSERT_EQUALS( queue.size(), 1 );

        TS_ASSERT_EQUALS( queue.dequeue(), s3 ); 
        TS_ASSERT_EQUALS( queue.size(), 0 );

        TS_ASSERT( queue.is_empty() );
    }

    void test_dequeue_to_empty() {
        queue = protocols::bootcamp::Queue();
        TS_TRACE( "Testing dequeue operation that fully empties the queue." );
        queue.enqueue("A");
        queue.enqueue("B");
        TS_ASSERT_EQUALS( queue.size(), 2 );

        queue.dequeue();
        TS_ASSERT_EQUALS( queue.size(), 1 );
        TS_ASSERT( !queue.is_empty() );

        queue.dequeue();
        TS_ASSERT_EQUALS( queue.size(), 0 );
        TS_ASSERT( queue.is_empty() );
    }

    void test_mixed_operations() {
        queue = protocols::bootcamp::Queue();
        TS_TRACE( "Testing mixed enqueue/dequeue operations." );
        
        queue.enqueue("X");
        TS_ASSERT_EQUALS( queue.size(), 1 );

        queue.enqueue("Y");
        TS_ASSERT_EQUALS( queue.size(), 2 );

        queue.dequeue();
        TS_ASSERT_EQUALS( queue.size(), 1 );
        TS_ASSERT_EQUALS( queue.dequeue(), "Y" ); 

        TS_ASSERT_EQUALS( queue.size(), 0 );
        
        queue.enqueue("Z");
        queue.enqueue("W");
        TS_ASSERT_EQUALS( queue.size(), 2 );

        TS_ASSERT_EQUALS( queue.dequeue(), "Z" );
    }

    void test_dequeue_from_empty() {
        queue = protocols::bootcamp::Queue();
        TS_TRACE( "Testing dequeue from an empty queue (expecting empty string or error)." );
        TS_ASSERT( queue.is_empty() );
        
        try {
            std::string result = queue.dequeue();
        } catch ( const std::exception& e ) {

        }

        queue = protocols::bootcamp::Queue();

        TS_ASSERT_EQUALS( queue.size(), 0 );
        TS_ASSERT( queue.is_empty() );
    }
private:
    protocols::bootcamp::Queue queue;

};
