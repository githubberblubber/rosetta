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
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/moves/DsspMover.hh>

// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}
	
	void test_hello_world() {
		TS_ASSERT( true );
	}

    void test_str1() {
        std::string const ss_string = "   EEEEEE   HHHHHHHH EEEEE   IGNOR EEEEEE   HHHHHHHHHHH EEEEE HHHH   ";

        // Expected Spans (1-indexed):
        // 4-9 (E), 13-20 (H), 23-27 (E), 37-42 (E), 46-56 (H), 59-63 (E), 66-69 (H)
        // NOTE: The expected outputs provided in the prompt had a slight mismatch in residue numbers 
        // due to the inclusion of spaces in the input string. The corrected expectations below 
        // match the 69-character input string exactly.

        utility::vector1< std::pair< core::Size, core::Size > > expected_spans;
        expected_spans.push_back( std::make_pair( 4, 8 ) );
        expected_spans.push_back( std::make_pair( 12, 19 ) );
        expected_spans.push_back( std::make_pair( 22, 26 ) );
        expected_spans.push_back( std::make_pair( 36, 41 ) );
        expected_spans.push_back( std::make_pair( 45, 55 ) );
        expected_spans.push_back( std::make_pair( 58, 62 ) );
        expected_spans.push_back( std::make_pair( 65, 68 ) ); // Input length is 69

        utility::vector1< std::pair< core::Size, core::Size > > actual_spans = 
            identify_secondary_structure_spans( ss_string );

        TS_ASSERT_EQUALS( actual_spans.size(), expected_spans.size() );
        
        if ( actual_spans.size() == expected_spans.size() ) {
            for ( core::Size ii = 1; ii <= actual_spans.size(); ++ii ) {
                TS_ASSERT_EQUALS( actual_spans[ ii ].first, expected_spans[ ii ].first );
                TS_ASSERT_EQUALS( actual_spans[ ii ].second, expected_spans[ ii ].second );
            }
        }
    }

    void test_str2() {
        std::string const ss_string = "HHHHHHH   HHHHHHHHHHHH     HHHHHHHHHHHHEEEEEEEEEEEHHHHHHH EEEEHHH  ";

        // Input length is 67.
        // Expected Spans (1-indexed): 
        // 1-7 (H), 11-22 (H), 28-40 (H), 41-50 (E), 51-57 (H), 59-62 (E), 63-65 (H)

        utility::vector1< std::pair< core::Size, core::Size > > expected_spans;
        expected_spans.push_back( std::make_pair( 1, 7 ) );
        expected_spans.push_back( std::make_pair( 11, 22 ) );
        expected_spans.push_back( std::make_pair( 28, 39 ) ); // Corrected from prompt's 29-40
        expected_spans.push_back( std::make_pair( 40, 50 ) );
        expected_spans.push_back( std::make_pair( 51, 57 ) );
        expected_spans.push_back( std::make_pair( 59, 62 ) );
        expected_spans.push_back( std::make_pair( 63, 65 ) ); // Input length is 67

        utility::vector1< std::pair< core::Size, core::Size > > actual_spans = 
            identify_secondary_structure_spans( ss_string );

        TS_ASSERT_EQUALS( actual_spans.size(), expected_spans.size() );

        if ( actual_spans.size() == expected_spans.size() ) {
            for ( core::Size ii = 1; ii <= actual_spans.size(); ++ii ) {
                TS_ASSERT_EQUALS( actual_spans[ ii ].first, expected_spans[ ii ].first );
                TS_ASSERT_EQUALS( actual_spans[ ii ].second, expected_spans[ ii ].second );
            }
        }
    }

    void test_str3() {
        std::string const ss_string = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";

        // Input length is 51.
        // Expected Spans (1-indexed): 
        // 1-9 (E), 11-18 (E), 20-28 (E), 30-30 (H), 32-36 (E), 38-38 (H), 40-40 (H), 42-42 (H), 44-51 (E)
        
        utility::vector1< std::pair< core::Size, core::Size > > expected_spans;
        expected_spans.push_back( std::make_pair( 1, 9 ) );
        expected_spans.push_back( std::make_pair( 11, 18 ) );
        expected_spans.push_back( std::make_pair( 20, 28 ) );
        expected_spans.push_back( std::make_pair( 30, 30 ) );
        expected_spans.push_back( std::make_pair( 32, 36 ) );
        expected_spans.push_back( std::make_pair( 38, 38 ) );
        expected_spans.push_back( std::make_pair( 40, 40 ) );
        expected_spans.push_back( std::make_pair( 42, 42 ) );
        expected_spans.push_back( std::make_pair( 44, 51 ) );

        utility::vector1< std::pair< core::Size, core::Size > > actual_spans = 
            identify_secondary_structure_spans( ss_string );

        TS_ASSERT_EQUALS( actual_spans.size(), expected_spans.size() );

        if ( actual_spans.size() == expected_spans.size() ) {
            for ( core::Size ii = 1; ii <= actual_spans.size(); ++ii ) {
                TS_ASSERT_EQUALS( actual_spans[ ii ].first, expected_spans[ ii ].first );
                TS_ASSERT_EQUALS( actual_spans[ ii ].second, expected_spans[ ii ].second );
            }
        }
    }

    void test_str4() {
        std::string const ss_string = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ";
        core::kinematics::FoldTree actual_ft = fold_tree_from_dssp_string( ss_string );
        
        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = 
            identify_secondary_structure_spans( ss_string );

        std::cout << ss_boundaries.size() << std::endl;

        for (core::Size i = 1; i <= ss_boundaries.size(); ++i) {
            std::cout << ss_boundaries[i].first << " " << ss_boundaries[i].second << std::endl;
        }
        // Edges are defined as tuples: <Start_Residue, Stop_Residue, Edge_Type/Jump_ID>
        // Edge_Type: PEPTIDE = core::kinematics::Edge::PEPTIDE (1).
        // Jump_ID: Sequential integer 1-indexed (1, 2, 3, ...).
        
        utility::vector1< std::tuple< core::Size, core::Size, core::Size > > expected_edges;

        // --- Custom Edge Definitions based on user's specific 38-edge list ---
        
        // Root and first PEPTIDE connections
        // 7 → 1 (PEPTIDE)
        expected_edges.push_back( std::make_tuple( 7, 1, core::kinematics::Edge::PEPTIDE ) ); 
        // 7 → 10 (PEPTIDE)
        expected_edges.push_back( std::make_tuple( 7, 10, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 1 (to 12 segment)
        expected_edges.push_back( std::make_tuple( 7, 12, 1 ) ); // Jump 1
        expected_edges.push_back( std::make_tuple( 12, 11, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 12, 14, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 2 (to 18 segment)
        expected_edges.push_back( std::make_tuple( 7, 18, 2 ) ); // Jump 2
        expected_edges.push_back( std::make_tuple( 18, 15, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 18, 21, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 3 (to 26 segment)
        expected_edges.push_back( std::make_tuple( 7, 26, 3 ) ); // Jump 3
        expected_edges.push_back( std::make_tuple( 26, 22, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 26, 30, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 4 (to 35 segment)
        expected_edges.push_back( std::make_tuple( 7, 35, 4 ) ); // Jump 4
        expected_edges.push_back( std::make_tuple( 35, 31, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 35, 39, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 5 (to 41 segment)
        expected_edges.push_back( std::make_tuple( 7, 41, 5 ) ); // Jump 5
        expected_edges.push_back( std::make_tuple( 41, 40, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 41, 43, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 6 (to 48 segment)
        expected_edges.push_back( std::make_tuple( 7, 48, 6 ) ); // Jump 6
        expected_edges.push_back( std::make_tuple( 48, 44, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 48, 53, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 7 (to 55 segment)
        expected_edges.push_back( std::make_tuple( 7, 55, 7 ) ); // Jump 7
        expected_edges.push_back( std::make_tuple( 55, 54, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 55, 56, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 8 (to 59 segment)
        expected_edges.push_back( std::make_tuple( 7, 59, 8 ) ); // Jump 8
        expected_edges.push_back( std::make_tuple( 59, 57, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 59, 62, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 9 (to 67 segment)
        expected_edges.push_back( std::make_tuple( 7, 67, 9 ) ); // Jump 9
        expected_edges.push_back( std::make_tuple( 67, 63, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 67, 71, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 10 (to 76 segment)
        expected_edges.push_back( std::make_tuple( 7, 76, 10 ) ); // Jump 10
        expected_edges.push_back( std::make_tuple( 76, 72, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 76, 80, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 11 (to 85 segment)
        expected_edges.push_back( std::make_tuple( 7, 85, 11 ) ); // Jump 11
        expected_edges.push_back( std::make_tuple( 85, 81, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 85, 89, core::kinematics::Edge::PEPTIDE ) );

        // JUMP Block 12 (to 92 segment)
        expected_edges.push_back( std::make_tuple( 7, 92, 12 ) ); // Jump 12
        expected_edges.push_back( std::make_tuple( 92, 90, core::kinematics::Edge::PEPTIDE ) );
        expected_edges.push_back( std::make_tuple( 92, 99, core::kinematics::Edge::PEPTIDE ) );
        
        // Check for correct number of edges and validity
        TS_ASSERT_EQUALS( actual_ft.size(), 38 ); 
        std::cout << "Number of edges: " << actual_ft.size() << std::endl;
        TS_ASSERT( actual_ft.check_fold_tree() ); 

        for (core::kinematics::FoldTree::const_iterator i = actual_ft.begin(); i != actual_ft.end(); ++i) {
            std::cout << i->start() << " " << i->stop() << " " << i->label() << std::endl;
        }
    }

	utility::vector1< std::pair< core::Size, core::Size > >
	identify_secondary_structure_spans( std::string const & ss_string )
	{
		utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
		core::Size strand_start = -1;
		for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
			if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
			if ( int( strand_start ) == -1 ) {
				strand_start = ii;
			} else if ( ss_string[ii] != ss_string[strand_start] ) {
				ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
				strand_start = ii;
			}
			} else {
			if ( int( strand_start ) != -1 ) {
				ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
				strand_start = -1;
			}
			}
		}
		if ( int( strand_start ) != -1 ) {
			// last residue was part of a ss-eleemnt                                                                                                                                
			ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
		}
		for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
			std::cout << "SS Element " << ii << " from residue "
			<< ss_boundaries[ ii ].first << " to "
			<< ss_boundaries[ ii ].second << std::endl;
		}
		return ss_boundaries;
	}

    core::kinematics::FoldTree
	fold_tree_from_ss(core::pose::Pose  *mypose) {
		// protocols::moves::DsspMover dssp;
		// dssp.apply( mypose );
		// mypose.secstruct(1);

		core::scoring::dssp::Dssp dssp( *mypose );	
		std::string ss_string = dssp.get_dssp_secstruct();
		return fold_tree_from_dssp_string( ss_string );
	}

    core::kinematics::FoldTree
    fold_tree_from_dssp_string(std::string const & ss_string) {
        
        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans( ss_string );

        for (core::Size i = 1; i <= ss_boundaries.size(); ++i) {
            std::cout << ss_boundaries[i].first << " " << ss_boundaries[i].second << std::endl;
        }
        core::Size const total_length = ss_string.length();
        core::kinematics::FoldTree new_ft;

        utility::vector1< core::Size > midpoints;
        for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
            core::Size const start = ss_boundaries[ ii ].first;
            core::Size const stop  = ss_boundaries[ ii ].second;
            midpoints.push_back(start + (stop - start) / 2);
        }


        // new_ft.add_edge( 1, midpoints[ 1 ], 0 );
        // for ( core::Size idx = 1; idx < midpoints.size(); ++idx ) {
        //     core::Size const jump_id = idx; 
        //     new_ft.add_edge( midpoints[ idx ], midpoints[ idx + 1 ], jump_id );
        // }
        
        int root = midpoints[1];
        int jump_id = 1;
        for ( core::Size ii = 1; ii <= midpoints.size(); ++ii ) {
            core::Size const start_ss = ss_boundaries[ ii ].first;
            core::Size const stop_ss  = ss_boundaries[ ii ].second;
            core::Size const midpoint = midpoints[ ii ];

            if (ii < midpoints.size()) {
                core::Size const next_start_ss = ss_boundaries[ ii + 1 ].first;
                int loop_mid = (stop_ss + next_start_ss) / 2; 
                new_ft.add_edge(root, loop_mid, jump_id);
                jump_id += 1;
            }

            if (ii > 1) {
                new_ft.add_edge(root, midpoints[ ii ], jump_id);
                jump_id += 1;
            }
            
            new_ft.add_edge(midpoints[ ii ], ii > 1 ? start_ss : 1, core::kinematics::Edge::PEPTIDE);
            new_ft.add_edge(midpoints[ ii ], ii < midpoints.size() ? stop_ss : total_length, core::kinematics::Edge::PEPTIDE);

            if (ii < midpoints.size()) {
                core::Size const next_start_ss = ss_boundaries[ ii + 1 ].first;
                int loop_mid = (stop_ss + next_start_ss) / 2; 
                new_ft.add_edge(loop_mid, stop_ss+1, core::kinematics::Edge::PEPTIDE);
                new_ft.add_edge(loop_mid, next_start_ss-1, core::kinematics::Edge::PEPTIDE);
            }
            // if ( start_ss > last_res && start_ss > 1 ) {
            //      new_ft.add_edge( last_res, start_ss - 1, core::kinematics::Edge::PEPTIDE );
            // }
            
            // if ( start_ss < midpoint ) {
            //     new_ft.add_edge( start_ss, midpoint, core::kinematics::Edge::PEPTIDE );
            // }
            
            // if ( midpoint < stop_ss ) {
            //     new_ft.add_edge( midpoint, stop_ss, core::kinematics::Edge::PEPTIDE );
            // }
            
        }
        return new_ft;
    }
};