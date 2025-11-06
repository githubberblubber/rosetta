// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/bootcamp/fold_tree_from_ss.cc
/// @brief  Declarations for utility functions to generate a FoldTree from secondary structure.
/// @author Jonathan Wu (jonathanw2856@gmail.com)

#include <protocols/bootcamp/fold_tree_from_ss.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>
#include <utility>

namespace protocols {
    namespace bootcamp {

    utility::vector1< std::pair< core::Size, core::Size > >
    identify_secondary_structure_spans( std::string const & ss_string )
    {
        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
        
        // Using core::Size(-1) as a safe sentinel value for 'uninitialized'
        core::Size strand_start = core::Size(-1); 
        
        // Loop is 0-indexed for ss_string
        for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
            bool is_ss_element = ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H' );

            if ( is_ss_element ) {
                // Safely check sentinel without casting to int
                if ( strand_start == core::Size(-1) ) { 
                    strand_start = ii;
                } else if ( ss_string[ii] != ss_string[strand_start] ) {
                    // End of the old SS element, start of a new one of different type.
                    // ii is the index of the first non-matching residue (0-indexed).
                    // ii is the 1-indexed residue number that ends the segment.
                    ss_boundaries.push_back( std::make_pair( strand_start + 1, ii ) );
                    strand_start = ii; 
                }
            } else { // Not an SS element
                // Safely check sentinel without casting to int
                if ( strand_start != core::Size(-1) ) { 
                    // End of the SS element span.
                    ss_boundaries.push_back( std::make_pair( strand_start + 1, ii ) );
                    strand_start = core::Size(-1); // Reset
                }
            }
        }

        // Handle the case where the SS element continues to the end of the string
        if ( strand_start != core::Size(-1) ) {
            ss_boundaries.push_back( std::make_pair( strand_start + 1, ss_string.size() ));
        }
        
        return ss_boundaries;
    }

    // =========================================================================
    // FUNCTION 1: fold_tree_from_ss (Wrapper function)
    // =========================================================================
    core::kinematics::FoldTree
    fold_tree_from_ss(core::pose::Pose const & mypose) {
        
        // 1. Calculate the secondary structure string using the DSSP scoring class
        core::scoring::dssp::Dssp dssp( mypose );  
        std::string ss_string = dssp.get_dssp_secstruct();
        
        // 2. Pass the string to the implementation function
        return fold_tree_from_dssp_string( ss_string );
    }

    // =========================================================================
    // FUNCTION 2: fold_tree_from_dssp_string (Implementation function)
    // =========================================================================
    core::kinematics::FoldTree
    fold_tree_from_dssp_string(std::string const & ss_string) {         
        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans( ss_string );

        core::Size const total_length = ss_string.length();
        core::kinematics::FoldTree new_ft;

        utility::vector1< core::Size > midpoints;
        for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
            core::Size const start = ss_boundaries[ ii ].first;
            core::Size const stop  = ss_boundaries[ ii ].second;
            midpoints.push_back(start + (stop - start) / 2);
        }
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
            
        }
        return new_ft;
    }

    } // namespace bootcamp
} // namespace protocols