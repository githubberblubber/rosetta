// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/bootcamp/fold_tree_from_ss.hh
/// @brief  Declarations for utility functions to generate a FoldTree from secondary structure.
/// @author Jonathan Wu (jonathanw2856@gmail.com)

#ifndef INCLUDED_protocols_bootcamp_fold_tree_from_ss_hh
#define INCLUDED_protocols_bootcamp_fold_tree_from_ss_hh

#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <utility>

namespace protocols {
    namespace bootcamp {

    /// @brief Identifies continuous spans of secondary structure elements ('E' or 'H')
    /// from a DSSP string and returns them as 1-indexed (start, stop) residue pairs.
    utility::vector1< std::pair< core::Size, core::Size > >
    identify_secondary_structure_spans( std::string const & ss_string );

    /// @brief Generates a FoldTree for a Pose by calculating secondary structure via DSSP.
    // Uses const reference as per standard Rosetta practice.
    core::kinematics::FoldTree
    fold_tree_from_ss( core::pose::Pose const & mypose ); 

    /// @brief Constructs a FoldTree based on the midpoints of secondary structure spans
    /// defined in the input DSSP string.
    core::kinematics::FoldTree
    fold_tree_from_dssp_string( std::string const & ss_string );

    } // namespace bootcamp
} // namespace protocols

#endif // INCLUDED_protocols_bootcamp_fold_tree_from_ss_hh