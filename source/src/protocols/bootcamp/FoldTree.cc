#include <iostream>
#include <string>
#include <utility>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

void find_ones_block( utility::vector1<int> const & bitstring ) {
    int start = 0;
    for ( uint ii = 1; ii <= bitstring.size(); ++ii ) {
      std::cout << "ii: " << ii << " " << bitstring[ ii ] << std::endl;
      if ( start != 0 ) {
        if ( bitstring[ ii ] != 1 ) {
          std::cout << "start: " << start << " stop: " << ii-1 << std::endl;
          start = 0;
        }
      } else {
        if ( bitstring[ ii ] == 1 ) {
          start = ii;
        }
      }
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
    ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
  }
  for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
    std::cout << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}

int main( int argc, char *argv[] ) {

    // Check for correct number of arguments
    if ( argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " <function> <input_string>" << std::endl;
        std::cerr << "  function: 'ss_spans' or 'ones_block'" << std::endl;
        return 1;
    }

    // The first argument (argv[1]) is the function name
    std::string function_name = argv[1];
    
    // The second argument (argv[2]) is the input data string
    std::string input_data = argv[2];

    if ( function_name == "ss_spans" ) {
        std::cout << "--- Identifying Secondary Structure Spans ---" << std::endl;
        identify_secondary_structure_spans( input_data );

    } else if ( function_name == "ones_block" ) {
        utility::vector1< int > bitstring;
        for ( char c : input_data ) {
            if ( c == '0' ) {
                bitstring.push_back( 0 );
            } else if ( c == '1' ) {
                bitstring.push_back( 1 );
            } else {
                std::cerr << "Error: bitstring input must only contain '0' or '1'." << std::endl;
                return 1;
            }
        }
        
        std::cout << "--- Finding Ones Blocks ---" << std::endl;
        find_ones_block( bitstring );
        
    } else {
        std::cerr << "Error: Unknown function name '" << function_name << "'" << std::endl;
        return 1;
    }

    return 0;
}