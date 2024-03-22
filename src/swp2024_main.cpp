#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> test {"AGTCGATGCTAGTCGAT"_dna4};
    seqan3::debug_stream << "Hello, World!" << std::endl;

    auto kmere = test | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << kmere << '\n';

    auto mini = test | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{3}}, 
                                                                    seqan3::window_size{5},
                                                                    seqan3::seed{0});
    seqan3::debug_stream << mini << '\n';
}