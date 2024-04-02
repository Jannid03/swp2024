#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>

#include <algorithm>
#include <ranges>
#include <map>

using namespace seqan3::literals;


//TO DO: 1. to map conversion
// Wenn nicht -> System (Input1 +, Input2 -, dann verrechnen)

//2. Jaccard Funktion
//3. Ausgabe
double jaccard_index (auto first, auto second) {
    double size_1 {0};
    double size_2 {0};
    double uni {0};
    std::map<size_t, size_t> dic;

    for (auto i : first) {
        auto search = dic.find(i);
        size_1++;

        if (search == dic.end()) {
            dic.insert({i, 1});
        }
        else {
            dic.at(i)++;
        }
    }

    for (auto i : second) {
        auto search = dic.find(i);

        if(search == dic.end()) {
            size_2++;
        }
        // else if (dic.at(i) == 0) {
        //     size_2++;
        // }
        else {
            dic.at(i)--;
            uni++;

            if (dic.at(i) == 0) {
                dic.erase(i);
            }
        }

    }

    return uni/(size_1+size_2);
}

int main()
{
    std::vector<seqan3::dna4> seq1 {"AGCTGTCGAAAGTCGAAAT"_dna4};
    std::vector<seqan3::dna4> seq2 {"CATGATGTCACTGATCGTA"_dna4};

    auto kmere_seq1 = seq1 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << kmere_seq1 << '\n';

    auto kmere_seq2 = seq2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << kmere_seq2 << '\n';

    double kmere_jac = jaccard_index(kmere_seq1, kmere_seq2);
    std::cout << kmere_jac << std::endl;

    auto mini_seq1 = seq1 | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{3}}, 
                                                                    seqan3::window_size{5},
                                                                    seqan3::seed{0});
    seqan3::debug_stream << mini_seq1 << '\n';

    auto mini_seq2 = seq2 | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{3}}, 
                                                                    seqan3::window_size{5},
                                                                    seqan3::seed{0});
    seqan3::debug_stream << mini_seq2 << '\n';

    double mini_jac = jaccard_index(mini_seq1, mini_seq2);
    std::cout << mini_jac << std::endl;

}