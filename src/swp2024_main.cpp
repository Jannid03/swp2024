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
    //std::ranges::sort (kmere);
    //seqan3::debug_stream << "Sortiert: " << *kmere_seq1.begin() << '\n';
    // auto kmere_seq1a = std::ranges::drop_view{kmere_seq1, 1};
    // seqan3::debug_stream << "Sortiert: " << kmere_seq1a << '\n';

    double i = jaccard_index(kmere_seq1, kmere_seq2);
    //double i = jaccard_index(std::vector<size_t> {1,2,2,3,5}, std::vector<size_t> {0,2,7,9});
    std::cout << i << std::endl;
    // auto mini = seq1 | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{3}}, 
    //                                                                 seqan3::window_size{5},
    //                                                                 seqan3::seed{0});
    // seqan3::debug_stream << mini << '\n';


    // //std::ranges::sort (mini);
    // seqan3::debug_stream << "Mini sortiert: " << *mini.begin() << '\n';
    // std::vector<int> sec;

    // std::ranges::set_intersection (kmere, mini, 
    //                         std::back_inserter(sec));
    // double sec_size = sec.size();

    // seqan3::debug_stream << sec_size << std::endl;

    /*std::vector<int> v1 {1,2,5,6};
    std::vector<int> v2 {1,3,4,5,9};

    std::vector<int> sec;

    std::set_intersection (v1.begin(), v1.end(), v2.begin(), v2.end(), 
                            std::back_inserter(sec));
    double sec_size = sec.size();

    sec.clear();
    std::set_union (v1.begin(), v1.end(), v2.begin(), v2.end(), 
                            std::back_inserter(sec));
    double uni_size = sec.size();

    seqan3::debug_stream << sec_size << std::endl;
    seqan3::debug_stream << uni_size << std::endl;
    seqan3::debug_stream << sec_size/uni_size << std::endl;*/
}