#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>

#include <sharg/all.hpp>
#include <sharg/validators.hpp>

#include <algorithm>
#include <ranges>
#include <map>

using namespace seqan3::literals;


//TO DO: 1. to map conversion
// Wenn nicht -> System (Input1 +, Input2 -, dann verrechnen)

//2. Jaccard Funktion
//3. Ausgabe

struct eingabe {
    std::string modus{"k"};
    std::vector<std::filesystem::path> file{};
};

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

void intialize_parser (sharg::parser & parser, eingabe & in) {
    parser.info.author = "Jannik Dubrau";
    parser.info.short_description = "Jaccard Index für K-Mere und Minimizer";

    sharg::value_list_validator mode_validator{"m", "k", "b"};
    parser.add_option(in.modus, 
                        sharg::config{  .short_id = 'm', 
                                        .long_id = "mode", 
                                        .description ="Auswahl der Methode",
                                        .validator = mode_validator});
    int test {0};
    parser.add_option(test, sharg::config {.short_id = 't', .long_id = "Test", .hidden = true});
    /*sharg::input_file_validator my_file_ext_validator{{"fa", "fasta"}};
    parser.add_option(in.file, sharg::config{   .short_id = 'i',
                                                .long_id = "input",
                                                .description = "Bitte zwei Eingabedateien eingeben",
                                                .validator = my_file_ext_validator});*/
    
}

void kmere (std::vector<seqan3::dna4> & seq1, std::vector<seqan3::dna4> & seq2) {

    auto kmere_seq1 = seq1 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << kmere_seq1 << '\n';

    auto kmere_seq2 = seq2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << kmere_seq2 << '\n';

    double kmere_jac = jaccard_index(kmere_seq1, kmere_seq2);
    std::cout << kmere_jac << std::endl;
}

void mini (std::vector<seqan3::dna4> & seq1, std::vector<seqan3::dna4> & seq2) {
    uint64_t seed = 0x8F3F73B5CF1C9ADE;
    auto mini_seq1 = seq1 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}})
                            | std::views::transform([seed](uint64_t i) {return i^seed;})
                            | seqan3::views::minimiser (3);
    seqan3::debug_stream << mini_seq1 << '\n';

    auto mini_seq2 = seq2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}})
                            | std::views::transform([seed](uint64_t i) {return i^seed;})
                            | seqan3::views::minimiser (3);
    seqan3::debug_stream << mini_seq2 << '\n';

    double mini_jac = jaccard_index(mini_seq1, mini_seq2);
    std::cout << mini_jac << std::endl;
}

void run_program (char modus) {


    std::vector<seqan3::dna4> seq1 {"AGCTGTCGAAAGTCGAAAT"_dna4};
    std::vector<seqan3::dna4> seq2 {"CATGATGTCACTGATCGTA"_dna4};

    if (modus == 'k') {

        kmere(seq1, seq2);
    }
    else if (modus == 'm') {

        mini(seq1, seq2);
    }
    else {
        kmere(seq1,seq2);
        mini(seq1,seq2);
    }
    
}

int main(int argc, char** argv) {

    sharg::parser myparser{"swp2024_main", argc, argv};
    eingabe in{};
    intialize_parser(myparser, in);

    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << ext.what() << "\n"; // customise your error message
        return -1;
    }
 
    run_program(in.modus[0]);

}