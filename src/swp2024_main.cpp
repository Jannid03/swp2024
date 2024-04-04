#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>

#include <sharg/all.hpp>
#include <sharg/validators.hpp>

#include <algorithm>
#include <ranges>
#include <map>
#include <limits.h>

using namespace seqan3::literals;


//TO DO: 
// 1. K mer Customizing
// 2. File input 
// 3. Window size custom

//Struct für die Eingabeparameter
struct eingabe {
    std::string modus{"k"};
    uint8_t k {5};
    size_t window {5};
    std::vector<std::filesystem::path> file{};
};

//Generelle Idee: Wir gehen die kmere der ersten Sequenz durch ujd merken uns, 
// welche wie oft vorkommen (map). Wir merken uns auch die Anzahl (size_1)
//Anschließend gehen wir den zweiten Vektor durch und suchen, ob es Vorkommen gibt
// Wenn ja -> -1 in Map und Anzahl der geteilten kmere erhöht (inter)
// Wenn nein -> nur size_2 erhöht, size_2 steht für exklusiven Teil der zweiten Sequenz 
// Also: (#Kmere in Sequenz 2 / #Kmere in Intersection von 1 und 2)
// inter -> #Kmere in Intersection
// size_1 -> #Kmere in Sequenz 1 (INKLUSIVE der kmere in der Intersection)
// size_2 -> #Kmere in Sequenz EXKLUSIVE der kmere in der Intersection
// Darauf folgt: size_1 + size_2 == Union von Seq 1 und Seq 2
// Anschließend Jaccard Berechnung: inter / (size_1 + size_2)
double jaccard_index (auto first, auto second) {
    //Intialisierung
    double size_1 {0};
    double size_2 {0};
    double inter {0};
    std::map<size_t, size_t> dic;

    //Durch erste Sequenz, erstellen des dictionary
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

    //Zweite Sequenz: Wenn nicht gefunden, size_2 um eins erhöht
    for (auto i : second) {
        auto search = dic.find(i);

        if(search == dic.end()) {
            size_2++;
        }
        //Wenn gefunden, dann inter erhöht und in Dcitionary um eins verringert
        else {
            dic.at(i)--;
            inter++;

            //wenn aufgebraucht, dann aus dictionary entfernt
            if (dic.at(i) == 0) {
                dic.erase(i);
            }
        }

    }

    //Jaccard Index berechnet und zurückgegeben
    return inter/(size_1+size_2);
}

void intialize_parser (sharg::parser & parser, eingabe & in) {
    parser.info.author = "Jannik Dubrau";
    parser.info.short_description = "Jaccard Index für K-Mere und Minimizer";
    parser.info.version = "1.0";
    
    //Einfügen der Optionene für den Parser
    sharg::value_list_validator mode_validator{"m", "k", "b"};
    parser.add_option(in.modus, 
                        sharg::config{  .short_id = 'm', 
                                        .long_id = "mode", 
                                        .description ="Auswahl der Methode",
                                        .validator = mode_validator});

    sharg::arithmetic_range_validator kmere_validator{1, 58};
    parser.add_option(in.k, sharg::config {.short_id = 'k', .long_id = "kmer",
                                            .description = "Länge der Kmere",
                                            .required = true,
                                            .validator = kmere_validator});
    
    parser.add_option(in.window, sharg::config {.short_id = 'w', .long_id = "window",
                                            .description = "Länge des Windows"});
    /*sharg::input_file_validator my_file_ext_validator{{"fa", "fasta"}};
    parser.add_option(in.file, sharg::config{   .short_id = 'i',
                                                .long_id = "input",
                                                .description = "Bitte zwei Eingabedateien eingeben",
                                                .validator = my_file_ext_validator});*/
    
}

void kmere (std::vector<seqan3::dna4> & seq1, std::vector<seqan3::dna4> & seq2, uint8_t km) {
    //Aufruf der Funktionene und der Berechnung
    auto kmere_seq1 = seq1 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{km}});
    seqan3::debug_stream << kmere_seq1 << '\n';

    auto kmere_seq2 = seq2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{km}});
    seqan3::debug_stream << kmere_seq2 << '\n';

    double kmere_jac = jaccard_index(kmere_seq1, kmere_seq2);
    std::cout << kmere_jac << std::endl;
}

void mini (std::vector<seqan3::dna4> & seq1, std::vector<seqan3::dna4> & seq2, uint8_t km, uint32_t w) {
    uint64_t seed = 0x8F3F73B5CF1C9ADE; //Höchster Seed, eventuell Anpassugn später
    uint32_t window = w - km + 1; //window size Berechnung für eingabe bei minimiser
    auto mini_seq1 = seq1 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{km}})
                            | std::views::transform([seed](uint64_t i) {return i^seed;})
                            | seqan3::views::minimiser (window);
    seqan3::debug_stream << mini_seq1 << '\n';

    auto mini_seq2 = seq2 | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{km}})
                            | std::views::transform([seed](uint64_t i) {return i^seed;})
                            | seqan3::views::minimiser (window);
    seqan3::debug_stream << mini_seq2 << '\n';

    //Aufruf der Berechnung
    double mini_jac = jaccard_index(mini_seq1, mini_seq2);
    std::cout << mini_jac << std::endl;
}

void run_program (eingabe & in) {
    //Modus wird abgelesen, da string eingabe aber char einfacher zu vergleichen -> modus[0]
    char modus = in.modus[0];
    std::vector<seqan3::dna4> seq1 {"AGCTGTCGAAAGTCGAAAT"_dna4};
    std::vector<seqan3::dna4> seq2 {"CATGATGTCACTGATCGTA"_dna4};

    //Modus wird abgefragt mit entsprechenden Aufrufen
    if (modus == 'k') {

        kmere(seq1, seq2, in.k);
    }
    else if (modus == 'm') {

        mini(seq1, seq2, in.k, in.window);
    }
    else {
        kmere(seq1,seq2, in.k);
        mini(seq1,seq2, in.k, in.window);
    }
    
}

int main(int argc, char** argv) {
    //Erstellen und Intialisieren des Parsers
    sharg::parser myparser{"swp2024_main", argc, argv};
    eingabe in{};
    intialize_parser(myparser, in);

    //Block aus Sharg tutorial übernommen
    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << ext.what() << "\n"; // customise your error message
        return -1;
    }
 
    //Tatsächliche Programmausführung, eingegebene Argumente werden übergeben.
    run_program(in);

}