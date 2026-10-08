#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
namespace fs = std::filesystem;


// ============================================================
// 去除字符串两端空格
// ============================================================
string trim(const string& s) {

    size_t start = 0;

    while (start < s.size() &&
           isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }

    size_t end = s.size();

    while (end > start &&
           isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }

    return s.substr(start, end - start);
}


// ============================================================
// CSV解析
// 支持：
// abc,def,"hello,world",123
// ============================================================
vector<string> parseCSVLine(const string& line) {

    vector<string> fields;
    string current;

    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {

        char c = line[i];

        if (c == '"') {

            if (in_quotes &&
                i + 1 < line.size() &&
                line[i + 1] == '"') {

                current += '"';
                ++i;

            } else {

                in_quotes = !in_quotes;
            }

        } else if (c == ',' && !in_quotes) {

            fields.push_back(current);
            current.clear();

        } else {

            current += c;
        }
    }

    fields.push_back(current);

    return fields;
}


// ============================================================
// CSV字段转义
// ============================================================
string escapeCSVField(const string& field) {

    bool need_quotes = false;

    for (char c : field) {

        if (c == ',' ||
            c == '"' ||
            c == '\n' ||
            c == '\r') {

            need_quotes = true;
            break;
        }
    }

    if (!need_quotes) {
        return field;
    }

    string result = "\"";

    for (char c : field) {

        if (c == '"') {
            result += "\"\"";
        } else {
            result += c;
        }
    }

    result += "\"";

    return result;
}


// ============================================================
// 写CSV行
// ============================================================
void writeCSVRow(
    ofstream& out,
    const vector<string>& fields
) {

    for (size_t i = 0; i < fields.size(); ++i) {

        if (i > 0) {
            out << ",";
        }

        out << escapeCSVField(fields[i]);
    }

    out << "\n";
}


// ============================================================
// FASTA每行固定长度
// ============================================================
void writeFastaSequence(
    ofstream& out,
    const string& sequence,
    size_t width = 80
) {

    for (size_t i = 0; i < sequence.size(); i += width) {

        out << sequence.substr(i, width)
            << "\n";
    }
}


// ============================================================
// 主函数
// ============================================================
int main() {

    // ========================================================
    // 输入文件1：
    // 病毒分类信息
    // ========================================================

    const string virus_file =
        "Viro3D_viruses_list_with_coverage.csv";


    // ========================================================
    // 输入文件2：
    // 已经筛选：
    // protein
    // pLDDT >= 80
    // 50 <= length <= 2000
    // ========================================================

    const string protein_file =
        "Viro3D_proteins_filtered.csv";


    // ========================================================
    // 目标Family
    // ========================================================

    const string TARGET_FAMILY =
        "Baculoviridae";


    // ========================================================
    // 输出目录
    // ========================================================

    const string output_dir =
        "processed/Baculoviridae";


    fs::create_directories(output_dir);


    const string output_csv =
        output_dir +
        "/Baculoviridae_proteins.csv";


    const string output_fasta =
        output_dir +
        "/Baculoviridae_proteins.fasta";


    const string output_summary =
        output_dir +
        "/Baculoviridae_summary.txt";


    // ========================================================
    // Part 1
    // 读取virus metadata
    //
    // 建立：
    // ICTV Sort -> ICTV Family
    // ========================================================

    ifstream virus_in(virus_file);

    if (!virus_in.is_open()) {

        cerr
            << "ERROR: cannot open virus metadata:\n"
            << virus_file
            << endl;

        return 1;
    }


    string line;

    if (!getline(virus_in, line)) {

        cerr
            << "ERROR: virus metadata is empty."
            << endl;

        return 1;
    }


    if (!line.empty() &&
        line.back() == '\r') {

        line.pop_back();
    }


    vector<string> virus_headers =
        parseCSVLine(line);


    unordered_map<string, size_t>
        virus_col;


    for (size_t i = 0;
         i < virus_headers.size();
         ++i) {

        virus_col[
            trim(virus_headers[i])
        ] = i;
    }


    if (
        virus_col.find("ICTV Sort")
        ==
        virus_col.end()
        ||
        virus_col.find("ICTV Family")
        ==
        virus_col.end()
    ) {

        cerr
            << "ERROR: ICTV Sort or ICTV Family column missing."
            << endl;

        return 1;
    }


    size_t idx_v_sort =
        virus_col["ICTV Sort"];

    size_t idx_v_family =
        virus_col["ICTV Family"];


    unordered_map<string, string>
        sort_to_family;


    while (getline(virus_in, line)) {

        if (!line.empty() &&
            line.back() == '\r') {

            line.pop_back();
        }


        if (line.empty()) {
            continue;
        }


        vector<string> fields =
            parseCSVLine(line);


        if (
            fields.size()
            <
            virus_headers.size()
        ) {

            continue;
        }


        string sort_id =
            trim(fields[idx_v_sort]);

        string family =
            trim(fields[idx_v_family]);


        if (!sort_id.empty()) {

            sort_to_family[
                sort_id
            ] = family;
        }
    }


    virus_in.close();


    cout
        << "Loaded ICTV Sort -> Family mappings: "
        << sort_to_family.size()
        << endl;


    // ========================================================
    // Part 2
    // 读取筛选后的protein文件
    // ========================================================

    ifstream protein_in(protein_file);

    if (!protein_in.is_open()) {

        cerr
            << "ERROR: cannot open protein file:\n"
            << protein_file
            << endl;

        return 1;
    }


    if (!getline(protein_in, line)) {

        cerr
            << "ERROR: protein file is empty."
            << endl;

        return 1;
    }


    if (!line.empty() &&
        line.back() == '\r') {

        line.pop_back();
    }


    vector<string> protein_headers =
        parseCSVLine(line);


    unordered_map<string, size_t>
        protein_col;


    for (size_t i = 0;
         i < protein_headers.size();
         ++i) {

        protein_col[
            trim(protein_headers[i])
        ] = i;
    }


    // ========================================================
    // 需要的列
    // ========================================================

    vector<string> required_columns = {

        "ICTV Sort",
        "ICTV Species",
        "Viro3D ID",
        "GenBank Protein ID",
        "GenBank Genome Accession",
        "Protein Sequence"
    };


    for (const string& col :
         required_columns) {

        if (
            protein_col.find(col)
            ==
            protein_col.end()
        ) {

            cerr
                << "ERROR: missing protein column: "
                << col
                << endl;

            return 1;
        }
    }


    size_t idx_sort =
        protein_col["ICTV Sort"];

    size_t idx_species =
        protein_col["ICTV Species"];

    size_t idx_viro3d_id =
        protein_col["Viro3D ID"];

    size_t idx_genbank_protein =
        protein_col["GenBank Protein ID"];

    size_t idx_genome =
        protein_col["GenBank Genome Accession"];

    size_t idx_sequence =
        protein_col["Protein Sequence"];


    // ========================================================
    // 打开输出
    // ========================================================

    ofstream csv_out(output_csv);

    ofstream fasta_out(output_fasta);


    if (!csv_out.is_open() ||
        !fasta_out.is_open()) {

        cerr
            << "ERROR: cannot create output files."
            << endl;

        return 1;
    }


    // CSV原始字段全部保留
    // 并额外增加 ICTV Family
    vector<string> new_headers =
        protein_headers;

    new_headers.push_back(
        "ICTV Family"
    );


    writeCSVRow(
        csv_out,
        new_headers
    );


    // ========================================================
    // 统计
    // ========================================================

    size_t total_input_proteins = 0;

    size_t selected_proteins = 0;

    size_t unmatched_family = 0;


    unordered_set<string>
        species_set;

    unordered_set<string>
        genome_set;

    unordered_set<string>
        viro3d_id_set;


    // ========================================================
    // 遍历protein
    // ========================================================

    while (getline(protein_in, line)) {

        if (!line.empty() &&
            line.back() == '\r') {

            line.pop_back();
        }


        if (line.empty()) {
            continue;
        }


        vector<string> fields =
            parseCSVLine(line);


        if (
            fields.size()
            <
            protein_headers.size()
        ) {

            cerr
                << "WARNING: malformed protein row skipped."
                << endl;

            continue;
        }


        total_input_proteins++;


        string sort_id =
            trim(fields[idx_sort]);


        auto it =
            sort_to_family.find(
                sort_id
            );


        if (
            it ==
            sort_to_family.end()
        ) {

            unmatched_family++;
            continue;
        }


        const string& family =
            it->second;


        // ====================================================
        // 只留下Baculoviridae
        // ====================================================

        if (
            family
            !=
            TARGET_FAMILY
        ) {

            continue;
        }


        selected_proteins++;


        string species =
            trim(
                fields[idx_species]
            );

        string genome =
            trim(
                fields[idx_genome]
            );

        string viro3d_id =
            trim(
                fields[idx_viro3d_id]
            );

        string genbank_protein =
            trim(
                fields[idx_genbank_protein]
            );

        string sequence =
            trim(
                fields[idx_sequence]
            );


        species_set.insert(
            species
        );

        genome_set.insert(
            genome
        );

        viro3d_id_set.insert(
            viro3d_id
        );


        // ====================================================
        // 输出CSV
        // ====================================================

        vector<string> output_fields =
            fields;

        output_fields.push_back(
            family
        );


        writeCSVRow(
            csv_out,
            output_fields
        );


        // ====================================================
        // 输出FASTA
        //
        // header:
        // >Viro3D_ID|GenBank_ID|Species
        // ====================================================

        fasta_out
            << ">"
            << viro3d_id
            << "|"
            << genbank_protein
            << "|"
            << species
            << "\n";


        writeFastaSequence(
            fasta_out,
            sequence
        );
    }


    protein_in.close();

    csv_out.close();

    fasta_out.close();


    // ========================================================
    // summary
    // ========================================================

    ofstream summary_out(
        output_summary
    );


    summary_out
        << "Baculoviridae subset summary\n";

    summary_out
        << "========================================\n\n";


    summary_out
        << "Source filtered protein file:\n"
        << protein_file
        << "\n\n";


    summary_out
        << "Target ICTV Family:\n"
        << TARGET_FAMILY
        << "\n\n";


    summary_out
        << "Total proteins in filtered Viro3D dataset: "
        << total_input_proteins
        << "\n";


    summary_out
        << "Selected Baculoviridae proteins: "
        << selected_proteins
        << "\n";


    summary_out
        << "Unique ICTV species: "
        << species_set.size()
        << "\n";


    summary_out
        << "Unique genome accessions: "
        << genome_set.size()
        << "\n";


    summary_out
        << "Unique Viro3D IDs: "
        << viro3d_id_set.size()
        << "\n";


    summary_out
        << "Proteins without Family mapping: "
        << unmatched_family
        << "\n";


    summary_out
        << "\nFiltering inherited from previous step:\n";

    summary_out
        << "  Peptide Category == protein\n";

    summary_out
        << "  ColabFold pLDDT >= 80\n";

    summary_out
        << "  Protein Length >= 50\n";

    summary_out
        << "  Protein Length <= 2000\n";


    summary_out.close();


    // ========================================================
    // 终端结果
    // ========================================================

    cout
        << endl
        << "========================================"
        << endl;

    cout
        << "Baculoviridae extraction summary"
        << endl;

    cout
        << "========================================"
        << endl;


    cout
        << "Input filtered proteins:       "
        << total_input_proteins
        << endl;

    cout
        << "Selected Baculoviridae:        "
        << selected_proteins
        << endl;

    cout
        << "Unique species:                "
        << species_set.size()
        << endl;

    cout
        << "Unique genomes:                "
        << genome_set.size()
        << endl;

    cout
        << "Unique Viro3D IDs:             "
        << viro3d_id_set.size()
        << endl;

    cout
        << "No Family mapping:             "
        << unmatched_family
        << endl;


    cout
        << endl
        << "Output CSV:"
        << endl
        << output_csv
        << endl;


    cout
        << endl
        << "Output FASTA:"
        << endl
        << output_fasta
        << endl;


    cout
        << endl
        << "Summary:"
        << endl
        << output_summary
        << endl;


    return 0;
}